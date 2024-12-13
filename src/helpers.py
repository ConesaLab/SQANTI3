import os
import sys
import subprocess
import gzip
import re
import bisect

from collections import defaultdict
from Bio import SeqIO
from bx.intervals import Interval

from .utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
from .utilities.cupcake.sequence.err_correct_w_genome import err_correct
from .utilities.cupcake.sequence.sam_to_gff3 import convert_sam_to_gff3

from .config import seqid_rex1, seqid_rex2, seqid_fusion
from .commands import GMAP_CMD, MINIMAP2_CMD, DESALT_CMD, ULTRA_CMD, GMST_CMD, GFFREAD_PROG
from .qc_classes import myQueryProteins

### Environment manipulation functions ###
def rename_isoform_seqids(input_fasta, force_id_ignore=False):
    """
    Rename input isoform fasta/fastq, which is usually mapped, collapsed Iso-Seq data with IDs like:

    PB.1.1|chr1:10-100|xxxxxx

    to just being "PB.1.1"

    :param input_fasta: Could be either fasta or fastq, autodetect. Can be gzipped.
    :return: output fasta with the cleaned up sequence ID, is_fusion flag
    """
    type = 'fasta'
    # gzip.open and open have different default open modes:
    # gzip.open uses "rb" (read in binary format)
    # open uses "rt" (read in text format)
    # This can be solved by making explicit the read text mode (which is required
    # by SeqIO.parse)
    open_function = gzip.open if input_fasta.endswith('.gz') else open
    with open_function(input_fasta, mode="rt") as h:
        if h.readline().startswith('@'): type = 'fastq'
    f = open(input_fasta[:input_fasta.rfind('.fast')]+'.renamed.fasta', mode='wt')
    for r in SeqIO.parse(open_function(input_fasta, "rt"), type):
        m1 = seqid_rex1.match(r.id)
        m2 = seqid_rex2.match(r.id)
        m3 = seqid_fusion.match(r.id)
        if not force_id_ignore and (m1 is None and m2 is None and m3 is None):
            print("Invalid input IDs! Expected PB.X.Y or PB.X.Y|xxxxx or PBfusion.X format but saw {0} instead. Abort!".format(r.id), file=sys.stderr)
            sys.exit(1)
        if r.id.startswith('PB.') or r.id.startswith('PBfusion.'):  # PacBio fasta header
            newid = r.id.split('|')[0]
        else:
            raw = r.id.split('|')
            if len(raw) > 4:  # RefSeq fasta header
                newid = raw[3]
            else:
                newid = r.id.split()[0]  # Ensembl fasta header
        f.write(">{0}\n{1}\n".format(newid, r.seq))
    f.close()
    return f.name


### Input/Output functions ###
def write_collapsed_GFF_with_CDS(isoforms_info, input_gff, output_gff):
    """
    Augment a collapsed GFF with CDS information
    *NEW* Also, change the "gene_id" field to use the classification result
    :param isoforms_info: dict of id -> QueryTranscript
    :param input_gff:  input GFF filename
    :param output_gff: output GFF filename
    """
    with open(output_gff, 'w') as f:
        reader = collapseGFFReader(input_gff)
        for r in reader:
            r.geneid = isoforms_info[r.seqid].geneName()  # set the gene name
            s = isoforms_info[r.seqid].CDS_genomic_start  # could be 'NA'
            e = isoforms_info[r.seqid].CDS_genomic_end    # could be 'NA'
            r.cds_exons = []
            if s!='NA' and e!='NA': # has ORF prediction for this isoform
                if r.strand == '+':
                    assert s < e
                    s = s - 1 # make it 0-based
                else:
                    assert e < s
                    s, e = e, s
                    s = s - 1 # make it 0-based
                # TODO: change the loop to a binary search (reduces complexity) 
                # TODO: Include more checks into the intervals, with an equal condition
                for i,exon in enumerate(r.ref_exons):
                    if exon.end > s: break
                r.cds_exons = [Interval(s, min(e,exon.end))]
                for exon in r.ref_exons[i+1:]:
                    if exon.start > e: break
                    r.cds_exons.append(Interval(exon.start, min(e, exon.end)))
            write_collapseGFF_format(f, r)

def get_corr_filenames(args, dir=None):
    d = dir if dir is not None else args.dir
    corrPathPrefix = os.path.join(d, args.output)
    corrGTF = corrPathPrefix +"_corrected.gtf"
    corrSAM = corrPathPrefix +"_corrected.sam"
    corrFASTA = corrPathPrefix +"_corrected.fasta"
    corrORF =  corrPathPrefix +"_corrected.faa"
    corrCDS_GTF_GFF = corrPathPrefix + "_corrected.gtf.cds.gff"
    return corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF

def get_isoform_hits_name(args, dir=None):
    d = dir if dir is not None else args.dir
    corrPathPrefix = os.path.join(d, args.output)
    isoform_hits_name = corrPathPrefix +"_isoform_hits.txt"
    return isoform_hits_name

def get_class_junc_filenames(args, dir=None):
    d = dir if dir is not None else args.dir
    outputPathPrefix = os.path.join(d, args.output)
    outputClassPath = outputPathPrefix + "_classification.txt"
    outputJuncPath = outputPathPrefix + "_junctions.txt"
    return outputClassPath, outputJuncPath

def get_omitted_name(args, dir=None):
    d = dir if dir is not None else args.dir
    corrPathPrefix = os.path.join(d, args.output)
    omitted_name = corrPathPrefix +"_omitted_due_to_min_ref_len.txt"
    return omitted_name


def write_junctionInfo(trec, junctions_by_chr, accepted_canonical_sites, indelInfo, genome_dict, fout, covInf=None, covNames=None, phyloP_reader=None):
    """
    :param trec: query isoform genePredRecord
    :param junctions_by_chr: dict of chr -> {'donors': <sorted list of donors>, 'acceptors': <sorted list of acceptors>, 'da_pairs': <sorted list of junctions>}
    :param accepted_canonical_sites: list of accepted canonical splice sites
    :param indelInfo: indels near junction information, dict of pbid --> list of junctions near indel (in Interval format)
    :param genome_dict: genome fasta dict
    :param fout: DictWriter handle
    :param covInf: (optional) junction coverage information, dict of (chrom,strand) -> (0-based start,1-based end) -> dict of {sample -> (unique, multi) read count}
    :param covNames: (optional) list of sample names for the junction coverage information
    :param phyloP_reader: (optional) dict of (chrom,0-based coord) --> phyloP score

    Write a record for each junction in query isoform
    """
    def find_closest_in_list(lst, pos):
        i = bisect.bisect_left(lst, pos)
        if i == 0:
            return lst[0]-pos
        elif i == len(lst):
            return lst[-1]-pos
        else:
            a, b = lst[i-1]-pos, lst[i]-pos
            if abs(a) < abs(b): return a
            else: return b

    # go through each trec junction
    for junction_index, (d, a) in enumerate(trec.junctions):
        # NOTE: donor just means the start, not adjusted for strand
        # Check if the chromosome of the transcript has any annotation by the reference
        # create a list in case there are chromosomes present in the input but not in the annotation dictionary junctions_by_chr
        missing_chr=[]
        junction_cat = "novel"
        if (trec.chrom in junctions_by_chr) and (trec.chrom not in missing_chr):
            # Find the closest junction start site
            min_diff_s = -find_closest_in_list(junctions_by_chr[trec.chrom]['donors'], d)
            # find the closest junction end site
            min_diff_e = find_closest_in_list(junctions_by_chr[trec.chrom]['acceptors'], a)
            if ((d,a) in junctions_by_chr[trec.chrom]['da_pairs']):
                junction_cat = "known"
        else:
            # if there is no record in the reference of junctions in this chromosome, minimum distances will be NA
            # add also new chromosome to the junctions_by_chr with one dummy SJ d=1, a=2
            if trec.chrom not in missing_chr:
                missing_chr.append(trec.chrom)
            min_diff_s = float("NaN")
            min_diff_e = float("NaN")

        splice_site = trec.get_splice_site(genome_dict, junction_index)

        indel_near_junction = "NA"
        if indelInfo is not None:
            indel_near_junction = "TRUE" if (trec.id in indelInfo and Interval(d,a) in indelInfo[trec.id]) else "FALSE"

        sample_cov = defaultdict(lambda: (0,0))  # sample -> (unique, multi) count for this junction
        if covInf is not None:
            sample_cov = covInf[(trec.chrom, trec.strand)][(d,a)]

        # if phyloP score dict exists, give the triplet score of (last base in donor exon), donor site -- similarly for acceptor
        phyloP_start, phyloP_end = 'NA', 'NA'
        if phyloP_reader is not None:
            phyloP_start = ",".join([str(x) for x in [phyloP_reader.get_pos(trec.chrom, d-1), phyloP_reader.get_pos(trec.chrom, d), phyloP_reader.get_pos(trec.chrom, d+1)]])
            phyloP_end = ",".join([str(x) for x in [phyloP_reader.get_pos(trec.chrom, a-1), phyloP_reader.get_pos(trec.chrom, a),
                                              phyloP_reader.get_pos(trec.chrom, a+1)]])

        qj = {'isoform': trec.id,
              'junction_number': "junction_"+str(junction_index+1),
              "chrom": trec.chrom,
              "strand": trec.strand,
              "genomic_start_coord": d+1,  # write out as 1-based start
              "genomic_end_coord": a,      # already is 1-based end
              "transcript_coord": "?????",  # this is where the exon ends w.r.t to id sequence, ToDo: implement later
              "junction_category": junction_cat,
              "start_site_category": "known" if min_diff_s==0 else "novel",
              "end_site_category": "known" if min_diff_e==0 else "novel",
              "diff_to_Ref_start_site": min_diff_s if min_diff_s==min_diff_s else "NA", # check if min_diff is actually nan
              "diff_to_Ref_end_site": min_diff_e if min_diff_e==min_diff_e else "NA",   # check if min_diff is actually nan
              "bite_junction": "TRUE" if ((min_diff_s<0 or min_diff_e<0) and not(min_diff_s>0 or min_diff_e>0)) else "FALSE",
              "splice_site": splice_site,
              "canonical": "canonical" if splice_site in accepted_canonical_sites else "non_canonical",
              "RTS_junction": "????", # First write ???? in _tmp, later is TRUE/FALSE
              "indel_near_junct": indel_near_junction,
              "phyloP_start": phyloP_start,
              "phyloP_end": phyloP_end,
              "sample_with_cov": sum([cov_uniq>0 for (cov_uniq,cov_multi) in sample_cov.values()]) if covInf is not None else "NA",
              "total_coverage_unique": sum([cov_uniq for (cov_uniq,cov_multi ) in sample_cov.values()]) if covInf is not None else "NA",
              "total_coverage_multi": sum([cov_multi for (cov_uniq,cov_multi ) in sample_cov.values()]) if covInf is not None else "NA"}

        if covInf is not None:
            for sample in covNames:
                cov_uniq, cov_multi = sample_cov[sample]
                qj[sample+'_unique'] = str(cov_uniq)
                qj[sample+'_multi'] = str(cov_multi)

        fout.writerow(qj)


### CMD related functions ###

def run_command(cmd, description="command execution"):
    """
    Executes a shell command and handles errors gracefully.
    
    :param cmd: The command to execute (string).
    :param description: A short description of the operation for better error messages (default: "command execution").
    :raises SystemExit: Exits the script if the command fails.
    """
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR during {description}: {cmd}", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)


### General functions ###
def correctionPlusORFpred(args, genome_dict, badstrandGTF):
    """
    Use the reference genome to correct the sequences (unless a pre-corrected GTF is given)
    """
    global corrORF
    global corrGTF
    global corrSAM
    global corrFASTA

    corrGTF, corrSAM, corrFASTA, corrORF , corrCDS_GTF_GFF = get_corr_filenames(args)
    p = os.path.splitext(os.path.basename(corrSAM))[0]
    n_cpu = max(1, args.cpus // args.chunks)

    # Step 1. IF GFF or GTF is provided, make it into a genome-based fasta
    #         IF sequence is provided, align as SAM then correct with genome
    if os.path.exists(corrFASTA):
        print("Error corrected FASTA {0} already exists. Using it...".format(corrFASTA), file=sys.stderr)
    else:
        if args.fasta:
            if os.path.exists(corrSAM):
                print("Aligned SAM {0} already exists. Using it...".format(corrSAM), file=sys.stderr)
            else:
                # Even though the speed does not change form the ifelse, this is cleaner
                match args.aligner_choice:
                    case "gmap":
                        print("****Aligning reads with GMAP...", file=sys.stdout)
                        cmd = GMAP_CMD.format(
                            cpus=n_cpu,
                            dir=os.path.dirname(args.gmap_index),
                            name=os.path.basename(args.gmap_index),
                            sense=args.sense,
                            i=args.isoforms,
                            o=corrSAM,
                        )
                    case "minimap2":
                        print("****Aligning reads with Minimap2...", file=sys.stdout)
                        cmd = MINIMAP2_CMD.format(
                            cpus=n_cpu,
                            sense=args.sense,
                            g=args.genome,
                            i=args.isoforms,
                            o=corrSAM,
                        )
                    case "deSALT":
                        print("****Aligning reads with deSALT...", file=sys.stdout)
                        cmd = DESALT_CMD.format(
                            cpus=n_cpu,
                            dir=args.gmap_index,
                            i=args.isoforms,
                            o=corrSAM,
                        )
                    case "uLTRA":
                        print("****Aligning reads with uLTRA...", file=sys.stdout)
                        cmd = ULTRA_CMD.format(
                            cpus=n_cpu,
                            prefix="../" + p,
                            g=args.genome,
                            a=args.annotation,
                            i=args.isoforms,
                            o_dir=args.dir + "/uLTRA_out/",
                        )
                    case _:
                        raise ValueError(f"Unsupported aligner choice: {args.aligner_choice}")

                run_command(cmd, description="aligning reads")

            # error correct the genome (input: corrSAM, output: corrFASTA)
            err_correct(args.genome, corrSAM, corrFASTA, genome_dict=genome_dict)
            # convert SAM to GFF --> GTF
            convert_sam_to_gff3(corrSAM, corrGTF+'.tmp', source=os.path.basename(args.genome).split('.')[0])  # convert SAM to GFF3
            cmd = "{p} {o}.tmp -T -o {o}".format(o=corrGTF, p=GFFREAD_PROG)
            # Try condition to better handle the error. Also, the exit code is corrected
            run_command(cmd, description="converting SAM to GTF")
        else:
            print("Skipping aligning of sequences because GTF file was provided.", file=sys.stdout)

            ind = 0
            with open(args.isoforms) as isoforms_gtf:
                for line in isoforms_gtf:
                    if line[0] != "#" and len(line.split("\t"))!=9:
                        sys.stderr.write("\nERROR: input isoforms file with not GTF format.\n")
                        sys.exit()
                    elif len(line.split("\t"))==9:
                        ind += 1
                if ind == 0:
                    print("WARNING: GTF has {0} no annotation lines.".format(args.isoforms), file=sys.stderr)


            # GFF to GTF (in case the user provides gff instead of gtf)
            corrGTF_tpm = corrGTF+".tmp"
            # Use the run command?
            try:
                subprocess.call([GFFREAD_PROG, args.isoforms , '-T', '-o', corrGTF_tpm])
            except (RuntimeError, TypeError, NameError):
                sys.stderr.write('ERROR: File %s without GTF/GFF format.\n' % args.isoforms)
                raise SystemExit(1)


            # check if gtf chromosomes inside genome file
            with open(corrGTF, 'w') as corrGTF_out:
                with open(corrGTF_tpm, 'r') as isoforms_gtf:
                    with (open(badstrandGTF, 'w')) as discard_gtf:
                        for line in isoforms_gtf:
                            if line[0] != "#":
                                chrom = line.split("\t")[0]
                                type = line.split("\t")[2]
                                strand = line.split("\t")[6]
                                if chrom not in list(genome_dict.keys()):
                                    sys.stderr.write("\nERROR: gtf \"%s\" chromosome not found in genome reference file.\n" % (chrom))
                                    sys.exit()
                                elif type in ('transcript', 'exon'):
                                    # In normal cirumstances, strand should be a string of values '-' or '+'
                                    # However, the strand can also be '.' , a dot, which means that
                                    # The strand is unknown. The frequence of these varies according to technologies 
                                    # and not taken into account downstream
                                    if(strand not in ['-','+']):
                                        print("WARNING: Discarding uknown strand transcript ")
                                        discard_gtf.write(line)
                                        continue
                                    corrGTF_out.write(line)
            os.remove(corrGTF_tpm)

            if not os.path.exists(corrSAM):
                sys.stdout.write("\nIndels will be not calculated since you ran SQANTI3 without alignment step (SQANTI3 with gtf format as transcriptome input).\n")

            # GTF to FASTA
            subprocess.call([GFFREAD_PROG, corrGTF, '-g', args.genome, '-w', corrFASTA])

    # ORF generation
    print("**** Predicting ORF sequences...", file=sys.stdout)

    gmst_dir = os.path.join(os.path.abspath(args.dir), "GMST")
    gmst_pre = os.path.join(gmst_dir, "GMST_tmp")
    if not os.path.exists(gmst_dir):
        os.makedirs(gmst_dir) 

    # sequence ID example: PB.2.1 gene_4|GeneMark.hmm|264_aa|+|888|1682
    gmst_rex = re.compile(r'(\S+\t\S+\|GeneMark.hmm)\|(\d+)_aa\|(\S)\|(\d+)\|(\d+)')
    orfDict = {}  # GMST seq id --> myQueryProteins object
    if args.skipORF:
        print("WARNING: Skipping ORF prediction because user requested it. All isoforms will be non-coding!", file=sys.stderr)
    elif os.path.exists(corrORF):
        print("ORF file {0} already exists. Using it....".format(corrORF), file=sys.stderr)
        for r in SeqIO.parse(open(corrORF), 'fasta'):
            # now process ORFs into myQueryProtein objects
            m = gmst_rex.match(r.description)
            if m is None:
                print("Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {0}! Abort!".format(r.description), file=sys.stderr)
                sys.exit(1)
            orf_length = int(m.group(2))
            cds_start = int(m.group(4))
            cds_end = int(m.group(5))
            orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, str(r.seq), proteinID=r.id)
    else:
        cur_dir = os.path.abspath(os.getcwd())
        os.chdir(args.dir)
        if args.orf_input is not None:
            print("Running ORF prediction of input on {0}...".format(args.orf_input))
            cmd = GMST_CMD.format(i=os.path.realpath(args.orf_input), o=gmst_pre)
        else:
            cmd = GMST_CMD.format(i=corrFASTA, o=gmst_pre)
        run_command(cmd, description="GMST ORF prediction")
        os.chdir(cur_dir)
        # Modifying ORF sequences by removing sequence before ATG
        with open(corrORF, "w") as f:
            for r in SeqIO.parse(open(gmst_pre+'.faa'), 'fasta'):
                m = gmst_rex.match(r.description)
                if m is None:
                    print("Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {0}! Abort!".format(r.description), file=sys.stderr)
                    sys.exit(1)
                id_pre = m.group(1)
                orf_length = int(m.group(2))
                orf_strand = m.group(3)
                cds_start = int(m.group(4))
                cds_end = int(m.group(5))
                pos = r.seq.find('M')
                if pos!=-1:
                    # must modify both the sequence ID and the sequence
                    orf_length -= pos
                    cds_start += pos*3
                    newid = "{0}|{1}_aa|{2}|{3}|{4}".format(id_pre, orf_length, orf_strand, cds_start, cds_end)
                    newseq = str(r.seq)[pos:]
                    orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, newseq, proteinID=newid)
                    f.write(">{0}\n{1}\n".format(newid, newseq))
                else:
                    new_rec = r
                    orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, str(r.seq), proteinID=r.id)
                    f.write(">{0}\n{1}\n".format(new_rec.description, new_rec.seq))

    if len(orfDict) == 0:
        print("WARNING: All input isoforms were predicted as non-coding", file=sys.stderr)

    return(orfDict)
