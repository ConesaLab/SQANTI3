import os
import sys
import subprocess
import glob
from collections import defaultdict
from csv import DictReader
from bx.intervals.intersection import IntervalTree
from statistics import mean
from Bio import SeqIO
import numpy as np

from .utilities.cupcake.sequence.STAR import STARJunctionReader
from .utilities.cupcake.io.GFF import collapseGFFReader

from .config import EXP_KALLISTO_HEADERS, EXP_RSEM_HEADERS, seqid_fusion
from .qc_classes import genePredReader, myQueryProteins
from .utils import mergeDict, flatten
#from .commands import GTF2GENEPRED_PROG

def reference_parser(annot,out_dir,out_pref,genome_chroms,gene_name=False,isoAnnot=False):
    """
    Parses the reference GTF file and generates various genomic interval data structures.

    Args:
        out_dir (str): Output directory where the reference annotation file will be saved.
        out_pref (str): Prefix for the output reference annotation file.
        gene_name (bool): Flag indicating whether to use gene names.
        isoAnnot (bool): Flag indicating whether to use isoform annotations.
        annot (str): Path to the input annotation GTF file.
        genome_chroms (list): List of chromosome names from the genome fasta, used for sanity checking.

    Returns:
        tuple: A tuple containing:
            - refs_1exon_by_chr (dict): Dictionary of single exon references by chromosome.
            - refs_exons_by_chr (dict): Dictionary of multi-exon references by chromosome.
            - junctions_by_chr (dict): Dictionary of junctions by chromosome.
        - junctions_by_gene (dict): Dictionary of junctions by gene.
        - known_5_3_by_gene (dict): Dictionary of known 5' and 3' ends by gene. 
    """
    from .commands import GTF2GENEPRED_PROG

    referenceFiles = os.path.join(out_dir, "refAnnotation_"+out_pref+".genePred")
    print("**** Parsing Reference Transcriptome....", file=sys.stdout)
    if os.path.exists(referenceFiles):
        print(f"{referenceFiles} already exists. Using it.", file=sys.stdout)
    else:
        # gtf to genePred
        gtf2genepred_args = [GTF2GENEPRED_PROG, annot, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons']

        if gene_name or isoAnnot:
            gtf2genepred_args.append('-geneNameAsName2')

        subprocess.call(gtf2genepred_args)

    ## parse reference annotation
    # 1. ignore all miRNAs (< 200 bp)
    # 2. separately store single exon and multi-exon references
    refs_1exon_by_chr = defaultdict(lambda: IntervalTree()) #IntervalTree is used to efficiently hangle genomic intervals
    refs_exons_by_chr = defaultdict(lambda: IntervalTree())
    # store donors as the exon end (1-based) and acceptor as the exon start (0-based)
    # will convert the sets to sorted list later
    junctions_by_chr = defaultdict(lambda: {'donors': set(), 'acceptors': set(), 'da_pairs': set()})
    # dict of gene name --> set of junctions (don't need to record chromosome)
    junctions_by_gene = defaultdict(lambda: set())
    # dict of gene name --> list of known begins and ends (begin always < end, regardless of strand)
    known_5_3_by_gene = defaultdict(lambda: {'begin':set(), 'end': set()})

    for r in genePredReader(referenceFiles):
        known_5_3_by_gene[r.gene]['begin'].add(r.txStart)
        known_5_3_by_gene[r.gene]['end'].add(r.txEnd)
        if r.exonCount == 1:
            refs_1exon_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
        else:
            refs_exons_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
            # only store junctions for multi-exon transcripts
            for d, a in r.junctions:
                junctions_by_chr[r.chrom]['donors'].add(d)
                junctions_by_chr[r.chrom]['acceptors'].add(a)
                junctions_by_chr[r.chrom]['da_pairs'].add((d,a))
                junctions_by_gene[r.gene].add((d,a))
            

    # check that all genes' chromosomes are in the genome file
    ref_chroms = set(refs_1exon_by_chr.keys()).union(list(refs_exons_by_chr.keys()))
    diff = ref_chroms.difference(genome_chroms)
    if len(diff) > 0:
        print(f"WARNING: ref annotation contains chromosomes not in genome: {','.join(diff)}\n", file=sys.stderr)

    # convert the content of junctions_by_chr to sorted list
    # This uses dictionary to iterate over the chormoosomes, keeping the keys, but sorting the values
    junctions_by_chr = {
        k: {
            key: sorted(list(value))
            for key, value in v.items()
        }
        for k, v in junctions_by_chr.items()
}
    return dict(refs_1exon_by_chr), dict(refs_exons_by_chr), dict(junctions_by_chr), dict(junctions_by_gene), dict(known_5_3_by_gene)


def isoforms_parser(queryFile):
    """
    Parse input isoforms (GTF) to dict (chr --> sorted list)
    """
    print("**** Parsing Isoforms....", file=sys.stderr)
    
    isoforms_list = defaultdict(lambda: []) # chr --> list to be sorted later

    for r in genePredReader(queryFile):
        isoforms_list[r.chrom].append(r)

    for k in isoforms_list:
        isoforms_list[k].sort(key=lambda r: r.txStart)

    return isoforms_list


def STARcov_parser(coverageFiles): # just valid with unstrand-specific RNA-seq protocols.
    """
    :param coverageFiles: comma-separated list of STAR junction output files or a directory containing junction files
    :return: list of samples, dict of (chrom,strand) --> (0-based start, 1-based end) --> {dict of sample -> (uniq,multi) reads supporting this junction}
    """
    if os.path.isdir(coverageFiles):
        cov_files = glob.glob(coverageFiles + "/*SJ.out.tab")
    elif coverageFiles.count(',') > 0: # In case the files are as a whole string and not a list
        cov_files = coverageFiles.split(",")
    else:
        cov_files = glob.glob(coverageFiles)

    print(f"Input pattern: {coverageFiles}.\nThe following files found and to be read as junctions:\n" + '\n'.join(cov_files), file=sys.stderr)

    cov_by_chrom_strand = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: (0,0))))
    undefined_strand_count = 0
    all_read = 0
    samples = []
    for file in cov_files:
        prefix = os.path.basename(file[:file.rfind('.')]) # use this as sample name
        samples.append(prefix)
        for r in STARJunctionReader(file):
            if r.strand == 'NA':
                # undefined strand, so we put them in BOTH strands otherwise we'll lose all non-canonical junctions from STAR
                cov_by_chrom_strand[(r.chrom, '+')][(r.start, r.end)][prefix] = (r.unique_count, r.multi_count)
                cov_by_chrom_strand[(r.chrom, '-')][(r.start, r.end)][prefix] = (r.unique_count, r.multi_count)
                undefined_strand_count += 1
            else:
                cov_by_chrom_strand[(r.chrom, r.strand)][(r.start, r.end)][prefix] = (r.unique_count, r.multi_count)
            all_read += 1
    print(f"{all_read} junctions read. {undefined_strand_count} junctions added to both strands because no strand information from STAR.", file=sys.stderr)

    return samples, cov_by_chrom_strand


def expression_parser(expressionFile):
    """
    Currently accepts expression format: Kallisto or RSEM
    It need the functions mergeDict, flatten to be defined before it
    :param expressionFile: Kallisto or RSEM
    :return: dict of PBID --> TPM
    Include the possibility of providing an expression matrix --> first column must be "ID"
    """
    # Similar to the STARcov_parser, the expression file can be a directory or a list of files
    if os.path.isdir(expressionFile)==True:
                exp_paths = [os.path.join(expressionFile,fn) for fn in next(os.walk(expressionFile))[2]]
    else:
                exp_paths = expressionFile.split(",")
    exp_all = {}
    ismatrix = False
    for exp_file in exp_paths:
        reader = DictReader(open(exp_file), delimiter='\t')
        # Finds the file format based on the header
        if all(k in reader.fieldnames for k in EXP_KALLISTO_HEADERS):
                print("Detected Kallisto expression format. Using 'target_id' and 'tpm' field.", file=sys.stderr)
                name_id, name_tpm = 'target_id', 'tpm'
        elif all(k in reader.fieldnames for k in EXP_RSEM_HEADERS):
                print("Detected RSEM expression format. Using 'transcript_id' and 'TPM' field.", file=sys.stderr)
                name_id, name_tpm = 'transcript_id', 'TPM'
        elif reader.fieldnames[0]=="ID":
                print("Detected expression matrix format")
                ismatrix = True
                name_id = 'ID'
        else:
                print(f"Expected Kallisto or RSEM file format from {expressionFile}. Abort!", file=sys.stderr)
        exp_sample = {}
        if ismatrix:
            for r in reader:
                exp_sample[r[name_id]] = np.average(list(map(float,list(r.values())[1:])))
        else:
            for r in reader:
                exp_sample[r[name_id]] = float(r[name_tpm])

        exp_all = mergeDict(exp_all, exp_sample)

    exp_dict = {}
    if len(exp_paths)>1:
        for k in exp_all:
            exp_all[k] = list(flatten(exp_all[k]))
            exp_dict[k] = mean(exp_all[k])
        return exp_dict
    else:
        exp_dict=exp_all
        return exp_dict


def get_fusion_component(fusion_gtf):
    """
    Parses a GTF file containing fusion gene information and calculates the 
    cumulative exon lengths for each isoform of each gene.

    Args:
        fusion_gtf (str): Path to the GTF file containing fusion gene data.

    Returns:
        dict: A dictionary where keys are fusion gene identifiers in the format 
              "PBfusion.<gene>.<isoform>" and values are tuples containing:
              - The cumulative number of exons up to and including the current isoform.
              - The cumulative length of exons up to and including the current isoform.
    """
    components = defaultdict(lambda: {})
    for r in collapseGFFReader(fusion_gtf):
        m = seqid_fusion.match(r.seqid)
        gene, iso = int(m.group(1)), int(m.group(2))
        components[gene][iso] = sum(e.end-e.start for e in r.ref_exons)

    result = {}
    for gene, comp in components.items():
        comp = list(comp.items())
        comp.sort(key=lambda x: x[0])  # now comp is (<isoform indx>, <length>)
        _iso, _len = comp[0]
        _acc = _len
        result[f"PBfusion.{gene}.{_iso}"] = (1, _len)
        for _iso, _len in comp[1:]:
            result[f"PBfusion.{gene}.{_iso}"] = (_acc+1, _acc+_len)
            _acc += _len
    return result



# TODO: make it less specific
def FLcount_parser(fl_count_filename):
    """
    :param fl_count_filename: could be a single sample or multi-sample (chained or demux) count file
    :return: list of samples, <dict>

    If single sample, returns True, dict of {pbid} -> {count}
    If multiple sample, returns False, dict of {pbid} -> {sample} -> {count}

    For multi-sample, acceptable formats are:
    //demux-based
    id,JL3N,FL1N,CL1N,FL3N,CL3N,JL1N
    PB.2.1,0,0,1,0,0,1
    PB.3.3,33,14,47,24,15,38
    PB.3.2,2,1,0,0,0,1

    //chain-based
    superPBID<tab>sample1<tab>sample2
    """
    fl_count_dict = {}
    samples = ['NA']
    flag_single_sample = True

    f = open(fl_count_filename)
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith('#'):
            # if it first thing is superPBID or id or pbid
            if line.startswith('pbid'):
                type = 'SINGLE_SAMPLE'
                sep  = '\t'
            elif line.startswith('superPBID'):
                type = 'MULTI_CHAIN'
                sep = '\t'
            elif line.startswith('id'):
                type = 'MULTI_DEMUX'
                sep = ','
            else:
                raise Exception("Unexpected count file format! Abort!")
            f.seek(cur_pos)
            break


    reader = DictReader(f, delimiter=sep)
    count_header = reader.fieldnames
    if type=='SINGLE_SAMPLE':
        if 'count_fl' not in count_header:
            print(f"Expected `count_fl` field in count file {fl_count_filename}. Abort!", file=sys.stderr)
            sys.exit(1)
        d = dict((r['pbid'], r) for r in reader)
    elif type=='MULTI_CHAIN':
        d = dict((r['superPBID'], r) for r in reader)
        flag_single_sample = False
    elif type=='MULTI_DEMUX':
        d = dict((r['id'], r) for r in reader)
        flag_single_sample = False
    else:
        print(f"Expected pbid or superPBID as a column in count file {fl_count_filename}. Abort!", file=sys.stderr)
        sys.exit(1)
    f.close()


    if flag_single_sample: # single sample
        for k,v in d.items():
            fl_count_dict[k] = int(v['count_fl'])
    else: # multi-sample
        for k,v in d.items():
            fl_count_dict[k] = {}
            samples = list(v.keys())
            for sample,count in v.items():
                if sample not in ('superPBID', 'id'):
                    if count=='NA':
                        fl_count_dict[k][sample] = 0
                    else:
                        try:
                            fl_count_dict[k][sample] = int(count)
                        except ValueError:
                            fl_count_dict[k][sample] = float(count)

    samples.sort()

    if type=='MULTI_CHAIN':
        samples.remove('superPBID')
    elif type=='MULTI_DEMUX':
        samples.remove('id')

    return samples, fl_count_dict


def parse_corrORF(corrORF,gmst_rex):
    orfDict = {}
    for r in SeqIO.parse(open(corrORF), 'fasta'):
        # now process ORFs into myQueryProtein objects
        m = gmst_rex.match(r.description)
        if m is None:
            print(f"Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {r.description}! Abort!", file=sys.stderr)
            sys.exit(1)
        orf_length = int(m.group(2))
        cds_start = int(m.group(4))
        cds_end = int(m.group(5))
        orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, str(r.seq), proteinID=r.id)
    return orfDict

def parse_GMST(corrORF,gmst_rex,gmst_pre):
    orfDict = {}
    # Modifying ORF sequences by removing sequence before ATG
    with open(corrORF, "w") as f:
        for r in SeqIO.parse(open(gmst_pre+'.faa'), 'fasta'):
            m = gmst_rex.match(r.description)
            if m is None:
                print(f"Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {r.description}! Abort!", file=sys.stderr)
                sys.exit(1)
            id_pre = m.group(1)
            orf_length = int(m.group(2))
            orf_strand = m.group(3)
            cds_start = int(m.group(4))
            cds_end = int(m.group(5))
            pos = r.seq.find('M')
            # TODO: Revisit when the predictor has been
            if pos!= -1: # 
                # must modify both the sequence ID and the sequence
                orf_length -= pos
                cds_start += pos*3
                newid = f"{id_pre}|{orf_length}_aa|{orf_strand}|{cds_start}|{cds_end}"
                newseq = str(r.seq)[pos:]
                orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, newseq, proteinID=r.id)
                f.write(f">{newid}\n{newseq}\n")
            else:
                orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, str(r.seq), proteinID=r.id)
                f.write(f">{r.description}\n{r.seq}\n")
    return orfDict