import os
import sys
import glob
import re

from collections import defaultdict
from csv import DictReader
from bx.intervals.intersection import IntervalTree
from statistics import mean
from Bio import SeqIO
import numpy as np

from src.utilities.cupcake.sequence.STAR import STARJunctionReader
from src.utilities.cupcake.io.GFF import collapseGFFReader

from src.config import EXP_KALLISTO_HEADERS, EXP_RSEM_HEADERS, seqid_fusion
from src.qc_classes import genePredReader, myQueryProteins
from src.utils import mergeDict, flatten
from src.module_logging import qc_logger
from src.commands import run_command
#from src.commands import GTF2GENEPRED_PROG

def reference_parser(annot,out_dir,out_pref,genome_chroms,gene_name=False,isoAnnot=False,logger=qc_logger):
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
    from src.commands import GTF2GENEPRED_PROG

    referenceFiles = os.path.join(out_dir, "refAnnotation_"+out_pref+".genePred")
    logger.info("**** Parsing Reference Transcriptome....")
    if os.path.exists(referenceFiles):
        logger.info(f"{referenceFiles} already exists. Using it.")
    else:
        # gtf to genePred
        cmd = f"{GTF2GENEPRED_PROG} {annot} {referenceFiles} -genePredExt -allErrors -ignoreGroupsWithoutExons"
        if gene_name or isoAnnot: #TODO: Discover why this flag was here or isoAnnot:
            if isoAnnot:
                qc_logger.warning("IsoAnnotLite needs the reference annotation to have the value 'gene_name' in the attributes column")
                qc_logger.warning("If the QC report has no genes, double check your reference GTF.")
            cmd += ' -geneNameAsName2'
        run_command(cmd,logger, f"{out_dir}/logs/GTF_to_genePred.log", "GTF to genePred conversion")

    ## parse reference annotation
    # 1. ignore all miRNAs (< 200 bp)
    # 2. separately store single exon and multi-exon references
    refs_1exon_by_chr = defaultdict(lambda: IntervalTree()) #IntervalTree is used to efficiently hangle genomic intervals
    refs_exons_by_chr = defaultdict(lambda: IntervalTree())
    # store donors as the exon end (1-based) and acceptor as the exon start (0-based)
    # will convert the sets to sorted list later
    junctions_by_chr = defaultdict(lambda: {'donors': set(), 
                                            'acceptors': set(),
                                            'da_pairs': {'+': set(), '-': set()}})
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
                junctions_by_chr[r.chrom]['da_pairs'][r.strand].add((d,a))
                junctions_by_gene[r.gene].add((d,a))

    # check that all genes' chromosomes are in the genome file
    ref_chroms = set(refs_1exon_by_chr.keys()).union(list(refs_exons_by_chr.keys()))
    diff = ref_chroms.difference(genome_chroms)
    if len(diff) > 0:
        logger.warning(f"Reference annotation contains chromosomes not in genome: {','.join(diff)}\n")

    # convert the content of junctions_by_chr to sorted list
    # This uses dictionary to iterate over the chromosomes, keeping the keys, but sorting the values
    junctions_by_chr = {
        k: {
            key: sorted(list(value)) if key != 'da_pairs' else {
                strand: sorted(list(pairs)) for strand, pairs in value.items()
            }
            for key, value in v.items()
        }
        for k, v in junctions_by_chr.items()
    }
    
    # TODO: Find a more efficient way to fix this
    for chr in junctions_by_chr.keys():
        junctions_by_chr[chr]['da_tree'] = IntervalTree()
        for strand in junctions_by_chr[chr]['da_pairs'].keys():
            for junction in junctions_by_chr[chr]['da_pairs'][strand]:
                junctions_by_chr[chr]['da_tree'].insert(junction[0], junction[1], (*junction,strand))
    return dict(refs_1exon_by_chr), dict(refs_exons_by_chr), dict(junctions_by_chr), dict(junctions_by_gene), dict(known_5_3_by_gene)


def isoforms_parser(queryFile):
    """
    Parse input isoforms (GTF) to dict (chr --> sorted list)
    """
    qc_logger.info("**** Parsing Isoforms.")
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

    qc_logger.info(f"Input pattern: {coverageFiles}.\nThe following files found and to be read as junctions:\n" + '\n'.join(cov_files))

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
    qc_logger.info(f"{all_read} junctions read. {undefined_strand_count} junctions added to both strands because no strand information from STAR.")

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
                qc_logger.info("Detected Kallisto expression format. Using 'target_id' and 'tpm' field.")
                name_id, name_tpm = 'target_id', 'tpm'
        elif all(k in reader.fieldnames for k in EXP_RSEM_HEADERS):
                qc_logger.info("Detected RSEM expression format. Using 'transcript_id' and 'TPM' field.")
                name_id, name_tpm = 'transcript_id', 'TPM'
        elif reader.fieldnames[0]=="ID":
                qc_logger.info("Detected expression matrix format")
                ismatrix = True
                name_id = 'ID'
        else:
                qc_logger.error(f"Expected Kallisto or RSEM file format from {expressionFile}. Abort!")
                sys.exit(1)
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
        # Fix negative strand coordinates
        if r.strand == '-': 
            r.start, r.end =  r.ref_exons[-1].start, r.ref_exons[0].end # Positions get shifted by 1
            r.ref_exons.reverse() # Fix the reference exons 
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
            qc_logger.error(f"Expected `count_fl` field in count file {fl_count_filename}. Abort!")
            sys.exit(1)
        d = dict((r['pbid'], r) for r in reader)
    elif type=='MULTI_CHAIN':
        d = dict((r['superPBID'], r) for r in reader)
        flag_single_sample = False
    elif type=='MULTI_DEMUX':
        d = dict((r['id'], r) for r in reader)
        flag_single_sample = False
    else:
        qc_logger.error(f"Expected pbid or superPBID as a column in count file {fl_count_filename}. Abort!")
        sys.exit(1)
    f.close()


    if flag_single_sample: # single sample
        for k,v in d.items():
            try:
                fl_count_dict[k] = int(v['count_fl'])
            except ValueError:
                fl_count_dict[k] = float(v['count_fl'])
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


def parse_corrORF(corrORF):
    orfDict = {}
    for r in SeqIO.parse(open(corrORF), 'fasta'):
        # now process ORFs into myQueryProtein objects
        pattern = re.compile(r'^\S+\s+(\S+)\|(\d+)_aa\|([+-])\|(\d+)\|(\d+)$')
        m = pattern.match(r.description)
        if m is None:
            qc_logger.error(f"Expected the CDS IDs to be of format '<pbid> cds_name|<size>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {r.description}! Abort!")
            sys.exit(1)
        orf_length = int(m.group(2))
        cds_start = int(m.group(4))
        cds_end = int(m.group(5))
        orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, str(r.seq), proteinID=r.id)
    return orfDict

def parse_td2_to_dict(td2_faa):
    """
    Parses the TD2 FASTA file to extract ORF information.

    Returns:
        orfDict (dict): Keys are ORF IDs, values are myQueryProteins objects.
        records (list): List of dicts with keys: record, id_pre, orf_length, orf_strand, cds_start, cds_end
    """
    orfDict = {}
    records = []

    for r in SeqIO.parse(open(td2_faa), 'fasta'):
        info = extract_variables(r.description)
        id_pre = info['id_pre']

        orfDict[id_pre] = myQueryProteins(
            info['cds_start'],
            info['cds_end'],
            info['orf_length'],
            str(r.seq),
            proteinID=id_pre
        )

        # Include the record object along with extracted info
        records.append({
            'record': r,
            **info
        })

    return orfDict, records

def write_corr_orf(corrORF, records):
    """
    Writes reformatted ORF information to a file.
    """
    with open(corrORF, "w") as f:
        for entry in records:
            r = entry['record']
            newid = f"{entry['id_pre']}\t{r.id}|{entry['orf_length']}_aa|{entry['orf_strand']}|{entry['cds_start']}|{entry['cds_end']}"
            f.write(f">{newid}\n{str(r.seq)}\n")


def parse_TD2(corrORF, td2_faa):
    """
    Parses the TD2 output and writes corrected FASTA entries to a file.

    Returns:
        orfDict (dict): Dictionary with ORF metadata.
    """
    orfDict, records = parse_td2_to_dict(td2_faa)
    write_corr_orf(corrORF, records)
    return orfDict
    
def extract_variables(s):
    # Extract the ID prefix and CDS coordinates with strand
    match = re.search(r'([^\s:]+):(\d+)-(\d+)\(([\+\-])\)', s)
    if not match:
        qc_logger.error(f"Failed to parse coordinates and strand from: {s}")
        sys.exit(1)

    id_pre = match.group(1)
    cds_start = int(match.group(2))
    cds_end = int(match.group(3))  # Add 1 as requested
    orf_strand = match.group(4)
    orf_length = abs(cds_end - cds_start +1) // 3  # Calculate ORF length in amino acids

    return {
        'id_pre': id_pre,
        'cds_start': cds_start,
        'cds_end': cds_end,
        'orf_strand': orf_strand,
        'orf_length': orf_length
        }   
    