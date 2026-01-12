from collections import defaultdict
from Bio  import SeqIO

from src.utilities.rt_switching import rts
from src.parsers import FLcount_parser, expression_parser
from src.utilities.short_reads import kallisto
from src.utils import pstdev
from src.module_logging import qc_logger

def process_rts(isoforms_info, outputJuncPath, genome, genome_dict=None,extension="_tmp"):
    """
    RTS_info: dict of (pbid) -> list of RT junction. if RTS_info[pbid] == [], means all junctions are non-RT.
    """
    if genome_dict is None:
        genome_dict = dict((r.name, r) for r in SeqIO.parse(open(genome), 'fasta'))

    RTS_info = rts([outputJuncPath+extension, genome, "-a"], genome_dict)
    for pbid in isoforms_info:
        if pbid in RTS_info:
            isoforms_info[pbid].RTS_stage = "TRUE"
        else:
            isoforms_info[pbid].RTS_stage = "FALSE"

    return isoforms_info, RTS_info

def classify_fsm(isoforms_info):
    """Classify Full Splice Match (FSM) for each isoform."""
    geneFSM_dict = defaultdict(list)
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()
        geneFSM_dict[gene].append(isoforms_info[iso].structural_category)
    
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()
        if len(geneFSM_dict[gene]) == 1:
            isoforms_info[iso].FSM_class = "A"
        elif "full-splice_match" in geneFSM_dict[gene]:
            isoforms_info[iso].FSM_class = "C"
        else:
            isoforms_info[iso].FSM_class = "B"
    return isoforms_info

def ratio_TSS_dict_reading(isoforms_info,ratio_TSS_dict):
    qc_logger.info('**** Adding TSS ratio data.')
    for iso in ratio_TSS_dict:
        if iso not in isoforms_info:
            qc_logger.warning(f"Isoform {iso} found in ratio TSS file but not in input FASTA/GTF")
    for iso in isoforms_info:
        if iso in ratio_TSS_dict:
            if str(ratio_TSS_dict[iso]['return_ratio']) == 'nan':
                isoforms_info[iso].ratio_TSS = None
            else:
                isoforms_info[iso].ratio_TSS = ratio_TSS_dict[iso]['return_ratio']
        else:
            qc_logger.warning(f"Isoform {iso} not found in ratio TSS file. Assign count as 1.")
            isoforms_info[iso].ratio_TSS = 1
    return isoforms_info

def full_length_quantification(fl_count, isoforms_info,fields_class_cur):

    qc_logger.info("**** Reading Full-length read abundance files.")
    fl_samples, fl_count_dict = FLcount_parser(fl_count)
    n = 0
    if len(fl_samples) == 1: # single sample from PacBio
        qc_logger.info("Single-sample PacBio FL count format detected.")
        for iso, obj in isoforms_info.items():
            count = fl_count_dict.get(iso)
            if count is not None:
                obj.FL = count
            else:
                n += 1
                obj.FL = 0
    else: # multi-sample
        qc_logger.info("Multi-sample PacBio FL count format detected.")
        fields_class_cur.extend(f"FL.{s}" for s in fl_samples)
        for iso, obj in isoforms_info.items():
            count = fl_count_dict.get(iso)
            if count is not None:
                obj.FL_dict = count
            else:
                n += 1
                obj.FL_dict = defaultdict(int)
    if n > 0:
        qc_logger.warning(f"{n} isoforms not found in FL count file. Assigned counts as 0.")
    return isoforms_info, fields_class_cur

def isoform_expression_info(isoforms_info,expression,short_reads,outdir,corrFASTA,cpus):
    if expression:
        qc_logger.info("**** Reading Isoform Expression Information.")
        exp_dict = expression_parser(expression)
        gene_exp_dict = {}
        for iso in isoforms_info:
            if iso not in exp_dict:
                exp_dict[iso] = 0
                qc_logger.warning(f"Isoform {iso} not found in expression matrix. Assigning TPM of 0.")
            gene = isoforms_info[iso].geneName()
            if gene not in gene_exp_dict:
                gene_exp_dict[gene] = exp_dict[iso]
            else:
                gene_exp_dict[gene] = gene_exp_dict[gene]+exp_dict[iso]
    else:
        if short_reads is not None:
            qc_logger.info("**** Running Kallisto to calculate isoform expressions. ")
            expression_files = kallisto(corrFASTA, short_reads, outdir, cpus)
            exp_dict = expression_parser(expression_files)
            gene_exp_dict = {}
            for iso in isoforms_info:
                if iso not in exp_dict:
                    exp_dict[iso] = 0
                    qc_logger.warning(f"Isoform {iso} not found in expression matrix. Assigning TPM of 0.")
                gene = isoforms_info[iso].geneName()
                if gene not in gene_exp_dict:
                    gene_exp_dict[gene] = exp_dict[iso]
                else:
                    gene_exp_dict[gene] = gene_exp_dict[gene]+exp_dict[iso]
        else:
            exp_dict = None
            gene_exp_dict = None
            qc_logger.info("Isoforms expression files not provided.")
    # Add expression information to isoforms_info
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()
        if exp_dict is not None and gene_exp_dict is not None:
            isoforms_info[iso].gene_exp = gene_exp_dict[gene]
            isoforms_info[iso].iso_exp  = exp_dict[iso]
    return isoforms_info


def isoforms_junctions(isoforms_info, reader):
    # Read the junction information to fill in several remaining unfilled fields in classification
    # (1) "canonical": is "canonical" if all junctions are canonical, otherwise "non_canonical"
    # (2) "bite": is TRUE if any of the junction "bite_junction" field is TRUE

    sj_covs_by_isoform = defaultdict(lambda: [])  # pbid --> list of total_cov for each junction so we can calculate SD later
    for r in reader:
        # only need to do assignment if:
        # (1) the .canonical field is still "NA"
        # (2) the junction is non-canonical
        assert r['canonical'] in ('canonical', 'non_canonical')
        if (isoforms_info[r['isoform']].all_canonical is None) or \
            (r['canonical'] == 'non_canonical'):
            isoforms_info[r['isoform']].all_canonical = r['canonical']

        if (isoforms_info[r['isoform']].bite is None) or (r['bite_junction'] == 'TRUE'):
            isoforms_info[r['isoform']].bite = r['bite_junction']

        if r['indel_near_junct'] == 'TRUE':
            if isoforms_info[r['isoform']].n_indels_junc is None:
                isoforms_info[r['isoform']].n_indels_junc = 0
            isoforms_info[r['isoform']].n_indels_junc += 1

        # min_cov: min( total_cov[j] for each junction j in this isoform )
        # min_cov_pos: the junction [j] that attributed to argmin(total_cov[j])
        # min_sample_cov: min( sample_cov[j] for each junction in this isoform )
        # sd_cov: sd( total_cov[j] for each junction j in this isoform )
        if r['sample_with_cov'] != 'NA':
            sample_with_cov = int(r['sample_with_cov'])
            if (isoforms_info[r['isoform']].min_sample_cov is None) or (isoforms_info[r['isoform']].min_sample_cov > sample_with_cov):
                isoforms_info[r['isoform']].min_sample_cov = sample_with_cov

        if r['total_coverage_unique'] != 'NA':
            total_cov = int(r['total_coverage_unique'])
            sj_covs_by_isoform[r['isoform']].append(total_cov)
            if (isoforms_info[r['isoform']].min_cov is None) or (isoforms_info[r['isoform']].min_cov > total_cov):
                isoforms_info[r['isoform']].min_cov = total_cov
                isoforms_info[r['isoform']].min_cov_pos = r['junction_number']


    for pbid, covs in sj_covs_by_isoform.items():
        isoforms_info[pbid].sd = pstdev(covs)
    return isoforms_info