import os
import sys
from bx.intervals import IntervalTree

from src.config import FIELDS_JUNC
from src.parsers import STARcov_parser
from src.utilities.short_reads import get_TSS_bed, get_bam_header, get_ratio_TSS, star
from src.utils import get_files_from_dir

def calc_overlap(s1, e1, s2, e2):
    """returns the overlap between two intervals"""
    if s1=='NA' or s2=='NA': return 0
    if s1 > s2:
        s1, e1, s2, e2 = s2, e2, s1, e1
    return max(0, min(e1,e2)-max(s1,s2) + 1) # Wouldn't this need a +1?

def calc_splicesite_agreement(query_exons, ref_exons):
    """Determines the number of splice sites that agree between query and reference"""
    q_sites = set(site for exon in query_exons for site in (exon.start, exon.end))
    return sum(1 for exon in ref_exons for site in (exon.start, exon.end) if site in q_sites)

def calc_exon_overlap(query_exons, ref_exons):
    """Determines the number of bases that overlap between query and reference"""
    query_ranges = set(b for e in query_exons for b in range(e.start, e.end))
    return sum(1 for e in ref_exons for b in range(e.start, e.end) if b in query_ranges)

def get_diff_tss_tts(trec, ref):
    """
    Calculate the differences between the Transcript Start Site (TSS) and 
    Transcript Termination Site (TTS) of two transcripts.

    For positive strand transcripts:
    - TSS is the difference between the reference's start and the transcript's start.
    - TTS is the difference between the transcript's end and the reference's end.
    Positive values indicate elongation, negative values indicate shortening.

    For negative strand transcripts:
    - The start and end are reversed (i.e., `start` corresponds to the `end` in positive strand, and vice versa).
    - TTS is calculated as the difference between the reference's start and the transcript's start.
    - TSS is the difference between the transcript's end and the reference's end.
    Again, positive values indicate elongation, negative values indicate shortening.

    Parameters:
    - trec: The transcript object representing the query transcript.
    - ref: The reference transcript object.

    Returns:
    - diff_tss (int): Difference in the transcript start sites (positive for elongation, negative for shortening).
    - diff_tts (int): Difference in the transcript termination sites (positive for elongation, negative for shortening).
    """
    if trec.strand == '+':
        diff_tss = ref.txStart - trec.txStart
        diff_tts = trec.txEnd  - ref.txEnd
    else:
        diff_tts = ref.txStart - trec.txStart
        diff_tss = trec.txEnd  - ref.txEnd
    return diff_tss, diff_tts


def get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene):
    """
    Calculate the nearest transcription start site (TSS) and transcription termination site (TTS) 
    differences for a given isoform hit relative to all isoforms of the gene.

    Args:
        isoform_hit: An object representing the isoform hit, which contains information about 
                        the genes it hits and their start/end sites.

    Modifies:
        isoform_hit: Updates the `tss_gene_diff` and `tts_gene_diff` attributes with the nearest 
                        TSS and TTS differences, respectively. If no valid difference is found, 
                        the attribute is set to 'NA'.
    """
    # now that we know the reference (isoform) it hits
    # add the nearest start/end site for that gene (all isoforms of the gene)
    nearest_start_diff, nearest_end_diff = float('inf'), float('inf')
    for ref_gene in isoform_hit.genes:
        for x in start_ends_by_gene[ref_gene]['begin']:
            d =  x - trec.txStart
            if abs(d) < abs(nearest_start_diff):
                nearest_start_diff = d
        for x in start_ends_by_gene[ref_gene]['end']:
            d = trec.txEnd - x
            if abs(d) < abs(nearest_end_diff):
                nearest_end_diff = d

    if trec.strand == '+':
        isoform_hit.tss_gene_diff = nearest_start_diff if nearest_start_diff!=float('inf') else 'NA'
        isoform_hit.tts_gene_diff = nearest_end_diff if nearest_end_diff!=float('inf') else 'NA'
    else:
        isoform_hit.tss_gene_diff = nearest_end_diff if nearest_start_diff!=float('inf') else 'NA'
        isoform_hit.tts_gene_diff = nearest_start_diff if nearest_end_diff!=float('inf') else 'NA'

def categorize_incomplete_matches(trec, ref):
    """
    intron_retention --- at least one trec exon covers at least two adjacent ref exons
    complete --- all junctions agree and is not IR
    5prime_fragment --- all junctions agree but trec has less 5' exons. The isoform is a 5' fragment of the reference transcript
    3prime_fragment --- all junctions agree but trec has less 3' exons. The isoform is a 3' fragment of the reference transcript
    internal_fragment --- all junctions agree but trec has less 5' and 3' exons
    """
    # check intron retention
    ref_exon_tree = IntervalTree()
    for i,e in enumerate(ref.exons): ref_exon_tree.insert(e.start, e.end, i)
    for e in trec.exons:
        if len(ref_exon_tree.find(e.start, e.end)) > 1: # multiple ref exons covered
            return "intron_retention"

    agree_front = trec.junctions[0]==ref.junctions[0]
    agree_end   = trec.junctions[-1]==ref.junctions[-1]
    if agree_front:
        if agree_end:
            return "complete"
        else: # front agrees, end does not
            return ("5prime_fragment" if trec.strand=='+' else '3prime_fragment')
    else:
        if agree_end: # front does not agree, end agrees
            return ("3prime_fragment" if trec.strand=='+' else '5prime_fragment')
        else:
            return "internal_fragment"


def SJ_coverage(short_reads,coverage_file,genome,outdir,cpus):
    """
    Calculate short-read coverage using STAR and parse the coverage information.
    Parameters:
    short_reads (str): Path to the short-read sequencing data file.
    coverage_file (str): Path to the precomputed coverage file. If provided, short-read mapping is skipped.
    genome (str): Path to the genome reference file.
    outdir (str): Directory where output files will be saved.
    cpus (int): Number of CPUs to use for STAR mapping.
    Returns:
    tuple: A tuple containing:
        - mapping_out (str): Path to the mapping output file.
        - star_index (str or None): Path to the STAR index directory, or None if coverage_file is provided.
        - SJcovNames (list): List of splice junction coverage names.
        - SJcovInfo (dict): Dictionary containing splice junction coverage information.
        - fields_junc_cur (list): List of fields for the junctions, including sample-specific coverage fields.
    """
    fields_junc_cur = FIELDS_JUNC # add the samples to the header
    mapping_out, star_index, SJcovInfo, SJcovNames = None, None, None, None
    if coverage_file is None:
            ## Short-read mapping
        if short_reads is not None:
            print("**** Running STAR for calculating Short-Read Coverage.", file=sys.stdout)
            mapping_out, star_index = star(genome,short_reads,outdir,cpus)
        else:
            print("No short-reads or coverage provided. Skipping short-read coverage calculation.", file=sys.stdout)
            return mapping_out, star_index, SJcovInfo, SJcovNames 
    else:
        mapping_out = coverage_file
        star_index = None
    
    SJcovNames, SJcovInfo = STARcov_parser(mapping_out)
    for name in SJcovNames:
            fields_junc_cur += [name + '_unique', name + '_multi']

    return mapping_out, star_index, SJcovNames, SJcovInfo, fields_junc_cur


def TSS_ratio_calculation(SR_bam,short_reads,star_out,star_index,corrGTF,ratio_TSS_metric):
    bams = None
    ratio_TSS_dict = None
    ## TSS ratio calculation
    if  SR_bam is not None:
        print("Using provided BAM files for calculating TSS ratio", file=sys.stdout)
        bams = get_files_from_dir(SR_bam,".bam")

        # TODO: where did this functions come from?
        chr_order = get_bam_header(bams[0])
        inside_bed, outside_bed = get_TSS_bed(corrGTF, chr_order)
        ratio_TSS_dict = get_ratio_TSS(inside_bed, outside_bed, bams, chr_order, ratio_TSS_metric)
    else:
        if short_reads is not None: # If short reads are provided, it looks for the STAR output
            print("Running calculation of TSS ratio", file=sys.stdout)
            chr_order = star_index + "/chrNameLength.txt"
            inside_bed, outside_bed = get_TSS_bed(corrGTF, chr_order)
            bams = get_files_from_dir(star_out,"SJ.out.tab")
            ratio_TSS_dict = get_ratio_TSS(inside_bed, outside_bed, bams, chr_order, ratio_TSS_metric)
        else:
            print('**** TSS ratio will not be calculated since SR information was not provided')
    return ratio_TSS_dict