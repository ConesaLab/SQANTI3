import os
import csv
import sys
from .parsers import get_fusion_component
from .qc_classes import CAGEPeak, PolyAPeak
from .utilities.cupcake.sequence.BED import LazyBEDPointReader
from src.config import FIELDS_JUNC
from src.parsers import STARcov_parser
from src.utilities.short_reads import get_TSS_bed, get_bam_header, get_ratio_TSS, star
from src.utils import get_files_from_dir

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


def read_fusion_components(isoforms, is_fusion):
    if is_fusion:
        return get_fusion_component(isoforms)
    return {}

def initialize_isoform_hits(outdir, prefix, isoform_hits):
    if isoform_hits:
        isoform_hits_name = os.path.join(outdir, prefix + '_isoform_hits')
        with open(isoform_hits_name + '_tmp', 'w') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow(['Isoform', 'Isoform_length', 'Isoform_exon_number', 'Hit', 'Hit_length', 'Hit_exon_number', 'Match', 'Diff_to_TSS', 'Diff_to_TTS', 'Matching_type'])
        return isoform_hits_name
    return None

def read_CAGE_peaks(CAGE_peak):
    if CAGE_peak:
        print("**** Reading CAGE Peak data.", file=sys.stdout)
        return CAGEPeak(CAGE_peak)
    return None

def read_polyA_peaks(polyA_peak):
    if polyA_peak:
        print("**** Reading polyA Peak data.", file=sys.stdout)
        return PolyAPeak(polyA_peak)
    return None

def read_polyA_motifs(polyA_motif_list):
    if polyA_motif_list:
        print("**** Reading PolyA motif list.", file=sys.stdout)
        motifs = []
        for line in open(polyA_motif_list):
            x = line.strip().upper().replace('U', 'A')
            if any(s not in ('A', 'T', 'C', 'G') for s in x):
                print("PolyA motif must be A/T/C/G only! Saw: {0}. Abort!".format(x), file=sys.stderr)
                sys.exit(1)
            motifs.append(x)
        return motifs
    return None

def read_phyloP_bed(phyloP_bed):
    if phyloP_bed:
        print("**** Reading PhyloP BED file.", file=sys.stdout)
        return LazyBEDPointReader(phyloP_bed)
    return None
