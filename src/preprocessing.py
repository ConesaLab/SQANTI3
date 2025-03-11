import os
import csv
import sys
from .parsers import get_fusion_component
from .qc_classes import CAGEPeak, PolyAPeak
from .classification_utils import SJ_coverage, TSS_ratio_calculation
from .utilities.cupcake.sequence.BED import LazyBEDPointReader
from .logging_config import qc_logger

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
        qc_logger.info("**** Reading CAGE Peak data.")
        return CAGEPeak(CAGE_peak)
    return None

def read_polyA_peaks(polyA_peak):
    if polyA_peak:
        qc_logger.info("**** Reading polyA Peak data.")
        return PolyAPeak(polyA_peak)
    return None

def read_polyA_motifs(polyA_motif_list):
    if polyA_motif_list:
        qc_logger.info("**** Reading PolyA motif list.")
        motifs = []
        for line in open(polyA_motif_list):
            x = line.strip().upper().replace('U', 'A')
            if any(s not in ('A', 'T', 'C', 'G') for s in x):
                qc_logger.error(f"PolyA motif must be A/T/C/G only! Saw: {x}. Abort!")
                sys.exit(1)
            motifs.append(x)
        return motifs
    return None

def read_phyloP_bed(phyloP_bed):
    if phyloP_bed:
        qc_logger.info("**** Reading PhyloP BED file.")
        return LazyBEDPointReader(phyloP_bed)
    return None
