import os
import csv
import sys

from src.utilities.cupcake.sequence.BED import LazyBEDPointReader

from src.parsers import get_fusion_component
from src.qc_classes import CAGEPeak, PolyAPeak

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
