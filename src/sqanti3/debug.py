import argparse
import bisect
import copy
import distutils.spawn
import glob
import itertools
import math
import os
import re
import shutil
import subprocess
import sys
import timeit
from collections import Counter, defaultdict, namedtuple
from collections.abc import Iterable
from csv import DictReader, DictWriter
from multiprocessing import Process

import numpy as np
from pygmst.pygmst import gmst
from Bio import SeqIO
from bx.intervals import Interval, IntervalTree
from cupcake.cupcake.tofu.compare_junctions import compare_junctions
from cupcake.sequence.BED import LazyBEDPointReader
from cupcake.sequence.err_correct_w_genome import err_correct
from cupcake.sequence.GFF import collapseGFFReader, write_collapseGFF_format
from cupcake.sequence.sam_to_gff3 import convert_sam_to_gff3
from cupcake.sequence.STAR import STARJunctionReader
from sqanti3.utilities.indels_annot import calc_indels_from_sam
from sqanti3.utilities.rt_switching import rts
from sqanti3.sqanti3_qc import (
    outputClassPath,
    outputJuncPath,
    FIELDS_CLASS,
    transcriptsKnownSpliceSites,
    novelIsoformsKnownGenes,
    associationOverlapping,
    write_junctionInfo,
    find_polyA_motif,
)


def isoformClassification(
    coverage,
    cage_peak,
    polyA_peak,
    polyA_motif_list,
    phyloP_bed,
    sites,
    window,
    novel_gene_prefix,
    isoforms_by_chr,
    refs_1exon_by_chr,
    refs_exons_by_chr,
    junctions_by_chr,
    junctions_by_gene,
    start_ends_by_gene,
    genome_dict,
    indelsJunc,
    orfDict,
):

    # read coverage files if provided

    if coverage is not None:
        print("**** Reading Splice Junctions coverage files.", file=sys.stdout)
        SJcovNames, SJcovInfo = STARcov_parser(coverage)
        fields_junc_cur = FIELDS_JUNC + SJcovNames  # add the samples to the header
    else:
        SJcovNames, SJcovInfo = None, None
        print("Splice Junction Coverage files not provided.", file=sys.stdout)
        fields_junc_cur = FIELDS_JUNC

    if cage_peak is not None:
        print("**** Reading CAGE Peak data.", file=sys.stdout)
        cage_peak_obj = CAGEPeak(cage_peak)
    else:
        cage_peak_obj = None

    if polyA_peak is not None:
        print("**** Reading polyA Peak data.", file=sys.stdout)
        polya_peak_obj = PolyAPeak(polyA_peak)
    else:
        polya_peak_obj = None

    if polyA_motif_list is not None:
        print("**** Reading PolyA motif list.", file=sys.stdout)
        polyA_motif_list = []
        for line in open(polyA_motif_list):
            x = line.strip().upper().replace("U", "A")
            if any(s not in ("A", "T", "C", "G") for s in x):
                print(
                    f"PolyA motif must be A/T/C/G only! Saw: {x}. Abort!",
                    file=sys.stderr,
                )
                sys.exit(-1)
            polyA_motif_list.append(x)
    else:
        polyA_motif_list = None

    if phyloP_bed is not None:
        print("**** Reading PhyloP BED file.", file=sys.stdout)
        phyloP_reader = LazyBEDPointReader(phyloP_bed)
    else:
        phyloP_reader = None

    # running classification
    print("**** Performing Classification of Isoforms....", file=sys.stdout)

    accepted_canonical_sites = list(sites.split(","))

    handle_class = open(outputClassPath + "_tmp", "w")
    fout_class = DictWriter(handle_class, fieldnames=FIELDS_CLASS, delimiter="\t")
    fout_class.writeheader()

    # outputJuncPath = outputPathPrefix+"_junctions.txt"
    handle_junc = open(outputJuncPath + "_tmp", "w")
    fout_junc = DictWriter(handle_junc, fieldnames=fields_junc_cur, delimiter="\t")
    fout_junc.writeheader()

    isoforms_info = {}
    novel_gene_index = 1

    for chrom, records in isoforms_by_chr.items():
        for rec in records:
            # Find best reference hit
            isoform_hit = transcriptsKnownSpliceSites(
                refs_1exon_by_chr,
                refs_exons_by_chr,
                start_ends_by_gene,
                rec,
                genome_dict,
                nPolyA=window,
            )

            if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
                # not FSM or ISM --> see if it is NIC, NNC, or fusion
                isoform_hit = novelIsoformsKnownGenes(
                    isoform_hit,
                    rec,
                    junctions_by_chr,
                    junctions_by_gene,
                    start_ends_by_gene,
                )
            elif isoform_hit.str_class in ("", "geneOverlap"):
                # possibly NNC, genic, genic intron, anti-sense, or intergenic
                isoform_hit = associationOverlapping(isoform_hit, rec, junctions_by_chr)

            # write out junction information
            write_junctionInfo(
                rec,
                junctions_by_chr,
                accepted_canonical_sites,
                indelsJunc,
                genome_dict,
                fout_junc,
                covInf=SJcovInfo,
                covNames=SJcovNames,
                phyloP_reader=phyloP_reader,
            )

            if isoform_hit.str_class in ("intergenic", "genic_intron"):
                # Liz: I don't find it necessary to cluster these novel genes. They should already be always non-overlapping.
                if (
                    novel_gene_prefix is not None
                ):  # used by splits to not have redundant novelGene IDs
                    isoform_hit.genes = [
                        "novelGene_"
                        + str(novel_gene_prefix)
                        + "_"
                        + str(novel_gene_index)
                    ]
                else:
                    isoform_hit.genes = ["novelGene_" + str(novel_gene_index)]
                isoform_hit.transcripts = ["novel"]
                novel_gene_index += 1

            # look at Cage Peak info (if available)
            if cage_peak_obj is not None:
                if rec.strand == "+":
                    within_cage, dist_cage = cage_peak_obj.find(
                        rec.chrom, rec.strand, rec.txStart
                    )
                else:
                    within_cage, dist_cage = cage_peak_obj.find(
                        rec.chrom, rec.strand, rec.txEnd
                    )
                isoform_hit.within_cage = within_cage
                isoform_hit.dist_cage = dist_cage

            # look at PolyA Peak info (if available)
            if polya_peak_obj is not None:
                if rec.strand == "+":
                    within_polya_site, dist_polya_site = polya_peak_obj.find(
                        rec.chrom, rec.strand, rec.txStart
                    )
                else:
                    within_polya_site, dist_polya_site = polya_peak_obj.find(
                        rec.chrom, rec.strand, rec.txEnd
                    )
                isoform_hit.within_polya_site = within_polya_site
                isoform_hit.dist_polya_site = dist_polya_site

            # polyA motif finding: look within 50 bp upstream of 3' end for the highest ranking polyA motif signal (user provided)
            if polyA_motif_list is not None:
                if rec.strand == "+":
                    polyA_motif, polyA_dist = find_polyA_motif(
                        str(genome_dict[rec.chrom][rec.txEnd - 50 : rec.txEnd].seq),
                        polyA_motif_list,
                    )
                else:
                    polyA_motif, polyA_dist = find_polyA_motif(
                        str(
                            genome_dict[rec.chrom][rec.txStart : rec.txStart + 50]
                            .reverse_complement()
                            .seq
                        ),
                        polyA_motif_list,
                    )
                isoform_hit.polyA_motif = polyA_motif
                isoform_hit.polyA_dist = polyA_dist

            # Fill in ORF/coding info and NMD detection
            if rec.id in orfDict:
                isoform_hit.coding = "coding"
                isoform_hit.ORFlen = orfDict[rec.id].orf_length
                isoform_hit.CDS_start = orfDict[rec.id].cds_start  # 1-based start
                isoform_hit.CDS_end = orfDict[rec.id].cds_end  # 1-based end

                m = {}  # transcript coord (0-based) --> genomic coord (0-based)
                if rec.strand == "+":
                    i = 0
                    for exon in rec.exons:
                        for c in range(exon.start, exon.end):
                            m[i] = c
                            i += 1
                else:  # - strand
                    i = 0
                    for exon in rec.exons:
                        for c in range(exon.start, exon.end):
                            m[rec.length - i - 1] = c
                            i += 1

                orfDict[rec.id].cds_genomic_start = (
                    m[orfDict[rec.id].cds_start - 1] + 1
                )  # make it 1-based
                orfDict[rec.id].cds_genomic_end = (
                    m[orfDict[rec.id].cds_end - 1] + 1
                )  # make it 1-based

                isoform_hit.CDS_genomic_start = orfDict[rec.id].cds_genomic_start
                isoform_hit.CDS_genomic_end = orfDict[rec.id].cds_genomic_end
                if (
                    orfDict[rec.id].cds_genomic_start is None
                ):  # likely SAM CIGAR mapping issue coming from aligner
                    continue  # we have to skip the NMD
                # NMD detection
                # if + strand, see if CDS stop is before the last junction
                if len(rec.junctions) > 0:
                    if rec.strand == "+":
                        dist_to_last_junc = (
                            orfDict[rec.id].cds_genomic_end - rec.junctions[-1][0]
                        )
                    else:  # - strand
                        dist_to_last_junc = (
                            rec.junctions[0][1] - orfDict[rec.id].cds_genomic_end
                        )
                    isoform_hit.is_NMD = "TRUE" if dist_to_last_junc < 0 else "FALSE"

            isoforms_info[rec.id] = isoform_hit
            fout_class.writerow(isoform_hit.as_dict())

    handle_class.close()
    handle_junc.close()
    return isoforms_info
