#!/usr/bin/env python
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila
# Modified by Liz (etseng@pacb.com) as SQANTI2/3 versions
# Modified by Fran (francisco.pardo.palacios@gmail.com) currently as SQANTI3 version (05/15/2020)

__author__ = "etseng@pacb.com"
# __version__ = "1.0.0"  # Python 3.7

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
import pygmst
from Bio import SeqIO
from bx.intervals import Interval, IntervalTree
from cupcake.cupcake.tofu.compare_junctions import compare_junctions
from cupcake.sequence.BED import LazyBEDPointReader
from cupcake.sequence.err_correct_w_genome import err_correct
from cupcake.sequence.GFF import collapseGFFReader, write_collapseGFF_format
from cupcake.sequence.sam_to_gff3 import convert_sam_to_gff3
from cupcake.sequence.STAR import STARJunctionReader
from pkg_resources import DistributionNotFound, get_distribution
from setuptools_scm import get_version

from sqanti3.utilities.indels_annot import calc_indels_from_sam
from sqanti3.utilities.rt_switching import rts

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

utilitiesPath = os.path.dirname(os.path.realpath(__file__)) + "/utilities/"
sys.path.insert(0, utilitiesPath)


GMAP_CMD = "gmap --cross-species -n 1 --max-intronlength-middle=2000000 --max-intronlength-ends=2000000 -L 3000000 -f samse -t {cpus} -D {dir} -d {name} -z {sense} {i} > {o}"
# MINIMAP2_CMD = "minimap2 -ax splice --secondary=no -C5 -O6,24 -B4 -u{sense} -t {cpus} {g} {i} > {o}"
MINIMAP2_CMD = (
    "minimap2 -ax splice --secondary=no -C5 -u{sense} -t {cpus} {g} {i} > {o}"
)
DESALT_CMD = "deSALT aln {dir} {i} -t {cpus} -x ccs -o {o}"

# GMSP_PROG = os.path.join(utilitiesPath, "gmst", "gmst.pl")
# GMST_CMD = "perl " + GMSP_PROG + " --strand direct --faa --fnn --output {o} {i}"

GTF2GENEPRED_PROG = distutils.spawn.find_executable("gtfToGenePred")
# GTF2GENEPRED_PROG = os.path.join(utilitiesPath, "gtfToGenePred")
GFFREAD_PROG = distutils.spawn.find_executable("gtfToGenePred")

if GTF2GENEPRED_PROG is None:
    print(f"Cannot find executable {GTF2GENEPRED_PROG}. Abort!", file=sys.stderr)
    sys.exit(-1)
if GFFREAD_PROG is None:
    print(f"Cannot find executable {GFFREAD_PROG}. Abort!", file=sys.stderr)
    sys.exit(-1)

seqid_rex1 = re.compile(r"PB\.(\d+)\.(\d+)$")
seqid_rex2 = re.compile(r"PB\.(\d+)\.(\d+)\|\S+")
seqid_fusion = re.compile(r"PBfusion\.(\d+)")

FIELDS_JUNC = [
    "isoform",
    "chrom",
    "strand",
    "junction_number",
    "genomic_start_coord",
    "genomic_end_coord",
    "transcript_coord",
    "junction_category",
    "start_site_category",
    "end_site_category",
    "diff_to_Ref_start_site",
    "diff_to_Ref_end_site",
    "bite_junction",
    "splice_site",
    "canonical",
    "RTS_junction",
    "indel_near_junct",
    "phyloP_start",
    "phyloP_end",
    "sample_with_cov",
    "total_coverage",
]  # +coverage_header

FIELDS_CLASS = [
    "isoform",
    "chrom",
    "strand",
    "length",
    "exons",
    "structural_category",
    "associated_gene",
    "associated_transcript",
    "ref_length",
    "ref_exons",
    "diff_to_TSS",
    "diff_to_TTS",
    "diff_to_gene_TSS",
    "diff_to_gene_TTS",
    "subcategory",
    "RTS_stage",
    "all_canonical",
    "min_sample_cov",
    "min_cov",
    "min_cov_pos",
    "sd_cov",
    "FL",
    "n_indels",
    "n_indels_junc",
    "bite",
    "iso_exp",
    "gene_exp",
    "ratio_exp",
    "FSM_class",
    "coding",
    "ORF_length",
    "CDS_length",
    "CDS_start",
    "CDS_end",
    "CDS_genomic_start",
    "CDS_genomic_end",
    "predicted_NMD",
    "perc_A_downstream_TTS",
    "seq_A_downstream_TTS",
    "dist_to_cage_peak",
    "within_cage_peak",
    "dist_to_polya_site",
    "within_polya_site",
    "polyA_motif",
    "polyA_dist",
]

RSCRIPTPATH = distutils.spawn.find_executable("Rscript")
RSCRIPT_REPORT = "SQANTI3_report.R"

if RSCRIPTPATH is None or os.system(RSCRIPTPATH + " --version") != 0:
    print("Rscript executable not found! Abort!", file=sys.stderr)
    sys.exit(-1)

SPLIT_ROOT_DIR = "splits/"


class genePredReader(object):
    def __init__(self, filename):
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return genePredRecord.from_line(line)


class genePredRecord(object):
    def __init__(
        self,
        id,
        chrom,
        strand,
        txStart,
        txEnd,
        cdsStart,
        cdsEnd,
        exonCount,
        exonStarts,
        exonEnds,
        gene=None,
    ):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.txStart = txStart  # 1-based start
        self.txEnd = txEnd  # 1-based end
        self.cdsStart = cdsStart  # 1-based start
        self.cdsEnd = cdsEnd  # 1-based end
        self.exonCount = exonCount
        self.exonStarts = exonStarts  # 0-based starts
        self.exonEnds = exonEnds  # 1-based ends
        self.gene = gene

        self.length = 0
        self.exons = []

        for s, e in zip(exonStarts, exonEnds):
            self.length += e - s
            self.exons.append(Interval(s, e))

        # junctions are stored (1-based last base of prev exon, 1-based first base of next exon)
        self.junctions = [
            (self.exonEnds[i], self.exonStarts[i + 1])
            for i in range(self.exonCount - 1)
        ]

    @property
    def segments(self):
        return self.exons

    @classmethod
    def from_line(cls, line):
        raw = line.strip().split("\t")
        return cls(
            id=raw[0],
            chrom=raw[1],
            strand=raw[2],
            txStart=int(raw[3]),
            txEnd=int(raw[4]),
            cdsStart=int(raw[5]),
            cdsEnd=int(raw[6]),
            exonCount=int(raw[7]),
            exonStarts=[
                int(x) for x in raw[8][:-1].split(",")
            ],  # exonStarts string has extra , at end
            exonEnds=[
                int(x) for x in raw[9][:-1].split(",")
            ],  # exonEnds string has extra , at end
            gene=raw[11] if len(raw) >= 12 else None,
        )

    def get_splice_site(self, genome_dict, i):
        """
        Return the donor-acceptor site (ex: GTAG) for the i-th junction
        :param i: 0-based junction index
        :param genome_dict: dict of chrom --> SeqRecord
        :return: splice site pattern, ex: "GTAG", "GCAG" etc
        """
        assert 0 <= i < self.exonCount - 1

        d = self.exonEnds[i]
        a = self.exonStarts[i + 1]

        seq_d = genome_dict[self.chrom].seq[d : d + 2]
        seq_a = genome_dict[self.chrom].seq[a - 2 : a]

        if self.strand == "+":
            return (str(seq_d) + str(seq_a)).upper()
        else:
            return (
                str(seq_a.reverse_complement()) + str(seq_d.reverse_complement())
            ).upper()


class myQueryTranscripts:
    def __init__(
        self,
        id,
        tss_diff,
        tts_diff,
        num_exons,
        length,
        str_class,
        subtype=None,
        genes=None,
        transcripts=None,
        chrom=None,
        strand=None,
        bite="NA",
        RT_switching="????",
        canonical="NA",
        min_cov="NA",
        min_cov_pos="NA",
        min_samp_cov="NA",
        sd="NA",
        FL="NA",
        FL_dict={},
        nIndels="NA",
        nIndelsJunc="NA",
        proteinID=None,
        ORFlen="NA",
        CDS_start="NA",
        CDS_end="NA",
        CDS_genomic_start="NA",
        CDS_genomic_end="NA",
        is_NMD="NA",
        isoExp="NA",
        geneExp="NA",
        coding="non_coding",
        refLen="NA",
        refExons="NA",
        refStart="NA",
        refEnd="NA",
        q_splicesite_hit=0,
        q_exon_overlap=0,
        FSM_class=None,
        percAdownTTS=None,
        seqAdownTTS=None,
        dist_cage="NA",
        within_cage="NA",
        dist_polya_site="NA",
        within_polya_site="NA",
        polyA_motif="NA",
        polyA_dist="NA",
    ):

        self.id = id
        self.tss_diff = tss_diff  # distance to TSS of best matching ref
        self.tts_diff = tts_diff  # distance to TTS of best matching ref
        self.tss_gene_diff = "NA"  # min distance to TSS of all genes matching the ref
        self.tts_gene_diff = "NA"  # min distance to TTS of all genes matching the ref
        self.genes = genes if genes is not None else []
        self.AS_genes = set()  # ref genes that are hit on the opposite strand
        self.transcripts = transcripts if transcripts is not None else []
        self.num_exons = num_exons
        self.length = length
        self.str_class = str_class  # structural classification of the isoform
        self.chrom = chrom
        self.strand = strand
        self.subtype = subtype
        self.RT_switching = RT_switching
        self.canonical = canonical
        self.min_samp_cov = min_samp_cov
        self.min_cov = min_cov
        self.min_cov_pos = min_cov_pos
        self.sd = sd
        self.proteinID = proteinID
        self.ORFlen = ORFlen
        self.CDS_start = CDS_start
        self.CDS_end = CDS_end
        self.coding = coding
        self.CDS_genomic_start = (
            CDS_genomic_start  # 1-based genomic coordinate of CDS start - strand aware
        )
        self.CDS_genomic_end = (
            CDS_genomic_end  # 1-based genomic coordinate of CDS end - strand aware
        )
        self.is_NMD = is_NMD  # (TRUE,FALSE) for NMD if is coding, otherwise "NA"
        self.FL = FL  # count for a single sample
        self.FL_dict = FL_dict  # dict of sample -> FL count
        self.nIndels = nIndels
        self.nIndelsJunc = nIndelsJunc
        self.isoExp = isoExp
        self.geneExp = geneExp
        self.refLen = refLen
        self.refExons = refExons
        self.refStart = refStart
        self.refEnd = refEnd
        self.q_splicesite_hit = q_splicesite_hit
        self.q_exon_overlap = q_exon_overlap
        self.FSM_class = FSM_class
        self.bite = bite
        self.percAdownTTS = percAdownTTS
        self.seqAdownTTS = seqAdownTTS
        self.dist_cage = dist_cage
        self.within_cage = within_cage
        self.within_polya_site = within_polya_site
        self.dist_polya_site = dist_polya_site  # distance to the closest polyA site (--polyA_peak, BEF file)
        self.polyA_motif = polyA_motif
        self.polyA_dist = polyA_dist  # distance to the closest polyA motif (--polyA_motif_list, 6mer motif list)

    def get_total_diff(self):
        return abs(self.tss_diff) + abs(self.tts_diff)

    def modify(self, ref_transcript, ref_gene, tss_diff, tts_diff, refLen, refExons):
        self.transcripts = [ref_transcript]
        self.genes = [ref_gene]
        self.tss_diff = tss_diff
        self.tts_diff = tts_diff
        self.refLen = refLen
        self.refExons = refExons

    def geneName(self):
        geneName = "_".join(set(self.genes))
        return geneName

    def ratioExp(self):
        if self.geneExp == 0 or self.geneExp == "NA":
            return "NA"
        else:
            ratio = float(self.isoExp) / float(self.geneExp)
        return ratio

    def CDSlen(self):
        if self.coding == "coding":
            return str(int(self.CDS_end) - int(self.CDS_start) + 1)
        else:
            return "NA"

    def __str__(self):
        tab = "\t"
        return (
            f"{self.chrom}{tab}"
            f"{self.strand}{tab}"
            f"{str(self.length)}{tab}"
            f"{str(self.num_exons)}{tab}"
            f"{str(self.str_class)}{tab}"
            f"{'_'.join(set(self.genes))}{tab}"
            f"{self.id}{tab}"
            f"{str(self.refLen)}{tab}"
            f"{str(self.refExons)}{tab}"
            f"{str(self.tss_diff)}{tab}"
            f"{str(self.tts_diff)}{tab}"
            f"{self.subtype}{tab}"
            f"{self.RT_switching}{tab}"
            f"{self.canonical}{tab}"
            f"{str(self.min_samp_cov)}{tab}"
            f"{str(self.min_cov)}{tab}"
            f"{str(self.min_cov_pos)}{tab}"
            f"{str(self.sd)}{tab}"
            f"{str(self.FL)}{tab}"
            f"{str(self.nIndels)}{tab}"
            f"{str(self.nIndelsJunc)}{tab}"
            f"{self.bite}{tab}"
            f"{str(self.isoExp)}{tab}"
            f"{str(self.geneExp)}{tab}"
            f"{str(self.ratioExp())}{tab}"
            f"{self.FSM_class}{tab}"
            f"{self.coding}{tab}"
            f"{str(self.ORFlen)}{tab}"
            f"{str(self.CDSlen())}{tab}"
            f"{str(self.CDS_start)}{tab}"
            f"{str(self.CDS_end)}{tab}"
            f"{str(self.CDS_genomic_start)}{tab}"
            f"{str(self.CDS_genomic_end)}{tab}"
            f"{str(self.is_NMD)}{tab}"
            f"{str(self.percAdownTTS)}{tab}"
            f"{str(self.seqAdownTTS)}{tab}"
            f"{str(self.dist_cage)}{tab}"
            f"{str(self.within_cage)}{tab}"
            f"{str(self.dist_polya_site)}{tab}"
            f"{str(self.within_polya_site)}{tab}"
            f"{str(self.polyA_motif)}{tab}"
            f"{str(self.polyA_dist)}{tab}"
        )

    def as_dict(self):
        d = {
            "isoform": self.id,
            "chrom": self.chrom,
            "strand": self.strand,
            "length": self.length,
            "exons": self.num_exons,
            "structural_category": self.str_class,
            "associated_gene": "_".join(set(self.genes)),
            "associated_transcript": "_".join(set(self.transcripts)),
            "ref_length": self.refLen,
            "ref_exons": self.refExons,
            "diff_to_TSS": self.tss_diff,
            "diff_to_TTS": self.tts_diff,
            "diff_to_gene_TSS": self.tss_gene_diff,
            "diff_to_gene_TTS": self.tts_gene_diff,
            "subcategory": self.subtype,
            "RTS_stage": self.RT_switching,
            "all_canonical": self.canonical,
            "min_sample_cov": self.min_samp_cov,
            "min_cov": self.min_cov,
            "min_cov_pos": self.min_cov_pos,
            "sd_cov": self.sd,
            "FL": self.FL,
            "n_indels": self.nIndels,
            "n_indels_junc": self.nIndelsJunc,
            "bite": self.bite,
            "iso_exp": self.isoExp,
            "gene_exp": self.geneExp,
            "ratio_exp": self.ratioExp(),
            "FSM_class": self.FSM_class,
            "coding": self.coding,
            "ORF_length": self.ORFlen,
            "CDS_length": self.CDSlen(),
            "CDS_start": self.CDS_start,
            "CDS_end": self.CDS_end,
            "CDS_genomic_start": self.CDS_genomic_start,
            "CDS_genomic_end": self.CDS_genomic_end,
            "predicted_NMD": self.is_NMD,
            "perc_A_downstream_TTS": self.percAdownTTS,
            "seq_A_downstream_TTS": self.seqAdownTTS,
            "dist_to_cage_peak": self.dist_cage,
            "within_cage_peak": self.within_cage,
            "dist_to_polya_site": self.dist_polya_site,
            "within_polya_site": self.within_polya_site,
            "polyA_motif": self.polyA_motif,
            "polyA_dist": self.polyA_dist,
        }
        for sample, count in self.FL_dict.items():
            d["FL." + sample] = count
        return d


class myQueryProteins:
    def __init__(self, cds_start, cds_end, orf_length, proteinID="NA"):
        self.orf_length = orf_length
        self.cds_start = cds_start  # 1-based start on transcript
        self.cds_end = cds_end  # 1-based end on transcript (stop codon), ORF is seq[cds_start-1:cds_end].translate()
        self.cds_genomic_start = (
            None  # 1-based genomic start of ORF, if - strand, is greater than end
        )
        self.cds_genomic_end = None  # 1-based genomic end of ORF
        self.proteinID = proteinID


def rewrite_sam_for_fusion_ids(sam_filename):
    seen_id_counter = Counter()

    f = open(sam_filename + ".tmp", "w")
    for line in open(sam_filename):
        if line.startswith("@"):
            f.write(line)
        else:
            raw = line.strip().split("\t")
            if not raw[0].startswith("PBfusion."):
                print(
                    f"Expecting fusion ID format `PBfusion.X` but saw {raw[0]} instead. Abort!",
                    file=sys.stderr,
                )
                sys.exit(-1)
            seen_id_counter[raw[0]] += 1
            raw[0] = raw[0] + "." + str(seen_id_counter[raw[0]])
            f.write("\t".join(raw) + "\n")
    f.close()
    os.rename(f.name, sam_filename)
    return sam_filename


def write_collapsed_GFF_with_CDS(isoforms_info, input_gff, output_gff):
    """
    Augment a collapsed GFF with CDS information
    *NEW* Also, change the "gene_id" field to use the classification result
    :param isoforms_info: dict of id -> QueryTranscript
    :param input_gff:  input GFF filename
    :param output_gff: output GFF filename
    """
    with open(output_gff, "w") as f:
        reader = collapseGFFReader(input_gff)
        for r in reader:
            r.geneid = isoforms_info[r.seqid].geneName()  # set the gene name

            s = isoforms_info[r.seqid].CDS_genomic_start  # could be 'NA'
            e = isoforms_info[r.seqid].CDS_genomic_end  # could be 'NA'
            r.cds_exons = []
            if s != "NA" and e != "NA":  # has ORF prediction for this isoform
                if r.strand == "+":
                    assert s < e
                    s = s - 1  # make it 0-based
                else:
                    assert e < s
                    s, e = e, s
                    s = s - 1  # make it 0-based
                for i, exon in enumerate(r.ref_exons):
                    if exon.end > s:
                        break
                r.cds_exons = [Interval(s, min(e, exon.end))]
                for exon in r.ref_exons[i + 1 :]:
                    if exon.start > e:
                        break
                    r.cds_exons.append(Interval(exon.start, min(e, exon.end)))
            write_collapseGFF_format(f, r)


def get_corr_filenames(args, dir=None):
    d = dir if dir is not None else args.dir
    corrPathPrefix = os.path.join(d, args.output)
    corrGTF = corrPathPrefix + "_corrected.gtf"
    corrSAM = corrPathPrefix + "_corrected.sam"
    corrFASTA = corrPathPrefix + "_corrected.fasta"
    corrORF = corrPathPrefix + "_corrected.faa"
    return corrGTF, corrSAM, corrFASTA, corrORF


def get_class_junc_filenames(args, dir=None):
    d = dir if dir is not None else args.dir
    outputPathPrefix = os.path.join(d, args.output)
    outputClassPath = outputPathPrefix + "_classification.txt"
    outputJuncPath = outputPathPrefix + "_junctions.txt"
    return outputClassPath, outputJuncPath


def correctionPlusORFpred(args, genome_dict):
    """
    Use the reference genome to correct the sequences (unless a pre-corrected GTF is given)
    """
    global corrORF
    global corrGTF
    global corrSAM
    global corrFASTA

    corrGTF, corrSAM, corrFASTA, corrORF = get_corr_filenames(args)

    n_cpu = max(1, args.cpus // args.chunks)

    # Step 1. IF GFF or GTF is provided, make it into a genome-based fasta
    #         IF sequence is provided, align as SAM then correct with genome
    if os.path.exists(corrFASTA):
        print(
            f"Error corrected FASTA {corrFASTA} already exists. Using it...",
            file=sys.stderr,
        )
    else:
        if not args.gtf:
            if os.path.exists(corrSAM):
                print(
                    f"Aligned SAM {corrSAM} already exists. Using it...",
                    file=sys.stderr,
                )
            else:
                if args.aligner_choice == "gmap":
                    print("****Aligning reads with GMAP...", file=sys.stdout)
                    cmd = GMAP_CMD.format(
                        cpus=n_cpu,
                        dir=os.path.dirname(args.gmap_index),
                        name=os.path.basename(args.gmap_index),
                        sense=args.sense,
                        i=args.isoforms,
                        o=corrSAM,
                    )
                elif args.aligner_choice == "minimap2":
                    print("****Aligning reads with Minimap2...", file=sys.stdout)
                    cmd = MINIMAP2_CMD.format(
                        cpus=n_cpu,
                        sense=args.sense,
                        g=args.genome,
                        i=args.isoforms,
                        o=corrSAM,
                    )
                elif args.aligner_choice == "deSALT":
                    print("****Aligning reads with deSALT...", file=sys.stdout)
                    cmd = DESALT_CMD.format(
                        cpus=n_cpu, dir=args.gmap_index, i=args.isoforms, o=corrSAM
                    )
                try:
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as error:
                    print(f"{error}\n {error.output}")
                    sys.exit(-1)

            # if is fusion - go in and change the IDs to reflect PBfusion.X.1, PBfusion.X.2...
            if args.is_fusion:
                corrSAM = rewrite_sam_for_fusion_ids(corrSAM)
            # error correct the genome (input: corrSAM, output: corrFASTA)
            err_correct(args.genome, corrSAM, corrFASTA, genome_dict=genome_dict)
            # convert SAM to GFF --> GTF
            print(f"corrSAM: {corrSAM}")
            print(f"corrFASTA: {corrFASTA}")
            print(f"corrGTF: {corrGTF}")
            print(f"args.genome: {args.genome}")

            convert_sam_to_gff3(
                corrSAM,
                corrGTF + ".tmp",
                source=os.path.basename(args.genome).split(".")[0],
            )  # convert SAM to GFF3
            cmd = f"{GFFREAD_PROG} {corrGTF}.tmp -T -o {corrGTF}"
            if subprocess.check_call(cmd, shell=True) != 0:
                print(f"ERROR running cmd: {cmd}", file=sys.stderr)
                sys.exit(-1)
        else:
            print(
                "Skipping aligning of sequences because GTF file was provided.",
                file=sys.stdout,
            )

            ind = 0
            with open(args.isoforms, "r") as isoforms_gtf:
                for line in isoforms_gtf:
                    if line[0] != "#" and len(line.split("\t")) != 9:
                        sys.stderr.write(
                            "\nERROR: input isoforms file with not GTF format.\n"
                        )
                        sys.exit()
                    elif len(line.split("\t")) == 9:
                        ind += 1
                if ind == 0:
                    print(
                        f"WARNING: GTF has {args.isoforms} no annotation lines.",
                        file=sys.stderr,
                    )

            # GFF to GTF (in case the user provides gff instead of gtf)
            corrGTF_tpm = corrGTF + ".tmp"
            try:
                subprocess.run([GFFREAD_PROG, args.isoforms, "-T", "-o", corrGTF_tpm])
            except (RuntimeError, TypeError, NameError):
                sys.stderr.write(
                    "ERROR: File %s without GTF/GFF format.\n" % args.isoforms
                )
                raise SystemExit(1)

            # check if gtf chromosomes inside genome file
            with open(corrGTF, "w") as corrGTF_out:
                with open(corrGTF_tpm, "r") as isoforms_gtf:
                    for line in isoforms_gtf:
                        if line[0] != "#":
                            chrom = line.split("\t")[0]
                            type = line.split("\t")[2]
                            if chrom not in list(genome_dict.keys()):
                                sys.stderr.write(
                                    '\nERROR: gtf "%s" chromosome not found in genome reference file.\n'
                                    % (chrom)
                                )
                                sys.exit()
                            elif type in ("transcript", "exon"):
                                corrGTF_out.write(line)
            os.remove(corrGTF_tpm)

            if not os.path.exists(corrSAM):
                sys.stdout.write(
                    "\nIndels will be not calculated since you ran SQANTI3 without alignment step (SQANTI3 with gtf format as transcriptome input).\n"
                )

            # GTF to FASTA
            subprocess.run([GFFREAD_PROG, corrGTF, "-g", args.genome, "-w", corrFASTA])

    # ORF generation
    print("**** Predicting ORF sequences...", file=sys.stdout)

    gmst_dir = os.path.join(os.path.abspath(args.dir), "GMST")
    gmst_pre = os.path.join(gmst_dir, "GMST_tmp")
    if not os.path.exists(gmst_dir):
        os.makedirs(gmst_dir)

    # sequence ID example: PB.2.1 gene_4|GeneMark.hmm|264_aa|+|888|1682
    gmst_rex = re.compile(r"(\S+\t\S+\|GeneMark.hmm)\|(\d+)_aa\|(\S)\|(\d+)\|(\d+)")
    orfDict = {}  # GMST seq id --> myQueryProteins object
    if args.skipORF:
        print(
            "WARNING: Skipping ORF prediction because user requested it. All isoforms will be non-coding!",
            file=sys.stderr,
        )
    elif os.path.exists(corrORF):
        print(f"ORF file {corrORF} already exists. Using it....", file=sys.stderr)
        for r in SeqIO.parse(open(corrORF), "fasta"):
            # now process ORFs into myQueryProtein objects
            m = gmst_rex.match(r.description)
            if m is None:
                print(
                    f"Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {r.description}! Abort!",
                    file=sys.stderr,
                )
                sys.exit(-1)
            orf_length = int(m.group(2))
            cds_start = int(m.group(4))
            cds_end = int(m.group(5))
            orfDict[r.id] = myQueryProteins(
                cds_start, cds_end, orf_length, proteinID=r.id
            )
    else:
        cur_dir = os.path.abspath(os.getcwd())
        os.chdir(args.dir)
        print(corrFASTA)
        pygmst.gmst(
            seqfile=corrFASTA, output=gmst_pre, faa=True, fnn=True, strand="direct"
        )
        # cmd = GMST_CMD.format(i=corrFASTA, o=gmst_pre)
        # try:
        #     subprocess.run(
        #         cmd,
        #         shell=True,
        #         check=True,
        #         stdout=subprocess.PIPE,
        #         stderr=subprocess.PIPE,
        #     )
        # except subprocess.CalledProcessError as error:
        #     print(error, subprocess.PIPE)
        #     sys.exit(-1)
        os.chdir(cur_dir)
        # Modifying ORF sequences by removing sequence before ATG
        with open(corrORF, "w") as f:
            for r in SeqIO.parse(open(gmst_pre + ".faa"), "fasta"):
                m = gmst_rex.match(r.description)
                if m is None:
                    print(
                        f"Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {r.description}! Abort!",
                        file=sys.stderr,
                    )
                    sys.exit(-1)
                id_pre = m.group(1)
                orf_length = int(m.group(2))
                orf_strand = m.group(3)
                cds_start = int(m.group(4))
                cds_end = int(m.group(5))
                pos = r.seq.find("M")
                if pos != -1:
                    # must modify both the sequence ID and the sequence
                    orf_length -= pos
                    cds_start += pos * 3
                    newid = (
                        f"{id_pre}|{orf_length}_aa|{orf_strand}|{cds_start}|{cds_end}"
                    )
                    newseq = str(r.seq)[pos:]
                    orfDict[r.id] = myQueryProteins(
                        cds_start, cds_end, orf_length, proteinID=newid
                    )
                    f.write(f">{newid}\n{newseq}\n")
                else:
                    new_rec = r
                    orfDict[r.id] = myQueryProteins(
                        cds_start, cds_end, orf_length, proteinID=r.id
                    )
                    f.write(f">{new_rec.description}\n{new_rec.seq}\n")

    if len(orfDict) == 0:
        print(
            "WARNING: All input isoforms were predicted as non-coding", file=sys.stderr
        )

    return orfDict


def reference_parser(args, genome_chroms):
    """
    Read the reference GTF file
    :param args:
    :param genome_chroms: list of chromosome names from the genome fasta, used for sanity checking
    :return: (refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene)
    """
    global referenceFiles

    referenceFiles = os.path.join(
        args.dir, "refAnnotation_" + args.output + ".genePred"
    )
    print("**** Parsing Reference Transcriptome....", file=sys.stdout)

    if os.path.exists(referenceFiles):
        print(f"{referenceFiles} already exists. Using it.", file=sys.stdout)
    else:
        # gtf to genePred
        if not args.genename:
            cmd = [
                GTF2GENEPRED_PROG,
                args.annotation,
                referenceFiles,
                "-genePredExt",
                "-allErrors",
                "-ignoreGroupsWithoutExons",
            ]
        else:
            cmd = [
                GTF2GENEPRED_PROG,
                args.annotation,
                referenceFiles,
                "-genePredExt",
                "-allErrors",
                "-ignoreGroupsWithoutExons",
                "-geneNameAsName2",
            ]
        try:
            subprocess.run(cmd, capture_output=True)
        except subprocess.CalledProcessError as error:
            print(error, file=subprocess.STDERR)
            sys.exit(-1)

    # parse reference annotation
    # 1. ignore all miRNAs (< 200 bp)
    # 2. separately store single exon and multi-exon references
    refs_1exon_by_chr = defaultdict(lambda: IntervalTree())  #
    refs_exons_by_chr = defaultdict(lambda: IntervalTree())
    # store donors as the exon end (1-based) and acceptor as the exon start (0-based)
    # will convert the sets to sorted list later
    junctions_by_chr = defaultdict(
        lambda: {"donors": set(), "acceptors": set(), "da_pairs": set()}
    )
    # dict of gene name --> set of junctions (don't need to record chromosome)
    junctions_by_gene = defaultdict(lambda: set())
    # dict of gene name --> list of known begins and ends (begin always < end, regardless of strand)
    known_5_3_by_gene = defaultdict(lambda: {"begin": set(), "end": set()})

    for r in genePredReader(referenceFiles):
        if r.length < args.min_ref_len:
            continue  # ignore miRNAs
        if r.exonCount == 1:
            refs_1exon_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
            known_5_3_by_gene[r.gene]["begin"].add(r.txStart)
            known_5_3_by_gene[r.gene]["end"].add(r.txEnd)
        else:
            refs_exons_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
            # only store junctions for multi-exon transcripts
            for d, a in r.junctions:
                junctions_by_chr[r.chrom]["donors"].add(d)
                junctions_by_chr[r.chrom]["acceptors"].add(a)
                junctions_by_chr[r.chrom]["da_pairs"].add((d, a))
                junctions_by_gene[r.gene].add((d, a))
            known_5_3_by_gene[r.gene]["begin"].add(r.txStart)
            known_5_3_by_gene[r.gene]["end"].add(r.txEnd)

    # check that all genes' chromosomes are in the genome file
    ref_chroms = set(refs_1exon_by_chr.keys()).union(list(refs_exons_by_chr.keys()))
    diff = ref_chroms.difference(genome_chroms)
    if len(diff) > 0:
        print(
            f"WARNING: ref annotation contains chromosomes not in genome: {','.join(diff)}\n",
            file=sys.stderr,
        )

    # convert the content of junctions_by_chr to sorted list
    for k in junctions_by_chr:
        junctions_by_chr[k]["donors"] = list(junctions_by_chr[k]["donors"])
        junctions_by_chr[k]["donors"].sort()
        junctions_by_chr[k]["acceptors"] = list(junctions_by_chr[k]["acceptors"])
        junctions_by_chr[k]["acceptors"].sort()
        junctions_by_chr[k]["da_pairs"] = list(junctions_by_chr[k]["da_pairs"])
        junctions_by_chr[k]["da_pairs"].sort()

    return (
        dict(refs_1exon_by_chr),
        dict(refs_exons_by_chr),
        dict(junctions_by_chr),
        dict(junctions_by_gene),
        dict(known_5_3_by_gene),
    )


def isoforms_parser(args):
    """
    Parse input isoforms (GTF) to dict (chr --> sorted list)
    """
    global queryFile
    queryFile = os.path.splitext(corrGTF)[0] + ".genePred"

    print("**** Parsing Isoforms....", file=sys.stderr)

    # gtf to genePred
    cmd = (
        GTF2GENEPRED_PROG
        + f" {corrGTF} {queryFile} -genePredExt -allErrors -ignoreGroupsWithoutExons"
    )
    if subprocess.check_call(cmd, shell=True) != 0:
        print(f"ERROR running cmd: {cmd}", file=sys.stderr)
        sys.exit(-1)

    isoforms_list = defaultdict(lambda: [])  # chr --> list to be sorted later

    for r in genePredReader(queryFile):
        isoforms_list[r.chrom].append(r)

    for k in isoforms_list:
        isoforms_list[k].sort(key=lambda r: r.txStart)

    return isoforms_list


def STARcov_parser(
    coverageFiles,
):  # just valid with unstrand-specific RNA-seq protocols.
    """
    :param coverageFiles: comma-separated list of STAR junction output files or a directory containing junction files
    :return: list of samples, dict of (chrom,strand) --> (0-based start, 1-based end) --> {dict of sample -> unique reads supporting this junction}
    """
    if os.path.isdir(coverageFiles) is True:
        cov_files = glob.glob(coverageFiles + "/*SJ.out.tab")
    else:
        cov_files = coverageFiles.split(",")

    newline = "\n"
    print(
        f"Input pattern: {coverageFiles}. The following files found and to be read as junctions:\n{newline.join(cov_files)}",
        file=sys.stderr,
    )

    cov_by_chrom_strand = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: 0))
    )
    undefined_strand_count = 0
    all_read = 0
    samples = []
    for file in cov_files:
        prefix = os.path.basename(file[: file.rfind(".")])  # use this as sample name
        samples.append(prefix)
        for r in STARJunctionReader(file):
            if r.strand == "NA":
                # undefined strand, so we put them in BOTH strands otherwise we'll lose all non-canonical junctions from STAR
                cov_by_chrom_strand[(r.chrom, "+")][(r.start, r.end)][prefix] = (
                    r.unique_count + r.multi_count
                )
                cov_by_chrom_strand[(r.chrom, "-")][(r.start, r.end)][prefix] = (
                    r.unique_count + r.multi_count
                )
                undefined_strand_count += 1
            else:
                cov_by_chrom_strand[(r.chrom, r.strand)][(r.start, r.end)][prefix] = (
                    r.unique_count + r.multi_count
                )
            all_read += 1
    print(
        f"{all_read} junctions read. {undefined_strand_count} junctions added to both strands because no strand information from STAR.",
        file=sys.stderr,
    )

    return samples, cov_by_chrom_strand


EXP_KALLISTO_HEADERS = ["target_id", "length", "eff_length", "est_counts", "tpm"]
EXP_RSEM_HEADERS = [
    "transcript_id",
    "length",
    "effective_length",
    "expected_count",
    "TPM",
]


def mergeDict(dict1, dict2):
    """ Merge dictionaries to collect info from several files"""
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
            dict3[key] = [value, dict1[key]]
    return dict3


def flatten(lis):

    for item in lis:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def expression_parser(expressionFile):
    """
    Currently accepts expression format: Kallisto or RSEM
    :param expressionFile: Kallisto or RSEM
    :return: dict of PBID --> TPM
    Include the possibility of providing an expression matrix --> first column must be "ID"
    """
    if os.path.isdir(expressionFile) is True:
        exp_paths = [
            os.path.join(expressionFile, fn) for fn in next(os.walk(expressionFile))[2]
        ]
    else:
        exp_paths = expressionFile.split(",")
    exp_all = {}
    ismatrix = False
    for exp_file in exp_paths:
        reader = DictReader(open(exp_file), delimiter="\t")
        if all(k in reader.fieldnames for k in EXP_KALLISTO_HEADERS):
            print(
                "Detected Kallisto expression format. Using 'target_id' and 'tpm' field.",
                file=sys.stderr,
            )
            name_id, name_tpm = "target_id", "tpm"
        elif all(k in reader.fieldnames for k in EXP_RSEM_HEADERS):
            print(
                "Detected RSEM expression format. Using 'transcript_id' and 'TPM' field.",
                file=sys.stderr,
            )
            name_id, name_tpm = "transcript_id", "TPM"
        elif reader.fieldnames[0] == "ID":
            print("Detected expression matrix format")
            ismatrix = True
            name_id = "ID"
        else:
            print(
                f"Expected Kallisto or RSEM file format from {expressionFile}. Abort!",
                file=sys.stderr,
            )
        exp_sample = {}
        if ismatrix:
            for r in reader:
                exp_sample[r[name_id]] = np.average(
                    list(map(float, list(r.values())[1:]))
                )
        else:
            for r in reader:
                exp_sample[r[name_id]] = float(r[name_tpm])

        exp_all = mergeDict(exp_all, exp_sample)

    exp_dict = {}
    if len(exp_paths) > 1:
        for k in exp_all:
            exp_all[k] = list(flatten(exp_all[k]))
            exp_dict[k] = np.mean(exp_all[k])
        return exp_dict
    else:
        exp_dict = exp_all
        return exp_dict


def transcriptsKnownSpliceSites(
    refs_1exon_by_chr, refs_exons_by_chr, start_ends_by_gene, trec, genome_dict, nPolyA
):
    """
    :param refs_1exon_by_chr: dict of single exon references (chr -> IntervalTree)
    :param refs_exons_by_chr: dict of multi exon references (chr -> IntervalTree)
    :param trec: id record (genePredRecord) to be compared against reference
    :param genome_dict: dict of genome (chrom --> SeqRecord)
    :param nPolyA: window size to look for polyA
    :return: myQueryTranscripts object that indicates the best reference hit
    """

    def calc_overlap(s1, e1, s2, e2):
        if s1 == "NA" or s2 == "NA":
            return 0
        if s1 > s2:
            s1, e1, s2, e2 = s2, e2, s1, e1
        return max(0, min(e1, e2) - max(s1, s2))

    def gene_overlap(ref1, ref2):
        if ref1 == ref2:
            return True  # same gene, diff isoforms
        # return True if the two reference genes overlap
        s1, e1 = (
            min(start_ends_by_gene[ref1]["begin"]),
            max(start_ends_by_gene[ref1]["end"]),
        )
        s2, e2 = (
            min(start_ends_by_gene[ref2]["begin"]),
            max(start_ends_by_gene[ref2]["end"]),
        )
        if s1 <= s2:
            return e1 <= s2
        else:
            return e2 <= s1

    def calc_splicesite_agreement(query_exons, ref_exons):
        q_sites = {}
        for e in query_exons:
            q_sites[e.start] = 0
            q_sites[e.end] = 0
        for e in ref_exons:
            if e.start in q_sites:
                q_sites[e.start] = 1
            if e.end in q_sites:
                q_sites[e.end] = 1
        return sum(q_sites.values())

    def calc_exon_overlap(query_exons, ref_exons):
        q_bases = {}
        for e in query_exons:
            for b in range(e.start, e.end):
                q_bases[b] = 0

        for e in ref_exons:
            for b in range(e.start, e.end):
                if b in q_bases:
                    q_bases[b] = 1
        return sum(q_bases.values())

    def get_diff_tss_tts(trec, ref):
        if trec.strand == "+":
            diff_tss = trec.txStart - ref.txStart
            diff_tts = ref.txEnd - trec.txEnd
        else:
            diff_tts = trec.txStart - ref.txStart
            diff_tss = ref.txEnd - trec.txEnd
        return diff_tss, diff_tts

    def get_gene_diff_tss_tts(isoform_hit):
        # now that we know the reference (isoform) it hits
        # add the nearest start/end site for that gene (all isoforms of the gene)
        nearest_start_diff, nearest_end_diff = float("inf"), float("inf")
        for ref_gene in isoform_hit.genes:
            for x in start_ends_by_gene[ref_gene]["begin"]:
                d = trec.txStart - x
                if abs(d) < abs(nearest_start_diff):
                    nearest_start_diff = d
            for x in start_ends_by_gene[ref_gene]["end"]:
                d = trec.txEnd - x
                if abs(d) < abs(nearest_end_diff):
                    nearest_end_diff = d

        if trec.strand == "+":
            isoform_hit.tss_gene_diff = (
                nearest_start_diff if nearest_start_diff != float("inf") else "NA"
            )
            isoform_hit.tts_gene_diff = (
                nearest_end_diff if nearest_end_diff != float("inf") else "NA"
            )
        else:
            isoform_hit.tss_gene_diff = (
                -nearest_end_diff if nearest_start_diff != float("inf") else "NA"
            )
            isoform_hit.tts_gene_diff = (
                -nearest_start_diff if nearest_end_diff != float("inf") else "NA"
            )

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
        for i, e in enumerate(ref.exons):
            ref_exon_tree.insert(e.start, e.end, i)
        for e in trec.exons:
            if (
                len(ref_exon_tree.find(e.start, e.end)) > 1
            ):  # multiple ref exons covered
                return "intron_retention"

        agree_front = trec.junctions[0] == ref.junctions[0]
        agree_end = trec.junctions[-1] == ref.junctions[-1]
        if agree_front:
            if agree_end:
                return "complete"
            else:  # front agrees, end does not
                return "5prime_fragment" if trec.strand == "+" else "3prime_fragment"
        else:
            if agree_end:  # front does not agree, end agrees
                return "3prime_fragment" if trec.strand == "+" else "5prime_fragment"
            else:
                return "internal_fragment"

    # Transcript information for a single query id and comparison with reference.

    # Intra-priming: calculate percentage of "A"s right after the end
    if trec.strand == "+":
        pos_TTS = trec.exonEnds[-1]
        seq_downTTS = str(
            genome_dict[trec.chrom].seq[pos_TTS : pos_TTS + nPolyA]
        ).upper()
    else:  # id on - strand
        pos_TTS = trec.exonStarts[0]
        seq_downTTS = str(
            genome_dict[trec.chrom].seq[pos_TTS - nPolyA : pos_TTS].reverse_complement()
        ).upper()

    percA = float(seq_downTTS.count("A")) / nPolyA * 100

    isoform_hit = myQueryTranscripts(
        id=trec.id,
        tts_diff="NA",
        tss_diff="NA",
        num_exons=trec.exonCount,
        length=trec.length,
        str_class="",
        chrom=trec.chrom,
        strand=trec.strand,
        subtype="no_subcategory",
        percAdownTTS=str(percA),
        seqAdownTTS=seq_downTTS,
    )

    # SPLICED TRANSCRIPTS
    cat_ranking = {
        "full-splice_match": 5,
        "incomplete-splice_match": 4,
        "anyKnownJunction": 3,
        "anyKnownSpliceSite": 2,
        "geneOverlap": 1,
        "": 0,
    }

    # if trec.id.startswith('PB.1961.2'):
    #    pdb.set_trace()
    if trec.exonCount >= 2:

        hits_by_gene = defaultdict(lambda: [])  # gene --> list of hits
        best_by_gene = {}  # gene --> best isoform_hit

        if trec.chrom in refs_exons_by_chr:
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)

        if len(hits_by_gene) == 0:
            return isoform_hit

        for ref_gene in hits_by_gene:
            isoform_hit = myQueryTranscripts(
                id=trec.id,
                tts_diff="NA",
                tss_diff="NA",
                num_exons=trec.exonCount,
                length=trec.length,
                str_class="",
                chrom=trec.chrom,
                strand=trec.strand,
                subtype="no_subcategory",
                percAdownTTS=str(percA),
                seqAdownTTS=seq_downTTS,
            )

            for ref in hits_by_gene[ref_gene]:
                if trec.strand != ref.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue

                # if trec.id.startswith('PB.102.9'):
                #    pdb.set_trace()
                if ref.exonCount == 1:  # mono-exonic reference, handle specially here
                    if (
                        calc_exon_overlap(trec.exons, ref.exons) > 0
                        and cat_ranking[isoform_hit.str_class]
                        < cat_ranking["geneOverlap"]
                    ):
                        isoform_hit = myQueryTranscripts(
                            trec.id,
                            "NA",
                            "NA",
                            trec.exonCount,
                            trec.length,
                            "geneOverlap",
                            subtype="mono-exon",
                            chrom=trec.chrom,
                            strand=trec.strand,
                            genes=[ref.gene],
                            transcripts=[ref.id],
                            refLen=ref.length,
                            refExons=ref.exonCount,
                            refStart=ref.txStart,
                            refEnd=ref.txEnd,
                            q_splicesite_hit=0,
                            q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                            percAdownTTS=str(percA),
                            seqAdownTTS=seq_downTTS,
                        )

                else:  # multi-exonic reference
                    match_type = compare_junctions(
                        trec,
                        ref,
                        internal_fuzzy_max_dist=0,
                        max_5_diff=999999,
                        max_3_diff=999999,
                    )

                    if match_type not in (
                        "exact",
                        "subset",
                        "partial",
                        "concordant",
                        "super",
                        "nomatch",
                    ):
                        raise Exception(f"Unknown match category {match_type}!")

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
                    # has_overlap = gene_overlap(isoform_hit.genes[-1], ref.gene) if len(isoform_hit.genes) >= 1 else Fals
                    # #############################
                    # SQANTI's full-splice_match
                    # #############################
                    if match_type == "exact":
                        subtype = "multi-exon"
                        # assign as a new hit if
                        # (1) no prev hits yet
                        # (2) this one is better (prev not FSM or is FSM but worse tss/tts)
                        if (
                            cat_ranking[isoform_hit.str_class]
                            < cat_ranking["full-splice_match"]
                            or abs(diff_tss) + abs(diff_tts)
                            < isoform_hit.get_total_diff()
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                diff_tss,
                                diff_tts,
                                trec.exonCount,
                                trec.length,
                                str_class="full-splice_match",
                                subtype=subtype,
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=[ref.id],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                percAdownTTS=str(percA),
                                seqAdownTTS=seq_downTTS,
                            )
                    # #######################################################
                    # SQANTI's incomplete-splice_match
                    # (only check if don't already have a FSM match)
                    # #######################################################
                    elif match_type == "subset":
                        subtype = categorize_incomplete_matches(trec, ref)
                        # assign as a new (ISM) hit if
                        # (1) no prev hit
                        # (2) prev hit not as good (is ISM with worse tss/tts or anyKnownSpliceSite)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking[
                            "incomplete-splice_match"
                        ] or (
                            isoform_hit.str_class == "incomplete-splice_match"
                            and abs(diff_tss) + abs(diff_tts)
                            < isoform_hit.get_total_diff()
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                diff_tss,
                                diff_tts,
                                trec.exonCount,
                                trec.length,
                                str_class="incomplete-splice_match",
                                subtype=subtype,
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=[ref.id],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                percAdownTTS=str(percA),
                                seqAdownTTS=seq_downTTS,
                            )
                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ("partial", "concordant", "super"):
                        q_sp_hit = calc_splicesite_agreement(trec.exons, ref.exons)
                        q_ex_overlap = calc_exon_overlap(trec.exons, ref.exons)
                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                        if (
                            cat_ranking[isoform_hit.str_class]
                            < cat_ranking["anyKnownJunction"]
                            or (
                                isoform_hit.str_class == "anyKnownJunction"
                                and q_sp_hit > isoform_hit.q_splicesite_hit
                            )
                            or (
                                isoform_hit.str_class == "anyKnownJunction"
                                and q_sp_hit == isoform_hit.q_splicesite_hit
                                and q_ex_overlap > isoform_hit.q_exon_overlap
                            )
                            or (
                                isoform_hit.str_class == "anyKnownJunction"
                                and q_sp_hit == isoform_hit.q_splicesite_hit
                                and q_exon_d
                                < abs(trec.exonCount - isoform_hit.refExons)
                            )
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                "NA",
                                "NA",
                                trec.exonCount,
                                trec.length,
                                str_class="anyKnownJunction",
                                subtype="no_subcategory",
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=["novel"],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                percAdownTTS=str(percA),
                                seqAdownTTS=seq_downTTS,
                            )
                    else:  # must be nomatch
                        assert match_type == "nomatch"
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if (
                            cat_ranking[isoform_hit.str_class]
                            < cat_ranking["anyKnownSpliceSite"]
                            and calc_splicesite_agreement(trec.exons, ref.exons) > 0
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                "NA",
                                "NA",
                                trec.exonCount,
                                trec.length,
                                str_class="anyKnownSpliceSite",
                                subtype="no_subcategory",
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=["novel"],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                percAdownTTS=str(percA),
                                seqAdownTTS=seq_downTTS,
                            )

                        if (
                            isoform_hit.str_class == ""
                        ):  # still not hit yet, check exonic overlap
                            if (
                                cat_ranking[isoform_hit.str_class]
                                < cat_ranking["geneOverlap"]
                                and calc_exon_overlap(trec.exons, ref.exons) > 0
                            ):
                                isoform_hit = myQueryTranscripts(
                                    trec.id,
                                    "NA",
                                    "NA",
                                    trec.exonCount,
                                    trec.length,
                                    str_class="geneOverlap",
                                    subtype="no_subcategory",
                                    chrom=trec.chrom,
                                    strand=trec.strand,
                                    genes=[ref.gene],
                                    transcripts=["novel"],
                                    refLen=ref.length,
                                    refExons=ref.exonCount,
                                    refStart=ref.txStart,
                                    refEnd=ref.txEnd,
                                    q_splicesite_hit=calc_splicesite_agreement(
                                        trec.exons, ref.exons
                                    ),
                                    q_exon_overlap=calc_exon_overlap(
                                        trec.exons, ref.exons
                                    ),
                                    percAdownTTS=str(percA),
                                    seqAdownTTS=seq_downTTS,
                                )

            best_by_gene[ref_gene] = isoform_hit
        # now we have best_by_gene:
        # start with the best scoring one (FSM is best) --> can add other genes if they don't overlap
        # if trec.id.startswith('PB.1252.'):
        #    pdb.set_trace()
        geneHitTuple = namedtuple(
            "geneHitTuple", ["score", "rStart", "rEnd", "rGene", "iso_hit"]
        )
        best_by_gene = [
            geneHitTuple(
                cat_ranking[iso_hit.str_class],
                iso_hit.refStart,
                iso_hit.refEnd,
                ref_gene,
                iso_hit,
            )
            for ref_gene, iso_hit in best_by_gene.items()
        ]
        best_by_gene = list(filter(lambda x: x.score > 0, best_by_gene))
        if len(best_by_gene) == 0:  # no hit
            return isoform_hit

        # sort matching genes by ranking, allow for multi-gene match as long as they don't overlap
        # cat_ranking = {'full-splice_match': 5, 'incomplete-splice_match': 4, 'anyKnownJunction': 3, 'anyKnownSpliceSite': 2,
        #                    'geneOverlap': 1, '': 0}

        best_by_gene.sort(
            key=lambda x: (
                x.score,
                x.iso_hit.q_splicesite_hit
                + (x.iso_hit.q_exon_overlap)
                * 1.0
                / sum(e.end - e.start for e in trec.exons)
                + calc_overlap(x.rStart, x.rEnd, trec.txStart, trec.txEnd)
                * 1.0
                / (x.rEnd - x.rStart)
                - abs(trec.exonCount - x.iso_hit.refExons),
            ),
            reverse=True,
        )  # sort by (ranking score, overlap)
        isoform_hit = best_by_gene[0].iso_hit
        cur_start, cur_end = best_by_gene[0].rStart, best_by_gene[0].rEnd
        for t in best_by_gene[1:]:
            if t.score == 0:
                break
            if calc_overlap(cur_start, cur_end, t.rStart, t.rEnd) <= 0:
                isoform_hit.genes.append(t.rGene)
                cur_start, cur_end = min(cur_start, t.rStart), max(cur_end, t.rEnd)

    # UNSPLICED TRANSCRIPTS
    else:  # single exon id
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                # see if there's already an existing match AND if so, if this one is better
                if isoform_hit.str_class == "":  # no match so far
                    isoform_hit = myQueryTranscripts(
                        trec.id,
                        diff_tss,
                        diff_tts,
                        trec.exonCount,
                        trec.length,
                        "full-splice_match",
                        subtype="mono-exon",
                        chrom=trec.chrom,
                        strand=trec.strand,
                        genes=[ref.gene],
                        transcripts=[ref.id],
                        refLen=ref.length,
                        refExons=ref.exonCount,
                        percAdownTTS=str(percA),
                        seqAdownTTS=seq_downTTS,
                    )
                elif abs(diff_tss) + abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.modify(
                        ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount
                    )

        if isoform_hit.str_class == "" and trec.chrom in refs_exons_by_chr:
            # no hits to single exon genes, let's see if it hits multi-exon genes
            # (1) if it overlaps with a ref exon and is contained in an exon, we call it ISM
            # (2) else, if it is completely within a ref gene start-end region, we call it NIC by intron retention
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if (
                    calc_exon_overlap(trec.exons, ref.exons) == 0
                ):  # no exonic overlap, skip!
                    continue
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                for e in ref.exons:
                    if e.start <= trec.txStart < trec.txEnd <= e.end:
                        isoform_hit.str_class = "incomplete-splice_match"
                        isoform_hit.subtype = "mono-exon"
                        isoform_hit.modify(
                            ref.id,
                            ref.gene,
                            diff_tss,
                            diff_tts,
                            ref.length,
                            ref.exonCount,
                        )
                        # this is as good a match as it gets, we can stop the search here
                        get_gene_diff_tss_tts(isoform_hit)
                        return isoform_hit

                # if we haven't exited here, then ISM hit is not found yet
                # instead check if it's NIC by intron retention
                # but we don't exit here since the next gene could be a ISM hit
                if ref.txStart <= trec.txStart < trec.txEnd <= ref.txEnd:
                    isoform_hit.str_class = "novel_in_catalog"
                    isoform_hit.subtype = "mono-exon"
                    # check for intron retention
                    if len(ref.junctions) > 0:
                        for (d, a) in ref.junctions:
                            if trec.txStart < d < a < trec.txEnd:
                                isoform_hit.subtype = "mono-exon_by_intron_retention"
                                break
                    isoform_hit.modify(
                        "novel", ref.gene, "NA", "NA", ref.length, ref.exonCount
                    )
                    get_gene_diff_tss_tts(isoform_hit)
                    return isoform_hit

                # if we get to here, means neither ISM nor NIC, so just add a ref gene and categorize further later
                isoform_hit.genes.append(ref.gene)

    get_gene_diff_tss_tts(isoform_hit)
    isoform_hit.genes.sort(key=lambda x: start_ends_by_gene[x]["begin"])
    return isoform_hit


def novelIsoformsKnownGenes(
    isoforms_hit, trec, junctions_by_chr, junctions_by_gene, start_ends_by_gene
):
    """
    At this point: definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    :return isoforms_hit: updated isoforms hit (myQueryTranscripts object)
    """

    def has_intron_retention():
        for e in trec.exons:
            m = bisect.bisect_left(
                junctions_by_chr[trec.chrom]["da_pairs"], (e.start, e.end)
            )
            if (
                m < len(junctions_by_chr[trec.chrom]["da_pairs"])
                and e.start
                <= junctions_by_chr[trec.chrom]["da_pairs"][m][0]
                < junctions_by_chr[trec.chrom]["da_pairs"][m][1]
                < e.end
            ):
                return True
        return False

    ref_genes = list(set(isoforms_hit.genes))

    # if trec.id.startswith('PB.37872'):
    # pdb.set_trace()
    #
    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update refLen or refExon
    # or tss_diff and tts_diff (always set to "NA" for non-FSM/ISM matches)
    #
    isoforms_hit.transcripts = ["novel"]
    if len(ref_genes) == 1:
        # hits exactly one gene, must be either NIC or NNC
        ref_gene_junctions = junctions_by_gene[ref_genes[0]]
        # 1. check if all donors/acceptor sites are known (regardless of which ref gene it came from)
        # 2. check if this query isoform uses a subset of the junctions from the single ref hit
        all_junctions_known = True
        all_junctions_in_hit_ref = True
        for d, a in trec.junctions:
            all_junctions_known = (
                all_junctions_known
                and (d in junctions_by_chr[trec.chrom]["donors"])
                and (a in junctions_by_chr[trec.chrom]["acceptors"])
            )
            all_junctions_in_hit_ref = all_junctions_in_hit_ref and (
                (d, a) in ref_gene_junctions
            )
        if all_junctions_known:
            isoforms_hit.str_class = "novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoforms_hit.subtype = "combination_of_known_junctions"
            else:
                isoforms_hit.subtype = "combination_of_known_splicesites"
        else:
            isoforms_hit.str_class = "novel_not_in_catalog"
            isoforms_hit.subtype = "at_least_one_novel_splicesite"
    else:  # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        # NOTE: some ref genes could be mono-exonic so no junctions
        all_ref_junctions = list(
            itertools.chain(
                junctions_by_gene[ref_gene]
                for ref_gene in ref_genes
                if ref_gene in junctions_by_gene
            )
        )

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = {
            i: all_ref_junctions.count(junc) for i, junc in enumerate(trec.junctions)
        }

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.values()) > 1:
            isoforms_hit.str_class = "moreJunctions"
        else:
            isoforms_hit.str_class = "fusion"
            isoforms_hit.subtype = "mono-exon" if trec.exonCount == 1 else "multi-exon"

    if has_intron_retention():
        isoforms_hit.subtype = "intron_retention"

    return isoforms_hit


def associationOverlapping(isoforms_hit, trec, junctions_by_chr):
    # at this point: definitely not FSM or ISM or NIC or NNC
    # possibly (in order of preference assignment):
    #  - antisense  (on opp strand of a known gene)
    #  - genic (overlaps a combination of exons and introns), ignore strand
    #  - genic_intron  (completely within an intron), ignore strand
    #  - intergenic (does not overlap a gene at all), ignore strand

    isoforms_hit.str_class = "intergenic"
    isoforms_hit.transcripts = ["novel"]
    isoforms_hit.subtype = "mono-exon" if trec.exonCount == 1 else "multi-exon"

    # if trec.id.startswith('PB.37872.'):
    #    pdb.set_trace()
    if len(isoforms_hit.genes) == 0:
        # completely no overlap with any genes on the same strand
        # check if it is anti-sense to a known gene, otherwise it's genic_intron or intergenic
        if len(isoforms_hit.AS_genes) == 0 and trec.chrom in junctions_by_chr:
            # no hit even on opp strand
            # see if it is completely contained within a junction
            da_pairs = junctions_by_chr[trec.chrom]["da_pairs"]
            i = bisect.bisect_left(da_pairs, (trec.txStart, trec.txEnd))
            while i < len(da_pairs) and da_pairs[i][0] <= trec.txStart:
                if da_pairs[i][0] <= trec.txStart <= trec.txStart <= da_pairs[i][1]:
                    isoforms_hit.str_class = "genic_intron"
                    break
                i += 1
        else:
            # hits one or more genes on the opposite strand
            isoforms_hit.str_class = "antisense"
            isoforms_hit.genes = [f"novelGene_{g}_AS" for g in isoforms_hit.AS_genes]
    else:
        # (Liz) used to put NNC here - now just genic
        isoforms_hit.str_class = "genic"
        # overlaps with one or more genes on the same strand
        # if trec.exonCount >= 2:
        #    # multi-exon and has a same strand gene hit, must be NNC
        #    isoforms_hit.str_class = "novel_not_in_catalog"
        #    isoforms_hit.subtype = "at_least_one_novel_splicesite"
        # else:
        #    # single exon, must be genic
        #    isoforms_hit.str_class = "genic"

    return isoforms_hit


def write_junctionInfo(
    trec,
    junctions_by_chr,
    accepted_canonical_sites,
    indelInfo,
    genome_dict,
    fout,
    covInf=None,
    covNames=None,
    phyloP_reader=None,
):
    """
    :param trec: query isoform genePredRecord
    :param junctions_by_chr: dict of chr -> {'donors': <sorted list of donors>, 'acceptors': <sorted list of acceptors>, 'da_pairs': <sorted list of junctions>}
    :param accepted_canonical_sites: list of accepted canonical splice sites
    :param indelInfo: indels near junction information, dict of pbid --> list of junctions near indel (in Interval format)
    :param genome_dict: genome fasta dict
    :param fout: DictWriter handle
    :param covInf: (optional) junction coverage information, dict of (chrom,strand) -> (0-based start,1-based end) -> dict of {sample -> unique read count}
    :param covNames: (optional) list of sample names for the junction coverage information
    :param phyloP_reader: (optional) dict of (chrom,0-based coord) --> phyloP score

    Write a record for each junction in query isoform
    """

    def find_closest_in_list(lst, pos):
        i = bisect.bisect_left(lst, pos)
        if i == 0:
            return lst[0] - pos
        elif i == len(lst):
            return lst[-1] - pos
        else:
            a, b = lst[i - 1] - pos, lst[i] - pos
            if abs(a) < abs(b):
                return a
            else:
                return b

    if trec.chrom not in junctions_by_chr:
        # nothing to do
        return

    # go through each trec junction
    for junction_index, (d, a) in enumerate(trec.junctions):
        # NOTE: donor just means the start, not adjusted for strand
        # find the closest junction start site
        min_diff_s = -find_closest_in_list(junctions_by_chr[trec.chrom]["donors"], d)
        # find the closest junction end site
        min_diff_e = find_closest_in_list(junctions_by_chr[trec.chrom]["acceptors"], a)

        splice_site = trec.get_splice_site(genome_dict, junction_index)

        indel_near_junction = "NA"
        if indelInfo is not None:
            indel_near_junction = (
                "TRUE"
                if (trec.id in indelInfo and Interval(d, a) in indelInfo[trec.id])
                else "FALSE"
            )

        sample_cov = defaultdict(lambda: 0)  # sample -> unique count for this junction
        if covInf is not None:
            sample_cov = covInf[(trec.chrom, trec.strand)][(d, a)]

        # if phyloP score dict exists, give the triplet score of (last base in donor exon), donor site -- similarly for acceptor
        phyloP_start, phyloP_end = "NA", "NA"
        if phyloP_reader is not None:
            phyloP_start = ",".join(
                [
                    str(x)
                    for x in [
                        phyloP_reader.get_pos(trec.chrom, d - 1),
                        phyloP_reader.get_pos(trec.chrom, d),
                        phyloP_reader.get_pos(trec.chrom, d + 1),
                    ]
                ]
            )
            phyloP_end = ",".join(
                [
                    str(x)
                    for x in [
                        phyloP_reader.get_pos(trec.chrom, a - 1),
                        phyloP_reader.get_pos(trec.chrom, a),
                        phyloP_reader.get_pos(trec.chrom, a + 1),
                    ]
                ]
            )

        qj = {
            "isoform": trec.id,
            "junction_number": "junction_" + str(junction_index + 1),
            "chrom": trec.chrom,
            "strand": trec.strand,
            "genomic_start_coord": d + 1,  # write out as 1-based start
            "genomic_end_coord": a,  # already is 1-based end
            "transcript_coord": "?????",  # this is where the exon ends w.r.t to id sequence, ToDo: implement later
            "junction_category": "known"
            if ((d, a) in junctions_by_chr[trec.chrom]["da_pairs"])
            else "novel",
            "start_site_category": "known" if min_diff_s == 0 else "novel",
            "end_site_category": "known" if min_diff_e == 0 else "novel",
            "diff_to_Ref_start_site": min_diff_s,
            "diff_to_Ref_end_site": min_diff_e,
            "bite_junction": "TRUE"
            if (
                (min_diff_s < 0 or min_diff_e < 0)
                and not (min_diff_s > 0 or min_diff_e > 0)
            )
            else "FALSE",
            "splice_site": splice_site,
            "canonical": "canonical"
            if splice_site in accepted_canonical_sites
            else "non_canonical",
            "RTS_junction": "????",  # First write ???? in _tmp, later is TRUE/FALSE
            "indel_near_junct": indel_near_junction,
            "phyloP_start": phyloP_start,
            "phyloP_end": phyloP_end,
            "sample_with_cov": sum(cov != 0 for cov in sample_cov.values())
            if covInf is not None
            else "NA",
            "total_coverage": sum(sample_cov.values()) if covInf is not None else "NA",
        }

        if covInf is not None:
            for sample in covNames:
                qj[sample] = sample_cov[sample]

        fout.writerow(qj)


def isoformClassification(
    args,
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

    if args.coverage is not None:
        print("**** Reading Splice Junctions coverage files.", file=sys.stdout)
        SJcovNames, SJcovInfo = STARcov_parser(args.coverage)
        fields_junc_cur = FIELDS_JUNC + SJcovNames  # add the samples to the header
    else:
        SJcovNames, SJcovInfo = None, None
        print("Splice Junction Coverage files not provided.", file=sys.stdout)
        fields_junc_cur = FIELDS_JUNC

    if args.cage_peak is not None:
        print("**** Reading CAGE Peak data.", file=sys.stdout)
        cage_peak_obj = CAGEPeak(args.cage_peak)
    else:
        cage_peak_obj = None

    if args.polyA_peak is not None:
        print("**** Reading polyA Peak data.", file=sys.stdout)
        polya_peak_obj = PolyAPeak(args.polyA_peak)
    else:
        polya_peak_obj = None

    if args.polyA_motif_list is not None:
        print("**** Reading PolyA motif list.", file=sys.stdout)
        polyA_motif_list = []
        for line in open(args.polyA_motif_list):
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

    if args.phyloP_bed is not None:
        print("**** Reading PhyloP BED file.", file=sys.stdout)
        phyloP_reader = LazyBEDPointReader(args.phyloP_bed)
    else:
        phyloP_reader = None

    # running classification
    print("**** Performing Classification of Isoforms....", file=sys.stdout)

    accepted_canonical_sites = list(args.sites.split(","))

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
                nPolyA=args.window,
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
                    args.novel_gene_prefix is not None
                ):  # used by splits to not have redundant novelGene IDs
                    isoform_hit.genes = [
                        "novelGene_"
                        + str(args.novel_gene_prefix)
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


def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    mean = sum(data) * 1.0 / n  # mean
    var = sum(pow(x - mean, 2) for x in data) / n  # variance
    return math.sqrt(var)  # standard deviation


def find_polyA_motif(genome_seq, polyA_motif_list):
    """
    :param genome_seq: genomic sequence to search polyA motifs from, must already be oriented
    :param polyA_motif_list: ranked list of motifs to find, report the top one found
    :return: polyA_motif, polyA_dist (how many bases upstream is this found)
    """
    for motif in polyA_motif_list:
        i = genome_seq.find(motif)
        if i >= 0:
            return motif, -(len(genome_seq) - i - len(motif) + 1)
    return "NA", "NA"


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
    samples = ["NA"]
    flag_single_sample = True

    f = open(fl_count_filename)
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith("#"):
            # if it first thing is superPBID or id or pbid
            if line.startswith("pbid"):
                type = "SINGLE_SAMPLE"
                sep = "\t"
            elif line.startswith("superPBID"):
                type = "MULTI_CHAIN"
                sep = "\t"
            elif line.startswith("id"):
                type = "MULTI_DEMUX"
                sep = ","
            else:
                raise Exception("Unexpected count file format! Abort!")
            f.seek(cur_pos)
            break

    reader = DictReader(f, delimiter=sep)
    count_header = reader.fieldnames
    if type == "SINGLE_SAMPLE":
        if "count_fl" not in count_header:
            print(
                f"Expected `count_fl` field in count file {fl_count_filename}. Abort!",
                file=sys.stderr,
            )
            sys.exit(-1)
        d = {r["pbid"]: r for r in reader}
    elif type == "MULTI_CHAIN":
        d = {r["superPBID"]: r for r in reader}
        flag_single_sample = False
    elif type == "MULTI_DEMUX":
        d = {r["id"]: r for r in reader}
        flag_single_sample = False
    else:
        print(
            f"Expected pbid or superPBID as a column in count file {fl_count_filename}. Abort!",
            file=sys.stderr,
        )
        sys.exit(-1)
    f.close()

    if flag_single_sample:  # single sample
        for k, v in d.items():
            fl_count_dict[k] = int(v["count_fl"])
    else:  # multi-sample
        for k, v in d.items():
            fl_count_dict[k] = {}
            samples = list(v.keys())
            for sample, count in v.items():
                if sample not in ("superPBID", "id"):
                    fl_count_dict[k][sample] = int(count) if count != "NA" else 0

    samples.sort()

    if type == "MULTI_CHAIN":
        samples.remove("superPBID")
    elif type == "MULTI_DEMUX":
        samples.remove("id")

    return samples, fl_count_dict


def run(args):
    global outputClassPath
    global outputJuncPath

    outputClassPath, outputJuncPath = get_class_junc_filenames(args)

    start3 = timeit.default_timer()

    print("**** Parsing provided files....", file=sys.stdout)
    print(f"Reading genome fasta {args.genome}....", file=sys.stdout)
    # NOTE: can't use LazyFastaReader because inefficient. Bring the whole genome in!
    genome_dict = {r.name: r for r in SeqIO.parse(open(args.genome), "fasta")}

    # correction of sequences and ORF prediction (if gtf provided instead of fasta file, correction of sequences will be skipped)
    orfDict = correctionPlusORFpred(args, genome_dict)

    # parse reference id (GTF) to dicts
    (
        refs_1exon_by_chr,
        refs_exons_by_chr,
        junctions_by_chr,
        junctions_by_gene,
        start_ends_by_gene,
    ) = reference_parser(args, list(genome_dict.keys()))

    # parse query isoforms
    isoforms_by_chr = isoforms_parser(args)

    # Run indel computation if sam exists
    # indelsJunc: dict of pbid --> list of junctions near indel (in Interval format)
    # indelsTotal: dict of pbid --> total indels count
    if os.path.exists(corrSAM):
        (indelsJunc, indelsTotal) = calc_indels_from_sam(corrSAM)
    else:
        indelsJunc = None
        indelsTotal = None

    # isoform classification + intra-priming + id and junction characterization
    isoforms_info = isoformClassification(
        args,
        isoforms_by_chr,
        refs_1exon_by_chr,
        refs_exons_by_chr,
        junctions_by_chr,
        junctions_by_gene,
        start_ends_by_gene,
        genome_dict,
        indelsJunc,
        orfDict,
    )

    print(f"Number of classified isoforms: {len(isoforms_info)}", file=sys.stdout)

    write_collapsed_GFF_with_CDS(isoforms_info, corrGTF, corrGTF + ".cds.gff")
    # os.rename(corrGTF+'.cds.gff', corrGTF)

    # RT-switching computation
    print("**** RT-switching computation....", file=sys.stderr)

    # RTS_info: dict of (pbid) -> list of RT junction. if RTS_info[pbid] == [], means all junctions are non-RT.
    RTS_info = rts([outputJuncPath + "_tmp", args.genome, "-a"], genome_dict)
    for pbid in isoforms_info:
        if pbid in RTS_info and len(RTS_info[pbid]) > 0:
            isoforms_info[pbid].RT_switching = "TRUE"
        else:
            isoforms_info[pbid].RT_switching = "FALSE"

    # FSM classification
    geneFSM_dict = defaultdict(lambda: [])
    for iso in isoforms_info:
        gene = isoforms_info[
            iso
        ].geneName()  # if multi-gene, returns "geneA_geneB_geneC..."
        geneFSM_dict[gene].append(isoforms_info[iso].str_class)

    fields_class_cur = FIELDS_CLASS
    # FL count file
    if args.fl_count:
        if not os.path.exists(args.fl_count):
            print(f"FL count file {args.fl_count} does not exist!", file=sys.stderr)
            sys.exit(-1)
        print("**** Reading Full-length read abundance files...", file=sys.stderr)
        fl_samples, fl_count_dict = FLcount_parser(args.fl_count)
        for pbid in fl_count_dict:
            if pbid not in isoforms_info:
                print(
                    f"WARNING: {pbid} found in FL count file but not in input fasta.",
                    file=sys.stderr,
                )
        if len(fl_samples) == 1:  # single sample from PacBio
            print("Single-sample PacBio FL count format detected.", file=sys.stderr)
            for iso in isoforms_info:
                if iso in fl_count_dict:
                    isoforms_info[iso].FL = fl_count_dict[iso]
                else:
                    print(
                        f"WARNING: {iso} not found in FL count file. Assign count as 0.",
                        file=sys.stderr,
                    )
                    isoforms_info[iso].FL = 0
        else:  # multi-sample
            print("Multi-sample PacBio FL count format detected.", file=sys.stderr)
            fields_class_cur = FIELDS_CLASS + ["FL." + s for s in fl_samples]
            for iso in isoforms_info:
                if iso in fl_count_dict:
                    isoforms_info[iso].FL_dict = fl_count_dict[iso]
                else:
                    print(
                        f"WARNING: {iso} not found in FL count file. Assign count as 0.",
                        file=sys.stderr,
                    )
                    isoforms_info[iso].FL_dict = defaultdict(lambda: 0)
    else:
        print("Full-length read abundance files not provided.", file=sys.stderr)

    # Isoform expression information
    if args.expression:
        print("**** Reading Isoform Expression Information.", file=sys.stderr)
        exp_dict = expression_parser(args.expression)
        gene_exp_dict = {}
        for iso in isoforms_info:
            if iso not in exp_dict:
                exp_dict[iso] = 0
                print(
                    f"WARNING: isoform {iso} not found in expression matrix. Assigning TPM of 0.",
                    file=sys.stderr,
                )
            gene = isoforms_info[iso].geneName()
            if gene not in gene_exp_dict:
                gene_exp_dict[gene] = exp_dict[iso]
            else:
                gene_exp_dict[gene] = gene_exp_dict[gene] + exp_dict[iso]
    else:
        exp_dict = None
        gene_exp_dict = None
        print("Isoforms expression files not provided.", file=sys.stderr)

    # Adding indel, FSM class and expression information
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()
        if exp_dict is not None and gene_exp_dict is not None:
            isoforms_info[iso].geneExp = gene_exp_dict[gene]
            isoforms_info[iso].isoExp = exp_dict[iso]
        if len(geneFSM_dict[gene]) == 1:
            isoforms_info[iso].FSM_class = "A"
        elif "full-splice_match" in geneFSM_dict[gene]:
            isoforms_info[iso].FSM_class = "C"
        else:
            isoforms_info[iso].FSM_class = "B"

    if indelsTotal is not None:
        for iso in isoforms_info:
            if iso in indelsTotal:
                isoforms_info[iso].nIndels = indelsTotal[iso]
            else:
                isoforms_info[iso].nIndels = 0

    # Read junction files and create attributes per id
    # Read the junction information to fill in several remaining unfilled fields in classification
    # (1) "canonical": is "canonical" if all junctions are canonical, otherwise "non_canonical"
    # (2) "bite": is TRUE if any of the junction "bite_junction" field is TRUE

    reader = DictReader(open(outputJuncPath + "_tmp"), delimiter="\t")
    fields_junc_cur = reader.fieldnames

    sj_covs_by_isoform = defaultdict(
        lambda: []
    )  # pbid --> list of total_cov for each junction so we can calculate SD later
    for r in reader:
        # only need to do assignment if:
        # (1) the .canonical field is still "NA"
        # (2) the junction is non-canonical
        assert r["canonical"] in ("canonical", "non_canonical")
        if (isoforms_info[r["isoform"]].canonical == "NA") or (
            r["canonical"] == "non_canonical"
        ):
            isoforms_info[r["isoform"]].canonical = r["canonical"]

        if (isoforms_info[r["isoform"]].bite == "NA") or (r["bite_junction"] == "TRUE"):
            isoforms_info[r["isoform"]].bite = r["bite_junction"]

        if r["indel_near_junct"] == "TRUE":
            if isoforms_info[r["isoform"]].nIndelsJunc == "NA":
                isoforms_info[r["isoform"]].nIndelsJunc = 0
            isoforms_info[r["isoform"]].nIndelsJunc += 1

        # min_cov: min( total_cov[j] for each junction j in this isoform )
        # min_cov_pos: the junction [j] that attributed to argmin(total_cov[j])
        # min_sample_cov: min( sample_cov[j] for each junction in this isoform )
        # sd_cov: sd( total_cov[j] for each junction j in this isoform )
        if r["sample_with_cov"] != "NA":
            sample_with_cov = int(r["sample_with_cov"])
            if (isoforms_info[r["isoform"]].min_samp_cov == "NA") or (
                isoforms_info[r["isoform"]].min_samp_cov > sample_with_cov
            ):
                isoforms_info[r["isoform"]].min_samp_cov = sample_with_cov

        if r["total_coverage"] != "NA":
            total_cov = int(r["total_coverage"])
            sj_covs_by_isoform[r["isoform"]].append(total_cov)
            if (isoforms_info[r["isoform"]].min_cov == "NA") or (
                isoforms_info[r["isoform"]].min_cov > total_cov
            ):
                isoforms_info[r["isoform"]].min_cov = total_cov
                isoforms_info[r["isoform"]].min_cov_pos = r["junction_number"]

    for pbid, covs in sj_covs_by_isoform.items():
        isoforms_info[pbid].sd = pstdev(covs)

    # Printing output file:
    print("**** Writing output files....", file=sys.stderr)

    # sort isoform keys
    iso_keys = list(isoforms_info.keys())
    iso_keys.sort(key=lambda x: (isoforms_info[x].chrom, isoforms_info[x].id))
    with open(outputClassPath, "w") as h:
        fout_class = DictWriter(h, fieldnames=fields_class_cur, delimiter="\t")
        fout_class.writeheader()
        for iso_key in iso_keys:
            fout_class.writerow(isoforms_info[iso_key].as_dict())

    # Now that RTS info is obtained, we can write the final junctions.txt
    with open(outputJuncPath, "w") as h:
        fout_junc = DictWriter(h, fieldnames=fields_junc_cur, delimiter="\t")
        fout_junc.writeheader()
        for r in DictReader(open(outputJuncPath + "_tmp"), delimiter="\t"):
            if r["isoform"] in RTS_info:
                if r["junction_number"] in RTS_info[r["isoform"]]:
                    r["RTS_junction"] = "TRUE"
                else:
                    r["RTS_junction"] = "FALSE"
            fout_junc.writerow(r)

    # Generating report
    if not args.skip_report:
        print("**** Generating SQANTI3 report....", file=sys.stderr)
        cmd = (
            RSCRIPTPATH
            + f" {utilitiesPath}/{RSCRIPT_REPORT} {outputClassPath} {outputJuncPath} {args.doc} {utilitiesPath}"
        )
        if subprocess.check_call(cmd, shell=True) != 0:
            print(f"ERROR running command: {cmd}", file=sys.stderr)
            sys.exit(-1)
    stop3 = timeit.default_timer()

    print("Removing temporary files....", file=sys.stderr)
    os.remove(outputClassPath + "_tmp")
    os.remove(outputJuncPath + "_tmp")

    print(f"SQANTI3 complete in {stop3 - start3} sec.", file=sys.stderr)

    # IsoAnnot Lite implementation
    ISOANNOT_PROG = os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")
    if args.isoAnnotLite:
        if args.gff3:
            ISOANNOT_CMD = (
                "python "
                + ISOANNOT_PROG
                + f" {corrGTF} {outputClassPath} {outputJuncPath} -gff3 {args.gff3} -d {args.dir} -o {args.output}"
            )
        else:
            ISOANNOT_CMD = (
                "python "
                + ISOANNOT_PROG
                + f" {corrGTF} {outputClassPath} {outputJuncPath} -d {args.dir} -o {args.output}"
            )
        if subprocess.check_call(ISOANNOT_CMD, shell=True) != 0:
            print(f"ERROR running command: {ISOANNOT_CMD}", file=sys.stderr)
            sys.exit(-1)


def rename_isoform_seqids(input_fasta, force_id_ignore=False):
    """
    Rename input isoform fasta/fastq, which is usually mapped, collapsed Iso-Seq data with IDs like:

    PB.1.1|chr1:10-100|xxxxxx

    to just being "PB.1.1"

    :param input_fasta: Could be either fasta or fastq, autodetect.
    :return: output fasta with the cleaned up sequence ID, is_fusion flag
    """
    type = "fasta"
    with open(input_fasta) as h:
        if h.readline().startswith("@"):
            type = "fastq"
    f = open(input_fasta[: input_fasta.rfind(".")] + ".renamed.fasta", "w")
    for r in SeqIO.parse(open(input_fasta), type):
        m1 = seqid_rex1.match(r.id)
        m2 = seqid_rex2.match(r.id)
        m3 = seqid_fusion.match(r.id)
        if not force_id_ignore and (m1 is None and m2 is None and m3 is None):
            print(
                f"Invalid input IDs! Expected PB.X.Y or PB.X.Y|xxxxx or PBfusion.X format but saw {r.id} instead. Abort!",
                file=sys.stderr,
            )
            sys.exit(-1)
        if r.id.startswith("PB.") or r.id.startswith(
            "PBfusion."
        ):  # PacBio fasta header
            newid = r.id.split("|")[0]
        else:
            raw = r.id.split("|")
            if len(raw) > 4:  # RefSeq fasta header
                newid = raw[3]
            else:
                newid = r.id.split()[0]  # Ensembl fasta header
        f.write(f">{newid}\n{r.seq}\n")
    f.close()
    return f.name


class CAGEPeak:
    def __init__(self, cage_bed_filename):
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(
            lambda: IntervalTree()
        )  # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            tss0 = int(raw[6])
            self.cage_peaks[(chrom, strand)].insert(start0, end1, (tss0, start0, end1))

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param start0: 0-based start of the 5' end to query
        :return: <True/False falls within a cage peak>, <nearest dist to TSS>
        dist to TSS is 0 if right on spot
        dist to TSS is + if downstream, - if upstream (watch for strand!!!)
        """
        within_peak, dist_peak = False, "NA"
        for (tss0, start0, end1) in self.cage_peaks[(chrom, strand)].find(
            query - search_window, query + search_window
        ):
            # Skip those cage peaks that are downstream the detected TSS because degradation just make the transcript shorter
            if strand == "+" and start0 > int(query) and end1 > int(query):
                continue
            if strand == "-" and start0 < int(query) and end1 < int(query):
                continue
            ##
            if not within_peak:
                within_peak, dist_peak = (
                    (start0 <= query < end1),
                    (query - tss0) * (-1 if strand == "-" else +1),
                )
            else:
                d = (query - tss0) * (-1 if strand == "-" else +1)
                if abs(d) < abs(dist_peak):
                    within_peak, dist_peak = (start0 <= query < end1), d
        return within_peak, dist_peak


class PolyAPeak:
    def __init__(self, polya_bed_filename):
        self.polya_bed_filename = polya_bed_filename
        self.polya_peaks = defaultdict(
            lambda: IntervalTree()
        )  # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.polya_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            self.polya_peaks[(chrom, strand)].insert(start0, end1, (start0, end1))

    def find(self, chrom, strand, query, search_window=100):
        """
        :param start0: 0-based start of the 5' end to query
        :return: <True/False falls within some distance to polyA>, distance to closest
        + if downstream, - if upstream (watch for strand!!!)
        """
        assert strand in ("+", "-")
        hits = self.polya_peaks[(chrom, strand)].find(
            query - search_window, query + search_window
        )
        if len(hits) == 0:
            return False, None
        else:
            s0, e1 = hits[0]
            min_dist = query - s0
            for s0, e1 in hits[1:]:
                d = query - s0
                if abs(d) < abs(min_dist):
                    min_dist = d
            if strand == "-":
                min_dist = -min_dist
            return True, min_dist


def split_input_run(args):
    if os.path.exists(SPLIT_ROOT_DIR):
        print(
            f"WARNING: {SPLIT_ROOT_DIR} directory already exists! Abort!",
            file=sys.stderr,
        )
        sys.exit(-1)
    else:
        os.makedirs(SPLIT_ROOT_DIR)

    if args.gtf:
        recs = [r for r in collapseGFFReader(args.isoforms)]
        n = len(recs)
        chunk_size = n // args.chunks + (n % args.chunks > 0)
        split_outs = []
        # pdb.set_trace()
        for i in range(args.chunks):
            if i * chunk_size >= n:
                break
            d = os.path.join(SPLIT_ROOT_DIR, str(i))
            os.makedirs(d)
            f = open(
                os.path.join(d, os.path.basename(args.isoforms) + ".split" + str(i)),
                "w",
            )
            for j in range(i * chunk_size, min((i + 1) * chunk_size, n)):
                write_collapseGFF_format(f, recs[j])
            f.close()
            split_outs.append((os.path.abspath(d), f.name))
    else:
        recs = [r for r in SeqIO.parse(open(args.isoforms), "fasta")]
        n = len(recs)
        chunk_size = n // args.chunks + (n % args.chunks > 0)
        split_outs = []
        for i in range(args.chunks):
            if i * chunk_size >= n:
                break
            d = os.path.join(SPLIT_ROOT_DIR, str(i))
            os.makedirs(d)
            f = open(
                os.path.join(d, os.path.basename(args.isoforms) + ".split" + str(i)),
                "w",
            )
            for j in range(i * chunk_size, min((i + 1) * chunk_size, n)):
                SeqIO.write(recs[j], f, "fasta")
            f.close()
            split_outs.append((os.path.abspath(d), f.name))

    pools = []
    for i, (d, x) in enumerate(split_outs):
        print(f"launching worker on on {x}....")
        args2 = copy.deepcopy(args)
        args2.isoforms = x
        args2.novel_gene_prefix = str(i)
        args2.dir = d
        args2.skip_report = True
        p = Process(target=run, args=(args2,))
        p.start()
        pools.append(p)

    for p in pools:
        p.join()
    return [d for (d, x) in split_outs]


def combine_split_runs(args, split_dirs):
    """
    Combine .faa, .fasta, .gtf, .classification.txt, .junctions.txt
    Then write out the PDF report
    """
    corrGTF, corrSAM, corrFASTA, corrORF = get_corr_filenames(args)
    outputClassPath, outputJuncPath = get_class_junc_filenames(args)

    if not args.skipORF:
        f_faa = open(corrORF, "w")
    f_fasta = open(corrFASTA, "w")
    f_gtf = open(corrGTF, "w")
    f_class = open(outputClassPath, "w")
    f_junc = open(outputJuncPath, "w")

    for i, split_d in enumerate(split_dirs):
        _gtf, _sam, _fasta, _orf = get_corr_filenames(args, split_d)
        _class, _junc = get_class_junc_filenames(args, split_d)
        if not args.skipORF:
            with open(_orf) as h:
                f_faa.write(h.read())
        with open(_gtf) as h:
            f_gtf.write(h.read())
        with open(_fasta) as h:
            f_fasta.write(h.read())
        with open(_class) as h:
            if i == 0:
                f_class.write(h.readline())
            else:
                h.readline()
            f_class.write(h.read())
        with open(_junc) as h:
            if i == 0:
                f_junc.write(h.readline())
            else:
                h.readline()
            f_junc.write(h.read())

    f_fasta.close()
    f_gtf.close()
    f_class.close()
    f_junc.close()
    if not args.skipORF:
        f_faa.close()

    if not args.skip_report:
        print("**** Generating SQANTI3 report....", file=sys.stderr)
        cmd = (
            RSCRIPTPATH
            + f" {utilitiesPath}/{RSCRIPT_REPORT} {outputClassPath} {outputJuncPath} {args.doc} {utilitiesPath}"
        )
        if subprocess.check_call(cmd, shell=True) != 0:
            print(f"ERROR running command: {cmd}", file=sys.stderr)
            sys.exit(-1)


def main():
    global utilitiesPath

    # arguments
    parser = argparse.ArgumentParser(
        description="Structural and Quality Annotation of Novel Transcript Isoforms"
    )
    parser.add_argument(
        "isoforms",
        help='\tIsoforms (FASTA/FASTQ or gtf format; By default "FASTA/FASTQ". For GTF, use --gtf',
    )
    parser.add_argument("annotation", help="\t\tReference annotation file (GTF format)")
    parser.add_argument("genome", help="\t\tReference genome (Fasta format)")
    parser.add_argument(
        "--min_ref_len",
        type=int,
        default=200,
        help="\t\tMinimum reference transcript length (default: 200 bp)",
    )
    parser.add_argument(
        "--force_id_ignore",
        action="store_true",
        default=False,
        help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)",
    )
    parser.add_argument(
        "--aligner_choice", choices=["minimap2", "deSALT", "gmap"], default="minimap2"
    )
    parser.add_argument(
        "--cage_peak", help="\t\tFANTOM5 Cage Peak (BED format, optional)"
    )
    parser.add_argument(
        "--polyA_motif_list", help="\t\tRanked list of polyA motifs (text, optional)"
    )
    parser.add_argument("--polyA_peak", help="\t\tPolyA Peak (BED format, optional)")
    parser.add_argument(
        "--phyloP_bed", help="\t\tPhyloP BED for conservation score (BED, optional)"
    )
    parser.add_argument(
        "--skipORF",
        default=False,
        action="store_true",
        help="\t\tSkip ORF prediction (to save time)",
    )
    parser.add_argument(
        "--is_fusion",
        default=False,
        action="store_true",
        help="\t\tInput are fusion isoforms, must supply GTF as input using --gtf",
    )
    parser.add_argument(
        "-g",
        "--gtf",
        help="\t\tUse when running SQANTI by using as input a gtf of isoforms",
        action="store_true",
    )
    parser.add_argument(
        "-e",
        "--expression",
        help="\t\tExpression matrix (supported: Kallisto tsv)",
        required=False,
    )
    parser.add_argument(
        "-x",
        "--gmap_index",
        help="\t\tPath and prefix of the reference index created by gmap_build. Mandatory if using GMAP unless -g option is specified.",
    )
    parser.add_argument(
        "-t",
        "--cpus",
        default=10,
        type=int,
        help="\t\tNumber of threads used during alignment by aligners. (default: 10)",
    )
    parser.add_argument(
        "-n",
        "--chunks",
        default=1,
        type=int,
        help="\t\tNumber of chunks to split SQANTI3 analysis in for speed up (default: 1).",
    )
    # parser.add_argument('-z', '--sense', help='\t\tOption that helps aligners know that the exons in you cDNA sequences are in the correct sense. Applicable just when you have a high quality set of cDNA sequences', required=False, action='store_true')
    parser.add_argument(
        "-o", "--output", help="\t\tPrefix for output files.", required=False
    )
    parser.add_argument(
        "-d",
        "--dir",
        help="\t\tDirectory for output files. Default: Directory where the script was run.",
        required=False,
    )
    parser.add_argument(
        "-c",
        "--coverage",
        help='\t\tJunction coverage files (provide a single file or a file pattern, ex: "mydir/*.junctions").',
        required=False,
    )
    parser.add_argument(
        "-s",
        "--sites",
        default="ATAC,GCAG,GTAG",
        help="\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.",
        required=False,
    )
    parser.add_argument(
        "-w",
        "--window",
        default="20",
        help="\t\tSize of the window in the genomic DNA screened for Adenine content downstream of TTS",
        required=False,
        type=int,
    )
    parser.add_argument(
        "--genename",
        help="\t\tUse gene_name tag from GTF to define genes. Default: gene_id used to define genes",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-fl",
        "--fl_count",
        help="\t\tFull-length PacBio abundance file",
        required=False,
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Display program version number.",
        action="version",
        version=f"{get_version(root='..', relative_to=__file__)}",
    )
    parser.add_argument(
        "--skip_report", action="store_true", default=False, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--isoAnnotLite",
        help="\t\tRun isoAnnot Lite to output a tappAS-compatible gff3 file",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--gff3",
        help="\t\tPrecomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes",
        required=False,
    )

    args = parser.parse_args()

    if args.is_fusion:
        print(
            "WARNING: Currently if --is_fusion is used, no ORFs will be predicted.",
            file=sys.stderr,
        )
        args.skipORF = True
        if not args.gtf:
            print(
                "ERROR: if --is_fusion is on, must supply GTF as input and use --gtf!",
                file=sys.stderr,
            )
            sys.exit(-1)

    if args.gff3 is not None:
        args.gff3 = os.path.abspath(args.gff3)
        if not os.path.isfile(args.gff3):
            print(
                f"ERROR: Precomputed tappAS GFF3 annoation file {args.genome} doesn't exist. Abort!",
                file=sys.stderr,
            )
            sys.exit(-1)

    if args.expression is not None:
        if os.path.isdir(args.expression) is True:
            print(
                f"Expression files located in {args.expression} folder", file=sys.stderr
            )
        else:
            for f in args.expression.split(","):
                if not os.path.exists(f):
                    print(f"Expression file {f} not found. Abort!", file=sys.stderr)
                    sys.exit(-1)

    # path and prefix for output files
    if args.output is None:
        args.output = os.path.splitext(os.path.basename(args.isoforms))[0]

    if args.dir is None:
        args.dir = os.getcwd()
    else:
        args.dir = os.path.abspath(args.dir)
        if os.path.isdir(args.dir):
            print(
                f"WARNING: output directory {args.dir} already exists. Overwriting!",
                file=sys.stderr,
            )
        else:
            os.makedirs(args.dir)

    args.genome = os.path.abspath(args.genome)
    if not os.path.isfile(args.genome):
        print(
            f"ERROR: genome fasta {args.genome} doesn't exist. Abort!", file=sys.stderr
        )
        sys.exit()

    args.isoforms = os.path.abspath(args.isoforms)
    if not os.path.isfile(args.isoforms):
        print(
            f"ERROR: Input isoforms {args.isoforms} doesn't exist. Abort!",
            file=sys.stderr,
        )
        sys.exit()

    if not args.gtf:
        if args.aligner_choice == "gmap":
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                print(
                    f"GMAP index {args.gmap_index} doesn't exist! Abort.",
                    file=sys.stderr,
                )
                sys.exit()
        elif args.aligner_choice == "deSALT":
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                print(
                    f"deSALT index {args.gmap_index} doesn't exist! Abort.",
                    file=sys.stderr,
                )
                sys.exit()

        print("Cleaning up isoform IDs...", file=sys.stderr)
        args.isoforms = rename_isoform_seqids(args.isoforms, args.force_id_ignore)
        print(
            f"Cleaned up isoform fasta file written to: {args.isoforms}",
            file=sys.stderr,
        )

    args.annotation = os.path.abspath(args.annotation)
    if not os.path.isfile(args.annotation):
        print(
            f"ERROR: Annotation {args.annotation} doesn't exist. Abort!",
            file=sys.stderr,
        )
        sys.exit()

    # if args.aligner_choice == "gmap":
    #    args.sense = "sense_force" if args.sense else "auto"
    # elif args.aligner_choice == "minimap2":
    #    args.sense = "f" if args.sense else "b"
    # (Liz) turned off option for --sense, always TRUE
    if args.aligner_choice == "gmap":
        args.sense = "sense_force"
    elif args.aligner_choice == "minimap2":
        args.sense = "f"
    # elif args.aligner_choice == "deSALT":  #deSALT does not support this yet
    #    args.sense = "--trans-strand"

    args.novel_gene_prefix = None
    # Print out parameters so can be put into report PDF later
    args.doc = os.path.join(os.path.abspath(args.dir), args.output + ".params.txt")
    print(f"Write arguments to {args.doc}...", file=sys.stdout)
    with open(args.doc, "w") as f:
        f.write("Version\t" + __version__ + "\n")
        f.write("Input\t" + os.path.basename(args.isoforms) + "\n")
        f.write("Annotation\t" + os.path.basename(args.annotation) + "\n")
        f.write("Genome\t" + os.path.basename(args.genome) + "\n")
        f.write("Aligner\t" + args.aligner_choice + "\n")
        f.write(
            "FLCount\t"
            + (os.path.basename(args.fl_count) if args.fl_count is not None else "NA")
            + "\n"
        )
        f.write(
            "Expression\t"
            + (
                os.path.basename(args.expression)
                if args.expression is not None
                else "NA"
            )
            + "\n"
        )
        f.write(
            "Junction\t"
            + (os.path.basename(args.coverage) if args.coverage is not None else "NA")
            + "\n"
        )
        f.write(
            "CagePeak\t"
            + (os.path.basename(args.cage_peak) if args.cage_peak is not None else "NA")
            + "\n"
        )
        f.write(
            "PolyA\t"
            + (
                os.path.basename(args.polyA_motif_list)
                if args.polyA_motif_list is not None
                else "NA"
            )
            + "\n"
        )
        f.write(
            "PolyAPeak\t"
            + (
                os.path.basename(args.polyA_peak)
                if args.polyA_peak is not None
                else "NA"
            )
            + "\n"
        )
        f.write("IsFusion\t" + str(args.is_fusion) + "\n")

    # Running functionality
    print("**** Running SQANTI3...", file=sys.stdout)
    if args.chunks == 1:
        run(args)
    else:
        split_dirs = split_input_run(args)
        combine_split_runs(args, split_dirs)
        shutil.rmtree(SPLIT_ROOT_DIR)


if __name__ == "__main__":
    main()
