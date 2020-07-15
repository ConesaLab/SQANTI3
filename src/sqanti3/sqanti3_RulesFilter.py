#!/usr/bin/env python
__author__ = "etseng@pacb.com"

# import argparse
import click
import distutils.spawn
import os
import subprocess
import sys
from argparse import Namespace
from csv import DictReader, DictWriter

from Bio import SeqIO
from cupcake.sequence.BioReaders import GMAPSAMReader
from cupcake.sequence.GFF import collapseGFFReader, write_collapseGFF_format

from typing import Optional
import logging

from sqanti3.__about__ import __version__

"""
Lightweight filtering of SQANTI by using .classification.txt output

Only keep Iso-Seq isoforms if:
The isoform is FSM, ISM, or NIC and (does not have intrapriming or has polyA_motif)
The isoform is NNC, does not have intrapriming/or polyA motif, not RT-switching, and all junctions are either all canonical or short-read-supported
The isoform is antisense, intergenic, genic, does not have intrapriming/or polyA motif, not RT-switching, and all junctions are either all canonical or short-read-supported
"""


utilitiesPath = os.path.dirname(os.path.realpath(__file__)) + "/utilities/"
RSCRIPTPATH = distutils.spawn.find_executable("Rscript")
RSCRIPT_REPORT = "SQANTI3_report.R"

# necessary for use in f-strings since \n and \t do not work
NEWLINE = "\n"
TAB = "\t"

CATEGORY_DICT = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
    "antisense": "AS",
    "intergenic": "intergenic",
    "genic_intron": "intron",
    "genic": "genic",
    "fusion": "fusion",
}


def sqanti_filter_lite(
    sqanti_class: str,
    isoforms: str,
    annotation: str,
    sam: Optional[str],
    faa: Optional[str],
    intrapriming: float,
    runAlength: int,
    max_dist_to_known_end: int,
    min_cov: int,
    filter_mono_exonic: bool,
    skipGTF: bool,
    skipFaFq: bool,
    skipJunction: bool,
) -> None:

    logger = logging.getLogger(__name__)

    fafq_type = "fasta"
    with open(isoforms) as h:
        if h.readline().startswith("@"):
            fafq_type = "fastq"

    prefix = sqanti_class[: sqanti_class.rfind(".")]

    with open(f"{prefix}.filtered_lite_reasons.txt", "w") as fcsv:
        header = (
            f"# classification: {sqanti_class}{NEWLINE}"
            f"# isoform: {isoforms}{NEWLINE}"
            f"# intrapriming cutoff: {intrapriming}{NEWLINE}"
            f"# min_cov cutoff: {min_cov}{NEWLINE}"
            f"filtered_isoform,reason{NEWLINE}"
        )
        fcsv.write(header)

        seqids_to_keep = set()
        total_count = 0
        for r in DictReader(open(sqanti_class), delimiter="\t"):
            total_count += 1
            filter_flag, filter_msg = False, ""
            percA = float(r["perc_A_downstream_TTS"]) / 100
            assert 0 <= percA <= 1
            runA = 0
            while runA < len(r["seq_A_downstream_TTS"]):
                if r["seq_A_downstream_TTS"][runA] != "A":
                    break
                runA += 1
            calc_min_cov = float(r["min_cov"]) if r["min_cov"] != "NA" else None
            num_exon = int(r["exons"])
            is_RTS = r["RTS_stage"] == "TRUE"
            is_canonical = r["all_canonical"] == "canonical"
            is_monoexonic = num_exon == 1

            cat = CATEGORY_DICT[r["structural_category"]]

            potential_intrapriming = (
                (percA >= intrapriming or runA >= runAlength)
                and r["polyA_motif"] == "NA"
                and (
                    r["diff_to_gene_TSS"] == "NA"
                    or abs(int(r["diff_to_gene_TTS"])) > max_dist_to_known_end
                )
            )

            if cat in ["FSM"]:
                if potential_intrapriming:
                    filter_flag, filter_msg = True, "IntraPriming"
                elif filter_mono_exonic and is_monoexonic:
                    filter_flag, filter_msg = True, "Mono-Exonic"
            else:
                if potential_intrapriming:
                    filter_flag, filter_msg = True, "IntraPriming"
                elif filter_mono_exonic and is_monoexonic:
                    filter_flag, filter_msg = True, "Mono-Exonic"
                elif is_RTS:
                    filter_flag, filter_msg = True, "RTSwitching"
                elif (not is_canonical) and (
                    calc_min_cov is None
                    or (calc_min_cov is not None and calc_min_cov < min_cov)
                ):
                    filter_flag, filter_msg = True, "LowCoverage/Non-Canonical"

            if not filter_flag:
                seqids_to_keep.add(r["isoform"])
            else:
                fcsv.write(f"{r['isoform']},{filter_msg}{NEWLINE}")

    logger.info(
        f"{total_count} isoforms read from {sqanti_class}. {len(seqids_to_keep)} to be kept."
    )

    with open(f"{prefix}.filtered_lite.{fafq_type}", "w") as fout:
        if not skipFaFq:
            for r in SeqIO.parse(open(isoforms), fafq_type):
                if r.id in seqids_to_keep:
                    SeqIO.write(r, fout, fafq_type)
            logger.info(f"Output written to: {fout.name}")

    # write out a new .classification.txt, .junctions.txt
    outputClassPath = f"{prefix}.filtered_lite_classification.txt"
    with open(outputClassPath, "w") as f:
        reader = DictReader(open(sqanti_class), delimiter="\t")
        writer = DictWriter(f, reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for r in reader:
            if r['isoform'] in seqids_to_keep:
                writer.writerow(r)
        logger.info(f"Output written to: {f.name}")

    if not skipJunction:
        outputJuncPath = f"{prefix}.filtered_lite_junctions.txt"
        with open(outputJuncPath, "w") as f:
            reader = DictReader(
                open(sqanti_class.replace("_classification", "_junctions")),
                delimiter="\t",
            )
            writer = DictWriter(f, reader.fieldnames, delimiter="\t")
            writer.writeheader()
            for r in reader:
                if r['isoform'] in seqids_to_keep:
                    writer.writerow(r)
            logger.info(f"Output written to: {f.name}")

    if not skipGTF:
        outputGTF = f"{prefix}.filtered_lite.gtf"
        with open(outputGTF, "w") as f:
            for r in collapseGFFReader(annotation):
                if r.seqid in seqids_to_keep:
                    write_collapseGFF_format(f, r)
            logger.info(f"Output written to: {f.name}")

    if sam is not None:
        outputSam = f"{prefix}.filtered_lite.sam"
        with open(outputSam, "w") as f:
            reader = GMAPSAMReader(sam, True)
            f.write(reader.header)
            for r in reader:
                if r.qID in seqids_to_keep:
                    f.write(f"{r.record_line}{NEWLINE}")
            logger.info(f"Output written to: {f.name}")

    if faa is not None:
        outputFAA = f"{prefix}.filtered_lite.faa"
        with open(outputFAA, "w") as f:
            for r in SeqIO.parse(open(faa), "fasta"):
                if r.id in seqids_to_keep:
                    f.write(f">{r.description}{NEWLINE}{r.seq}{NEWLINE}")
        logger.info(f"Output written to: {f.name}")

    logger.info("Generating SQANTI3 report...")
    cmd = (
        f"{RSCRIPTPATH} {utilitiesPath}/{RSCRIPT_REPORT} {outputClassPath} {outputJuncPath} {'mock'} {utilitiesPath}"
    )
    if subprocess.check_call(cmd, shell=True) != 0:
        logger.error(f"Running command failed: {cmd}")
        sys.exit(-1)


@click.command()
@click.argument("sqanti_class")
@click.argument("isoforms")
@click.argument("annotation")
@click.option(
    "--sam",
    help="(Optional) SAM alignment of the input fasta/fastq",
    type=str,
    default=None,
)
@click.option(
    "--faa", help="(Optional) ORF prediction faa file to be filtered by SQANTI3"
)
@click.option(
    "-a",
    "--intrapriming",
    help="Adenine percentage at genomic 3' end to flag an isoform as intra-priming",
    type=click.FloatRange(min=0.25, max=1.0),
    default=0.6,
    show_default=True,
)
@click.option(
    "-r",
    "--runAlength",
    help="Continuous run-A length at genomic 3' end to flag an isoform as intra-priming",
    type=click.IntRange(min=4, max=21),
    default=6,
    show_default=True,
)
@click.option(
    "-m",
    "--max_dist_to_known_end",
    help="Maximum distance to an annotated 3' end to preserve as a valid 3' end and not filter out",
    type=int,
    default=50,
    show_default=True,
)
@click.option(
    "-c",
    "--min_cov",
    help="Minimum junction coverage for each isoform (only used if min_cov field is not 'NA')",
    type=int,
    default=3,
    show_default=True,
)
@click.option(
    "--filter_mono_exonic",
    help="Filter out all mono-exonic transcripts",
    type=bool,
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--skipGTF",
    help="Skip output of GTF",
    type=bool,
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--skipFaFq",
    help="Skip output of isoform fasta/fastq",
    type=bool,
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--skipJunction",
    help="Skip output of junctions file",
    type=bool,
    default=False,
    show_default=True,
    is_flag=True,
)
@click.version_option()
@click.help_option(show_default=False)
def main(
    sqanti_class: str,
    isoforms: str,
    annotation: str,
    sam: Optional[str] = None,
    faa: Optional[str] = None,
    intrapriming: float = 0.6,
    runalength: int = 6,
    max_dist_to_known_end: int = 50,
    min_cov: int = 3,
    filter_mono_exonic: bool = False,
    skipgtf: bool = False,
    skipfafq: bool = False,
    skipjunction: bool = False,
) -> None:
    """"Filtering of Isoforms based on SQANTI3 attributes
    
    \b
    Parameters:
    -----------
    sqanti_class:
        SQANTI classification output file
    isoforms:
        fasta/fastq isoform file to be filtered by SQANTI3
    annotation:
        GTF matching the input fasta/fastq
    """
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(filename="sqanti3_rulesfilter.log", mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    st = logging.StreamHandler()
    st.setLevel(logging.INFO)
    st.setFormatter(formatter)
    logger.addHandler(st)

    if RSCRIPTPATH is None or os.system(f"{RSCRIPTPATH} --version") != 0:
        logger.error("Rscript executable not found! Abort!")
        sys.exit(-1)

    if intrapriming < 0.25 or intrapriming > 1.0:
        logger.error(
            f"--intrapriming must be between 0.25-1, instead given {intrapriming}! Abort!"
        )
        sys.exit(-1)

    if runalength < 4 or runalength > 20:
        logger.error(
            f"--runAlength must be between 4-20, instead given {runalength}! Abort!"
        )
        sys.exit(-1)

    sqanti_class = os.path.abspath(sqanti_class)
    if not os.path.isfile(sqanti_class):
        logger.error(f"{sqanti_class} doesn't exist. Abort!")
        sys.exit(-1)

    if not os.path.exists(isoforms):
        logger.error(f"{isoforms} doesn't exist. Abort!")
        sys.exit(-1)

    if not os.path.exists(annotation):
        logger.error(f"{annotation} doesn't exist. Abort!")
        sys.exit(-1)

    if sam is not None and not os.path.exists(sam):
        logger.error(f"{sam} doesn't exist. Abort!")
        sys.exit(-1)

    if faa is not None and not os.path.exists(faa):
        logger.error(f"{faa} doesn't exist. Abort!")
        sys.exit(-1)

    logger.info("Running SQANTI2 filtering...")

    sqanti_filter_lite(
        sqanti_class=sqanti_class,
        isoforms=isoforms,
        annotation=annotation,
        sam=sam,
        faa=faa,
        intrapriming=intrapriming,
        runAlength=runalength,
        max_dist_to_known_end=max_dist_to_known_end,
        min_cov=min_cov,
        filter_mono_exonic=filter_mono_exonic,
        skipGTF=skipgtf,
        skipFaFq=skipfafq,
        skipJunction=skipjunction,
    )


if __name__ == "__main__":
    main()
