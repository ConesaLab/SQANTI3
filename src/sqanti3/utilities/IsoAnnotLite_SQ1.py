# Script to generate a GFF3 file from SQANTI3 output and using a tappAS GFF3 as reference.

import logging
import math
import os
import sys
import time
import re
from typing import Optional, Tuple, Dict, List
from tqdm import tqdm

# import argparse
import click
import gtfparse
import pandas as pd

# import bisect

# Global Variables
version = 1.5
CLASS_COLUMN_USED = [0, 1, 2, 3, 5, 6, 7, 30, 32, 33]
CLASS_COLUMN_NAME = [
    "isoform",
    "chrom",
    "strand",
    "length",
    "structural_category",
    "associated_gene",
    "associated_transcript",
    "ORF_length",
    "CDS_start",
    "CDS_end",
]


class transcriptAnnotation:
    def __init__(self, line: str):
        line_items = line.strip("\n").split("\t")
        self.isoform = line_items[0]
        self.chrom = line_items[1]
        self.strand = line_items[2]
        self.length = line_items[3]
        self.exons = line_items[4]
        self.structural_category = line_items[5]
        self.associated_gene = line_items[6]
        self.associated_transcript = line_items[7]
        self.ref_length = line_items[8]
        self.ref_exons = line_items[9]
        self.diff_to_TSS = line_items[10]
        self.diff_to_TTS = line_items[11]
        self.diff_to_gene_TSS = line_items[12]
        self.diff_to_gene_TTS = line_items[13]
        self.subcategory = line_items[14]
        self.RTS_stage = line_items[15]
        self.all_canonical = line_items[16]
        self.min_sample_cov = line_items[17]
        self.min_cov = line_items[18]
        self.min_cov_pos = line_items[19]
        self.sd_cov = line_items[20]
        self.FL = line_items[21]
        self.n_indels = line_items[22]
        self.n_indels_junc = line_items[23]
        self.bite = line_items[24]
        self.iso_exp = line_items[25]
        self.gene_exp = line_items[26]
        self.ratio_exp = line_items[27]
        self.FSM_class = line_items[28]
        self.coding = line_items[29]
        self.ORF_length = line_items[30]
        self.CDS_length = line_items[31]
        self.CDS_start = line_items[32]
        self.CDS_end = line_items[33]
        self.CDS_genomic_start = line_items[34]
        self.CDS_genomic_end = line_items[35]
        self.predicted_NMD = line_items[36]
        self.perc_A_downstream_TTS = line_items[37]
        self.seq_A_downstream_TTS = line_items[38]
        self.dist_to_cage_peak = line_items[39]
        self.within_cage_peak = line_items[40]
        self.dist_to_polya_site = line_items[41]
        self.within_polya_site = line_items[42]
        self.polyA_motif = line_items[43]
        self.polyA_dist = line_items[44]

    def __str__(self):
        _ = (
            f"{self.isoform}\t"
            f"{self.chrom}\t"
            f"{self.strand}\t"
            f"{self.length}\t"
            f"{self.exons}\t"
            f"{self.structural_category}\t"
            f"{self.associated_gene}\t"
            f"{self.associated_transcript}\t"
            f"{self.ref_length}\t"
            f"{self.ref_exons}\t"
            f"{self.diff_to_TSS}\t"
            f"{self.diff_to_TTS}\t"
            f"{self.diff_to_gene_TSS}\t"
            f"{self.diff_to_gene_TTS}\t"
            f"{self.subcategory}\t"
            f"{self.RTS_stage}\t"
            f"{self.all_canonical}\t"
            f"{self.min_sample_cov}\t"
            f"{self.min_cov}\t"
            f"{self.min_cov_pos}\t"
            f"{self.sd_cov}\t"
            f"{self.FL}\t"
            f"{self.n_indels}\t"
            f"{self.n_indels_junc}\t"
            f"{self.bite}\t"
            f"{self.iso_exp}\t"
            f"{self.gene_exp}\t"
            f"{self.ratio_exp}\t"
            f"{self.FSM_class}\t"
            f"{self.coding}\t"
            f"{self.ORF_length}\t"
            f"{self.CDS_length}\t"
            f"{self.CDS_start}\t"
            f"{self.CDS_end}\t"
            f"{self.CDS_genomic_start}\t"
            f"{self.CDS_genomic_end}\t"
            f"{self.predicted_NMD}\t"
            f"{self.perc_A_downstream_TTS}\t"
            f"{self.seq_A_downstream_TTS}\t"
            f"{self.dist_to_cage_peak}\t"
            f"{self.within_cage_peak}\t"
            f"{self.dist_to_polya_site}\t"
            f"{self.within_polya_site}\t"
            f"{self.polyA_motif}\t"
            f"{self.polyA_dist}\t"
        )
        return _

    def __repr__(self):
        _ = (
            f"isoform: {self.isoform}\n"
            f"chrom: {self.chrom}\n"
            f"strand: {self.strand}\n"
            f"length: {self.length}\n"
            f"exons: {self.exons}\n"
            f"structural_category: {self.structural_category}\n"
            f"associated_gene: {self.associated_gene}\n"
            f"associated_transcript: {self.associated_transcript}\n"
            f"ref_length: {self.ref_length}\n"
            f"ref_exons: {self.ref_exons}\n"
            f"diff_to_TSS: {self.diff_to_TSS}\n"
            f"diff_to_TTS: {self.diff_to_TTS}\n"
            f"diff_to_gene_TSS: {self.diff_to_gene_TSS}\n"
            f"diff_to_gene_TTS: {self.diff_to_gene_TTS}\n"
            f"subcategory: {self.subcategory}\n"
            f"RTS_stage: {self.RTS_stage}\n"
            f"all_canonical: {self.all_canonical}\n"
            f"min_sample_cov: {self.min_sample_cov}\n"
            f"min_cov: {self.min_cov}\n"
            f"min_cov_pos: {self.min_cov_pos}\n"
            f"sd_cov: {self.sd_cov}\n"
            f"FL: {self.FL}\n"
            f"n_indels: {self.n_indels}\n"
            f"n_indels_junc: {self.n_indels_junc}\n"
            f"bite: {self.bite}\n"
            f"iso_exp: {self.iso_exp}\n"
            f"gene_exp: {self.gene_exp}\n"
            f"ratio_exp: {self.ratio_exp}\n"
            f"FSM_class: {self.FSM_class}\n"
            f"coding: {self.coding}\n"
            f"ORF_length: {self.ORF_length}\n"
            f"CDS_length: {self.CDS_length}\n"
            f"CDS_start: {self.CDS_start}\n"
            f"CDS_end: {self.CDS_end}\n"
            f"CDS_genomic_start: {self.CDS_genomic_start}\n"
            f"CDS_genomic_end: {self.CDS_genomic_end}\n"
            f"predicted_NMD: {self.predicted_NMD}\n"
            f"perc_A_downstream_TTS: {self.perc_A_downstream_TTS}\n"
            f"seq_A_downstream_TTS: {self.seq_A_downstream_TTS}\n"
            f"dist_to_cage_peak: {self.dist_to_cage_peak}\n"
            f"within_cage_peak: {self.within_cage_peak}\n"
            f"dist_to_polya_site: {self.dist_to_polya_site}\n"
            f"within_polya_site: {self.within_polya_site}\n"
            f"polyA_motif: {self.polyA_motif}\n"
            f"polyA_dist: {self.polyA_dist}\n"
        )
        return _

    def to_dict(self):
        return {
            "isoform": self.isoform,
            "chrom": self.chrom,
            "strand": self.strand,
            "length": self.length,
            "exons": self.exons,
            "structural_category": self.structural_category,
            "associated_gene": self.associated_gene,
            "associated_transcript": self.associated_transcript,
            "ref_length": self.ref_length,
            "ref_exons": self.ref_exons,
            "diff_to_TSS": self.diff_to_TSS,
            "diff_to_TTS": self.diff_to_TTS,
            "diff_to_gene_TSS": self.diff_to_gene_TSS,
            "diff_to_gene_TTS": self.diff_to_gene_TTS,
            "subcategory": self.subcategory,
            "RTS_stage": self.RTS_stage,
            "all_canonical": self.all_canonical,
            "min_sample_cov": self.min_sample_cov,
            "min_cov": self.min_cov,
            "min_cov_pos": self.min_cov_pos,
            "sd_cov": self.sd_cov,
            "FL": self.FL,
            "n_indels": self.n_indels,
            "n_indels_junc": self.n_indels_junc,
            "bite": self.bite,
            "iso_exp": self.iso_exp,
            "gene_exp": self.gene_exp,
            "ratio_exp": self.ratio_exp,
            "FSM_class": self.FSM_class,
            "coding": self.coding,
            "ORF_length": self.ORF_length,
            "CDS_length": self.CDS_length,
            "CDS_start": self.CDS_start,
            "CDS_end": self.CDS_end,
            "CDS_genomic_start": self.CDS_genomic_start,
            "CDS_genomic_end": self.CDS_genomic_end,
            "predicted_NMD": self.predicted_NMD,
            "perc_A_downstream_TTS": self.perc_A_downstream_TTS,
            "seq_A_downstream_TTS": self.seq_A_downstream_TTS,
            "dist_to_cage_peak": self.dist_to_cage_peak,
            "within_cage_peak": self.within_cage_peak,
            "dist_to_polya_site": self.dist_to_polya_site,
            "within_polya_site": self.within_polya_site,
            "polyA_motif": self.polyA_motif,
            "polyA_dist": self.polyA_dist,
        }

    def to_list(self):
        return [_ for _ in self.to_dict().values()]

    def to_Series(self):
        return pd.Series(self.to_dict())


class gtf_fields:
    def __init__(
        self,
        seqname,
        source,
        feature,
        start,
        end,
        score=".",
        strand="",
        frame=".",
        attribute="",
    ):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = attribute

    def __repr__(self):
        return (
            f"seqname: {self.seqname}\n"
            f"source: {self.source}\n"
            f"feature: {self.feature}\n"
            f"start: {self.start}\n"
            f"end: {self.end}\n"
            f"score: {self.score}\n"
            f"strand: {self.strand}\n"
            f"frame: {self.frame}\n"
            f"attributes: {self.attribute}\n"
        )

    def __str__(self):
        return (
            f"{self.seqname}\t"
            f"{self.source}\t"
            f"{self.feature}\t"
            f"{self.start}\t"
            f"{self.end}\t"
            f"{self.score}\t"
            f"{self.strand}\t"
            f"{self.frame}\t"
            f"{self.attribute}\n"
        )

    def to_dict(self):
        _ = {
            "seqname": self.seqname,
            "source": self.source,
            "feature": self.feature,
            "start": self.start,
            "end": self.end,
            "score": self.score,
            "strand": self.strand,
            "frame": self.frame,
            "attribute": self.attribute,
        }
        return _

    def to_Series(self):
        _ = pd.Series(data=self.to_dict())
        return _

    def write(self, filename, mode="a"):
        with open(filename, mode) as output:
            output.write(str(self))


# Functions
def createGTFFromSqanti(file_exons, file_trans, file_junct, filename):
    res = open(filename, "w+")

    desc = ""

    dc_coding = {}
    dc_gene = {}
    dc_SQstrand = {}
    f = open(file_trans)

    # check header
    global CLASS_COLUMN_USED
    global CLASS_COLUMN_NAME

    header = next(f)
    fields = header.split("\t")
    index = 0
    for column in CLASS_COLUMN_NAME:  # check all the columns we used
        if (
            column not in fields[CLASS_COLUMN_USED[index]]
        ):  # if now in the correct possition...
            logging.info(
                f"File classification does not have the correct structure. "
                f" The column '{column}' is not in the possition "
                f"{str(CLASS_COLUMN_USED[index])}"
                f" in the classification file. We have found the column '"
                f"{str(fields[CLASS_COLUMN_USED[index]])}'."
            )
            sys.exit()
        else:
            index = index + 1

    # add transcript, gene and CDS
    for line in f:
        fields = transcriptAnnotation(line)

        ### Create transcript entry
        transcript_entry = gtf_fields(
            seqname=fields.isoform,
            source="tappAS",
            feature="transcript",
            start="1",
            end=fields.length,
            score=".",
            strand=fields.strand,
            frame=".",
            attribute=f"ID={fields.associated_transcript};"
            f"primary_class={fields.structural_category}\n",
        )

        dc_SQstrand.update(
            {str(transcript_entry.transcript): transcript_entry.strand}
        )  # saving strand

        transcript_entry.write(res)

        ### Create gene entry
        gene_entry = gtf_fields(
            seqname=fields.isoform,
            source="tappAS",
            feature="gene",
            start="1",
            end=fields.length,
            score=".",
            strand=fields.strand,
            frame=".",
            attribute=f"ID={fields.associated_gene};"
            f"Name={fields.associated_gene};"
            f"Desc={fields.associated_gene}\n",
        )

        gene_entry.write(res)

        ### Create a CDS entry

        # if the CDS as a 'start' value, create a protein entry too
        if fields.CDS_start != "NA":
            cds_entry = gtf_fields(
                seqname=fields.isoform,
                source="tappAS",
                feature="CDS",
                start=fields.CDS_start,
                end=fields.CDS_end,
                score=".",
                strand=fields.strand,
                frame=".",
                attribute=f"ID=Protein_{fields.isoform};"
                f"Name=Protein_{fields.isoform};"
                f"Desc=Protein_{fields.isoform}\n",
            )

            protein_entry = gtf_fields(
                seqname=cds_entry.seqname,
                source="tappAS",
                feature="protein",
                start="1",
                end=str(
                    int(math.ceil((int(cds_entry.end) - int(cds_entry.start) - 1) / 3))
                ),
                score=".",
                strand=cds_entry.strand,
                frame=".",
                attribute=f"ID=Protein_{cds_entry.seqname};"
                f"Name=Protein_{cds_entry.seqname};"
                f"Desc=Protein_{cds_entry.seqname}\n",
            )
            cds_entry.write(res)
            protein_entry.write(res)
        else:
            cds_entry = gtf_fields(
                seqname=fields.isoform,
                source="tappAS",
                feature="CDS",
                start=".",
                end=".",
                score=".",
                strand=fields.strand,
                frame=".",
                attribute=f"ID=Protein_{fields.isoform};"
                f"Name=Protein_{fields.isoform};"
                f"Desc=Protein_{fields.isoform}\n",
            )
            cds_entry.write(res)

        ### Add an entry to the gene dictionary
        # Gene
        genomic_entry = gtf_fields(
            seqname=fields.associated_gene,
            source="tappAS",
            feature="genomic",
            start="1",
            end="1",
            score=".",
            strand=fields.strand,
            frame=".",
            attribute=f"Chr={fields.chrom};\n",
        )

        category = fields.structural_category
        # Strip the Ensembl version number
        transAssociated = re.sub(
            pattern=r"\.\d+$", repl="", string=fields.associated_transcript
        )

        if genomic_entry.seqname not in dc_gene:
            dc_gene.update(
                {
                    str(genomic_entry.seqname): [
                        fields.associated_gene,
                        fields.structural_category,
                        transAssociated,
                    ]
                }
            )
        else:
            dc_gene.update(
                {
                    str(genomic_entry.seqname): dc_gene[genomic_entry.seqname]
                    + [
                        fields.associated_gene,
                        fields.structural_category,
                        transAssociated,
                    ]
                }
            )

        ### Coding Dictionary
        CDSstart = fields.CDS_start
        CDSend = fields.CDS_end
        orf = fields.FSM_class

        if not dc_coding.get(cds_entry.seqname):
            dc_coding.update({str(cds_entry.seqname): [CDSstart, CDSend, orf]})
        else:
            dc_coding.update(
                {
                    str(cds_entry.seqname): dc_coding.get(cds_entry.seqname)
                    + [CDSstart, CDSend, orf]
                }
            )

        genomic_entry.write(res)

        # Write TranscriptAttributes
        if CDSstart != "NA":
            # 3'UTR
            tp_utr = gtf_fields(
                seqname=fields.associated_gene,
                source="TranscriptAttributes",
                feature="3UTR_Length",
                start=int(CDSend) + 1,
                end=fields.length,
                score=".",
                strand=fields.strand,
                frame=".",
                attribute="ID=3UTR_Length; Name=3UTR_Length; Desc=3UTR_Length\n",
            )

            tp_utr.write(res)

            # 5'UTR
            fp_utr = gtf_fields(
                seqname=fields.associated_gene,
                source="TranscriptAttributes",
                feature="5UTR_Length",
                start=1,
                end=int(fields.ORF_length) - 1 + 1,
                score=".",
                strand=fields.strand,
                frame=".",
                attribute="ID=5UTR_Length;Name=5UTR_Length;Desc=5UTR_Length\n",
            )

            fp_utr.write(res)

            # CDS
            cds = gtf_fields(
                seqname=fields.associated_gene,
                source="TranscriptAttributes",
                feature="CDS",
                start=CDSstart,
                end=CDSend,
                score=".",
                strand=fields.strand,
                frame=".",
                attribute="ID=CDS;Name=CDS;Desc=CDS\n",
            )

            res.write(str(cds))

            # polyA
            polyA = gtf_fields(
                seqname=fields.associated_gene,
                source="TranscriptAttributes",
                feature="polyA_Site",
                start=fields.length,
                end=fields.length,
                score=".",
                strand=fields.strand,
                frame=".",
                attribute="ID=polyA_Site;Name=polyA_Site;Desc=polyA_Site\n",
            )

            res.write(str(polyA))

    f.close()

    exons = gtfparse.parse_gtf_and_expand_attributes(gtf)
    exons = exons[exons[:, "feature"] == "exon"]
    dc_exons = {}
    # add exons
    for exon in exons:
        fields = gtf_fields(
            seqname=exon.seqname,
            source=exon.source,
            feature=exon.feature,
            start=exon.start,
            end=exon.end,
            score=".",
            strand=exon.strand,
            frame=".",
            attribute=f"Chr={exon.seqname}",
        )

        # Exons Dictionary
        if not dc_exons.get(exon.gene_id):
            dc_exons.update({exon.gene_id: [[exon.start, exon.end]]})
        else:
            dc_exons.update(
                {exon.gene_id: dc_exons.get(exon.gene_id) + [[exon.start, exon.end]]}
            )

        res.write(str(fields))

    # add junctions
    junctions = pd.read_csv(file_junct, delimiter="\t")
    # header

    for line in junctions.iterrows():

        # Junctions file can have a dvierse number of columns, not only 19 but 0-14 are allways the same
        junct = gtf_fields(
            seqname=line[1].isoform,
            source="tappAS",
            feature="splice_junction",
            start=line[1].genomic_start_coord,
            end=line[1].genomic_end_coord,
            score=".",
            strand=line[1].strand,
            frame=".",
            attribute=f"ID={line[1].junction_number}_{line[1].canonical};Chr={line[1].chrom}\n",
        )

        res.write(str(junct))

    res.close()

    return dc_exons, dc_coding, dc_gene, dc_SQstrand


def readGFF(gff3):
    f = open(gff3)
    # create dictionary for each transcript and dictionary for exons
    dc_GFF3 = {}
    dc_GFF3exonsTrans = {}
    dc_GFF3transExons = {}
    dc_GFF3coding = {}
    dc_GFF3strand = {}
    for line in f:
        fields = line.split("\t")
        if len(fields) == 9:
            # feature (transcript, gene, exons...)
            transcript = fields[0]
            feature = fields[2]
            start = fields[3]
            end = fields[4]
            strand = fields[6]

            if not strand == ".":
                dc_GFF3strand.update({str(transcript): strand})  # saving strand

            text = fields[8].split(" ")
            if not text[-1].endswith("\n"):
                line = line + "\n"

            if feature == "exon":
                if not dc_GFF3transExons.get(str(transcript)):
                    dc_GFF3transExons.update(
                        {str(transcript): [[int(start), int(end)]]}
                    )
                else:
                    dc_GFF3transExons.update(
                        {
                            str(transcript): dc_GFF3transExons.get(str(transcript))
                            + [[int(start), int(end)]]
                        }
                    )

                if not dc_GFF3exonsTrans.get(int(start)):
                    dc_GFF3exonsTrans.update({int(start): [transcript]})
                else:
                    dc_GFF3exonsTrans.update(
                        {int(start): dc_GFF3exonsTrans.get(int(start)) + [transcript]}
                    )
            elif feature == "CDS":
                if not dc_GFF3coding.get(str(transcript)):
                    dc_GFF3coding.update(
                        {str(transcript): [int(start), int(end)]}
                    )  # not int bc can be NA
                else:
                    dc_GFF3coding.update(
                        {
                            str(transcript): dc_GFF3coding.get(str(transcript))
                            + [int(start), int(end)]
                        }
                    )

            elif feature in [
                "splice_junction",
                "transcript",
                "gene",
                "protein",
                "genomic",
            ]:
                continue

            else:
                if not dc_GFF3.get(transcript):
                    dc_GFF3.update({str(transcript): [[start, end, line]]})
                else:
                    dc_GFF3.update(
                        {
                            str(transcript): dc_GFF3.get(transcript)
                            + [[start, end, line]]
                        }
                    )
        else:
            print("File GFF3 doesn't have the correct number of columns (9).")

    sorted(dc_GFF3exonsTrans.keys())
    return dc_GFF3, dc_GFF3exonsTrans, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand


def unique(list1):
    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)


def transformTransFeaturesToGenomic(
    dc_GFF3, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand
):
    newdc_GFF3 = {}
    bnegative = False

    for trans in dc_GFF3transExons.keys():
        if not dc_GFF3.get(trans):
            continue
        annot = dc_GFF3.get(trans)
        for values in annot:
            bProt = False
            line = values[2]
            fields = line.split("\t")
            text = fields[8].split(" ")
            strand = dc_GFF3strand.get(trans)

            start = 0
            end = 0
            startG = 0
            endG = 0
            # Transcript calculate normal - include CDS
            if text[-1].endswith("T\n") and not fields[3] == ".":
                start = int(fields[3])
                end = int(fields[4])
            elif text[-1].endswith("P\n") and not fields[3] == ".":
                start = int(fields[3])
                end = int(fields[4])
                bProt = True
            else:
                if not newdc_GFF3.get(trans):
                    newdc_GFF3.update({str(trans): [values]})
                    continue
                else:
                    newdc_GFF3.update({str(trans): newdc_GFF3.get(trans) + [values]})
                    continue

            totalDiff = end - start
            if not bProt:
                allExons = dc_GFF3transExons.get(trans)
            else:
                allExons = dc_GFF3coding.get(trans)
                if not allExons:
                    continue

            if strand == "+":
                allExons = sorted(allExons)
            else:
                allExons = sorted(allExons, reverse=True)

            bstart = False
            bend = False
            for exon in allExons:
                if totalDiff < 0:
                    bnegative = True
                    break

                # START already found
                if bstart:
                    if exon[0] + totalDiff - 1 <= exon[1]:  # pos ends here
                        end = exon[0] + totalDiff
                        endG = end
                        bend = True
                    else:  # pos ends in other exon and we add the final exon
                        totalDiff = totalDiff - (exon[1] - exon[0] + 1)

                # Search for START
                if exon[1] - exon[0] + 1 >= start and not bstart:  # pos starts here
                    start = exon[0] + int(start) - 1
                    startG = start
                    bstart = True
                    if start + totalDiff - 1 <= exon[1]:  # pos ends here
                        end = start + totalDiff
                        endG = end
                        bend = True
                    else:  # pos ends in other exon and we add the final exon
                        totalDiff = totalDiff - (exon[1] - start + 1)
                else:
                    # not in first exon, update the start and end pos substrating exon length
                    start = start - (exon[1] - exon[0] + 1)
                    end = end - (exon[1] - exon[0] + 1)

                if bend:
                    if not bProt:
                        if strand == "+":
                            newline = (
                                fields[0]
                                + "\t"
                                + fields[1]
                                + "\t"
                                + fields[2]
                                + "\t"
                                + str(startG)
                                + "\t"
                                + str(endG)
                                + "\t"
                                + fields[5]
                                + "\t"
                                + fields[6]
                                + "\t"
                                + fields[7]
                                + "\t"
                                + fields[8]
                            )
                        else:
                            newline = (
                                fields[0]
                                + "\t"
                                + fields[1]
                                + "\t"
                                + fields[2]
                                + "\t"
                                + str(endG)
                                + "\t"
                                + str(startG)
                                + "\t"
                                + fields[5]
                                + "\t"
                                + fields[6]
                                + "\t"
                                + fields[7]
                                + "\t"
                                + fields[8]
                            )
                        if not newdc_GFF3.get(trans):
                            newdc_GFF3.update({str(trans): [[startG, endG, newline]]})
                            break
                        else:
                            newdc_GFF3.update(
                                {
                                    str(trans): newdc_GFF3.get(trans)
                                    + [[startG, endG, newline]]
                                }
                            )
                            break
                    else:
                        if strand == "-":
                            aux = startG
                            startG = endG
                            endG = aux
                        if not newdc_GFF3.get(trans):
                            newdc_GFF3.update({str(trans): [[startG, endG, values[2]]]})
                            break
                        else:
                            newdc_GFF3.update(
                                {
                                    str(trans): newdc_GFF3.get(trans)
                                    + [[startG, endG, values[2]]]
                                }
                            )
                            break
            if bnegative:
                break

    return newdc_GFF3


def transformTransFeaturesToLocale(dc_GFF3, dc_SQexons):
    dc_newGFF3 = {}
    for trans in dc_GFF3.keys():
        annot = dc_GFF3.get(trans)
        line = annot[0]
        line = line.split("\t")
        strand = line[6]

        exons = dc_SQexons.get(trans)
        if strand == "+":
            exons = sorted(exons)
        else:
            exons = sorted(exons, reverse=True)
        start = 0
        end = 0
        for line in annot:
            fields = line.split("\t")
            text = fields[8].split(" ")
            if fields[1] == "tappAS":
                continue

            elif text[-1].endswith("T\n"):

                if strand == "+":
                    startG = fields[3]
                    endG = fields[4]
                else:
                    startG = fields[4]
                    endG = fields[3]
                bstart = False
                bend = False
                distance = 0  # other exons
                for ex in exons:
                    if not startG == "." or not endG == ".":
                        # SEARCH FOR START
                        if (
                            int(startG) >= int(ex[0])
                            and int(startG) <= int(ex[1])
                            and not bstart
                        ):  # start
                            start = (int(startG) - int(ex[0]) + 1) + distance
                            bstart = True
                            if int(endG) >= int(ex[0]) and int(endG) <= int(ex[1]):
                                end = start + (int(endG) - int(startG) + 1) - 1
                                bend = True
                                break
                            else:
                                distance = int(ex[1]) - int(startG) + 1
                                continue

                        elif not bstart:
                            distance = distance + (int(ex[1]) - int(ex[0]) + 1)

                        # SEARCH FOR END
                        if bstart:
                            if int(endG) >= int(ex[0]) and int(endG) <= int(ex[1]):
                                end = (
                                    (int(endG) - int(ex[0]) + 1) + distance + start - 1
                                )
                                bend = True
                                break
                            else:
                                distance = distance + (int(ex[1]) - int(ex[0]) + 1)
                    else:
                        start = startG
                        end = endG
                        bend = True
                        break
                if bend:  # to be sure in full-spliced match cases
                    newline = (
                        fields[0]
                        + "\t"
                        + fields[1]
                        + "\t"
                        + fields[2]
                        + "\t"
                        + str(start)
                        + "\t"
                        + str(end)
                        + "\t"
                        + fields[5]
                        + "\t"
                        + fields[6]
                        + "\t"
                        + fields[7]
                        + "\t"
                        + fields[8]
                    )
                    if not dc_newGFF3.get(trans):
                        dc_newGFF3.update({str(trans): [newline]})
                    else:
                        dc_newGFF3.update(
                            {str(trans): dc_newGFF3.get(trans) + [newline]}
                        )
            else:
                if not dc_newGFF3.get(trans):
                    dc_newGFF3.update({str(trans): [line]})
                else:
                    dc_newGFF3.update({str(trans): dc_newGFF3.get(trans) + [line]})
    return dc_newGFF3


def transformProtFeaturesToLocale(dc_GFF3, dc_SQexons, dc_SQcoding):
    dc_newGFF3 = {}
    for trans in dc_GFF3.keys():
        annot = dc_GFF3.get(trans)
        line = annot[0]
        line = line.split("\t")
        strand = line[6]

        exons = dc_SQexons.get(trans)
        if strand == "+":
            exons = sorted(exons)
        else:
            exons = sorted(exons, reverse=True)

        annot = dc_GFF3.get(trans)

        start = 0
        end = 0
        if not dc_SQcoding.get(trans):
            continue
        startcoding = dc_SQcoding.get(trans)[0]
        startcoding = startcoding[0]
        if startcoding == "NA":
            continue

        for line in annot:
            fields = line.split("\t")
            text = fields[8].split(" ")
            if fields[1] == "tappAS":
                continue
            if text[-1].endswith("P\n"):
                startG = fields[3]
                endG = fields[4]

                bstart = False
                CDSstart = False
                distance = 0  # other exons
                for ex in exons:
                    if not startG == "." or not endG == ".":
                        if not CDSstart:  # CDS start
                            if (
                                int(startcoding) >= int(ex[0])
                                and int(startcoding) <= int(ex[1])
                                and not bstart
                            ):  # start

                                start = (
                                    int(startG) - int(ex[0]) + 1 + distance
                                )  # CDSstart
                                CDSstart = True
                            else:
                                distance = distance + int(ex[1]) - int(startG) + 1
                        if (
                            int(startG) >= int(ex[0])
                            and int(startG) <= int(ex[1])
                            and not bstart
                            and CDSstart
                        ):  # start

                            start = (
                                int(startG) - int(start) + 1 + distance
                            )  # diff between genomic pos and CDSstart
                            bstart = True
                            if int(endG) >= int(ex[0]) and int(endG) <= int(ex[1]):
                                end = start + int(endG) - int(startG) + 1
                                break
                            else:
                                distance = distance + int(ex[1]) - int(startG) + 1
                        else:
                            distance = int(ex[1]) - int(ex[0]) + 1
                        if bstart and CDSstart:
                            if int(endG) >= int(ex[0]) and int(endG) <= int(ex[1]):
                                end = int(endG) - int(ex[0]) + 1 + distance
                                break
                            else:
                                distance = distance + (int(ex[1]) - int(ex[0]) + 1)
                    else:
                        start = startG
                        end = endG
                newline = (
                    fields[0]
                    + "\t"
                    + fields[1]
                    + "\t"
                    + fields[2]
                    + "\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t"
                    + fields[5]
                    + "\t"
                    + fields[6]
                    + "\t"
                    + fields[7]
                    + "\t"
                    + fields[8]
                )
                if not dc_newGFF3.get(trans):
                    dc_newGFF3.update({str(trans): [newline]})
                else:
                    dc_newGFF3.update({str(trans): dc_newGFF3.get(trans) + [newline]})
            else:
                if not dc_newGFF3.get(trans):
                    dc_newGFF3.update({str(trans): [line]})
                else:
                    dc_newGFF3.update({str(trans): dc_newGFF3.get(trans) + [line]})
    return dc_newGFF3


def transformCDStoGenomic(dc_SQcoding, dc_SQexons, dc_SQstrand):
    newdc_coding = {}
    bnegative = False

    for trans in dc_SQcoding.keys():
        newCDS = []
        aux = []
        CDS = dc_SQcoding.get(trans)

        if CDS[0] == "NA":
            if not newdc_coding.get(str(trans)):
                newdc_coding.update({str(trans): [CDS]})
            else:
                newdc_coding.update({str(trans): newdc_coding.get(str(trans)) + [CDS]})
            continue

        totalDiff = int(CDS[1]) - int(CDS[0])

        allExons = dc_SQexons.get(trans)
        if not allExons:
            continue

        if dc_SQstrand.get(trans) == "+":
            allExons = sorted(allExons)
        else:
            allExons = sorted(allExons, reverse=True)
        bstart = False
        bend = False
        start = 0
        end = 0
        for exon in allExons:
            if totalDiff < 0:
                print("The difference can't be negative.")
                bnegative = True
                break

            # START already found
            if bstart:
                if exon[0] + totalDiff - 1 <= exon[1]:  # CDS ends here
                    end = exon[0] + totalDiff - 1
                    aux = [[exon[0], end]]
                    newCDS = newCDS + aux
                    bend = True
                else:  # CDS ends in other exon and we add the final exon
                    aux = [[exon[0], exon[1]]]
                    newCDS = newCDS + aux
                    totalDiff = totalDiff - (exon[1] - exon[0] + 1)

            # Search for START
            if exon[1] - exon[0] + 1 >= int(CDS[0]) and not bstart:  # CDS starts here
                start = exon[0] + int(CDS[0]) - 1
                bstart = True
                if start + totalDiff - 1 <= exon[1]:  # CDS ends here
                    end = start + totalDiff - 1
                    aux = [[start, end]]
                    newCDS = newCDS + aux
                    bend = True
                else:  # CDS ends in other exon and we add the final exon
                    aux = [[start, exon[1]]]
                    newCDS = newCDS + aux
                    totalDiff = totalDiff - (exon[1] - start + 1)

            if bend:
                if not newdc_coding.get(str(trans)):
                    newdc_coding.update({str(trans): newCDS})
                else:
                    newdc_coding.update(
                        {str(trans): newdc_coding.get(str(trans)) + newCDS}
                    )
                break
        if bnegative:
            break

    return newdc_coding


def checkSameCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, strand):
    coding = True
    semicoding = True
    total_semi = 0
    total_annot = 0
    if dc_SQcoding.get(transSQ) and dc_GFF3coding.get(transGFF3):
        if not dc_SQcoding.get(transSQ)[0][0] == "NA":
            # Tenemos rango de intervalos en los exones:
            #   Si coinciden todos es coding
            #   Si coinciden todos menos sub exons (inicio o final) es semicoding
            allExonsGFF3 = dc_GFF3coding.get(transGFF3)
            if strand == "+":
                allExonsGFF3 = sorted(allExonsGFF3)
            else:
                allExonsGFF3 = sorted(allExonsGFF3, reverse=True)
            for ex in allExonsGFF3:
                allExonsSQ = dc_SQcoding.get(transSQ)
                if strand == "+":
                    allExonsSQ = sorted(allExonsSQ)
                else:
                    allExonsSQ = sorted(allExonsSQ, reverse=True)

                if ex in allExonsSQ:
                    total_annot = total_annot + 1
                    continue
                else:
                    coding = False
                    semicoding = False  # Check if we found semicoding
                    for exSQ in allExonsSQ:
                        if ex[0] <= exSQ[0] and exSQ[1] <= ex[1]:  # Region inside
                            total_semi = total_semi + 1
                            semicoding = True
                            break
                        elif (
                            exSQ[0] <= ex[0] and exSQ[1] <= ex[1]
                        ):  # or region bigger by left
                            total_semi = total_semi + 1
                            semicoding = True
                            break
                        elif (
                            ex[0] <= exSQ[0] and ex[1] <= exSQ[1]
                        ):  # or region bigger by right
                            total_semi = total_semi + 1
                            semicoding = True
                            break
                        elif (
                            exSQ[0] <= ex[0] and ex[1] <= exSQ[1]
                        ):  # or region bigger by both sides
                            total_semi = total_semi + 1
                            semicoding = True
                            break

        if total_annot == len(dc_GFF3coding.get(transGFF3)) and not dc_GFF3coding.get(
            transGFF3
        ):
            coding = True
        elif total_annot > 0 or total_semi > 0:
            semicoding = True
        else:
            coding = False
            semicoding = False
    return coding, semicoding


def checkFeatureInCDS(
    dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, start, end, strand
):
    bstart = False
    if (
        not dc_SQcoding.get(transSQ)[0] == "NA"
        and dc_SQcoding.get(transSQ)
        and dc_GFF3coding.get(transGFF3)
    ):
        # Tenemos rango de intervalos en los exones:
        #   Si coinciden todos es coding
        #   Si coinciden todos menos sub exons (inicio o final) es semicoding
        allExonsGFF3 = dc_GFF3coding.get(transGFF3)
        allExonsSQ = dc_SQcoding.get(transSQ)

        if strand == "+":
            allExonsGFF3 = sorted(allExonsGFF3)
            allExonsGFF3 = sorted(allExonsSQ)
        else:
            allExonsGFF3 = sorted(allExonsGFF3, reverse=True)
            allExonsSQ = sorted(allExonsSQ, reverse=True)

        for ex in allExonsGFF3:

            #########
            #  END  #
            #########
            if bstart:
                if ex[0] <= end and end <= ex[1]:  # end in exon
                    if ex in allExonsSQ:  # if exon exist
                        return True
                    else:
                        for exSQ in allExonsSQ:  # we just need end subexon
                            if (
                                exSQ[0] == ex[0] and end <= exSQ[1]
                            ):  # and feature in range
                                return True
                        return False  # doesnt find the feture in same exon
                else:  # in next exon
                    if not ex in allExonsSQ:
                        return (
                            False
                        )  # end in another exons and we don't have that intermediate in SQ
                    else:
                        continue

            #########
            # START #
            #########
            if ex[0] <= start and start <= ex[1] and not bstart:  # start in exon
                if ex[0] <= end and end <= ex[1]:  # end in exon
                    if ex in allExonsSQ:  # if exon exist
                        return True
                    else:
                        for exSQ in allExonsSQ:  # we just need start and end in subexon
                            if (
                                exSQ[0] <= start
                                and start <= exSQ[1]
                                and exSQ[0] <= end
                                and end <= exSQ[1]
                            ):  # and feature in range
                                return True
                        return False  # doesnt find the feture in same exon
                else:  # we need an exSQ that ends in same position to continue
                    for exSQ in allExonsSQ:
                        if (
                            exSQ[0] <= start and ex[1] == exSQ[1]
                        ):  # at begining just start but end the same
                            bstart = True
    return False


def checkFeatureInTranscript(
    dc_SQexons, dc_GFF3transExons, transSQ, transGFF3, start, end, strand
):
    bstart = False
    bnotMiddleExon = False
    if dc_SQexons.get(transSQ) and dc_GFF3transExons.get(transGFF3):
        # Tenemos rango de intervalos en los exones:
        #   Si coinciden todos es coding
        #   Si coinciden todos menos sub exons (inicio o final) es semicoding
        allExonsGFF3 = dc_GFF3transExons.get(transGFF3)
        allExonsSQ = dc_SQexons.get(transSQ)
        if strand == "+":
            allExonsGFF3 = sorted(allExonsGFF3)
            allExonsSQ = sorted(allExonsSQ)
        else:
            allExonsGFF3 = sorted(allExonsGFF3, reverse=True)
            allExonsSQ = sorted(allExonsSQ, reverse=True)

        for ex in allExonsGFF3:
            if ex[0] <= start and start <= ex[1] and not bstart:  # Annot in exon
                for exSQ in allExonsSQ:  # Look for Start
                    if ex[0] <= end and end <= ex[1]:  # also end it's here
                        if ex in allExonsSQ:  # if exon exist
                            return True
                        elif (
                            exSQ[0] <= start
                            and start <= exSQ[1]
                            and exSQ[0] <= end
                            and end <= exSQ[1]
                        ):  # case when we have the end at same exon but with different length (although same genomic positions
                            return True

                    elif (
                        exSQ[0] <= start and start <= exSQ[1] and ex[1] == exSQ[1]
                    ):  # end in another exon, we need same ending
                        bstart = True
            elif bstart and ex[0] <= end and end <= ex[1]:  # End Annot in exon
                for exSQ in allExonsSQ:  # Look for End
                    if ex in allExonsSQ:  # if exon exist
                        return True
                    elif (
                        exSQ[0] <= end and end <= exSQ[1] and ex[0] == exSQ[0]
                    ):  # end in another exon, we need same exon start
                        return True
                    else:  # we need same exon
                        if not ex[0] == exSQ[0] and not ex[1] == exSQ[1]:
                            bnotMiddleExon = True  # We don't found the middle Exons
                            break
            else:
                continue

            if bnotMiddleExon:
                break

    return False


def mappingFeatures(
    dc_SQexons,
    dc_SQcoding,
    dc_SQtransGene,
    dc_GFF3exonsTrans,
    dc_GFF3transExons,
    dc_GFF3,
    dc_GFF3coding,
    filename,
):
    f = open(filename, "a+")
    print("\n")
    transcriptsAnnotated = 0
    totalAnotations = 0
    featuresAnnotated = 0
    for transSQ in dc_SQexons.keys():

        # Be carefully - not all tranSQ must be in SQtransGene
        if not dc_SQtransGene.get(str(transSQ)):
            continue

        perct = transcriptsAnnotated / len(dc_SQexons) * 100
        print("\t" + "%.2f" % perct + " % of transcripts annotated...", end="\r")

        #######################
        # IF FULL-SPLICED-MATCH#
        #######################
        infoGenomic = dc_SQtransGene.get(transSQ)
        transGFF3 = infoGenomic[2]

        ###########################
        # IF NOT FULL-SPLICED-MATCH#
        ###########################
        val = ""
        if dc_GFF3.get(transGFF3):  # Novel Transcript won't be annoted
            val = dc_GFF3.get(transGFF3)
        elif dc_GFF3.get(transSQ):
            transGFF3 = transSQ
            val = dc_GFF3.get(transGFF3)
        else:
            continue

        line = val[0][2].split("\t")
        strand = line[6]
        # Check if we had same CDS to add Protein information
        coding, semicoding = checkSameCDS(
            dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, strand
        )

        for values in dc_GFF3.get(transGFF3):
            fields = values[2].split("\t")
            text = fields[8].split(" ")
            strand = fields[6]
            if fields[1] == "tappAS":
                continue
            totalAnotations = totalAnotations + 1
            ####################
            # PROTEIN ANNOTATION#
            ####################
            if (
                text[-1].endswith("P\n")
                or text[-1].endswith("G\n")
                or text[-1].endswith("N\n")
            ):  # protein
                if coding:
                    index = values[2].find("\t")
                    if values[2].endswith("\n"):
                        featuresAnnotated = featuresAnnotated + 1
                        f.write(transSQ + values[2][index:])  # write line
                    else:
                        featuresAnnotated = featuresAnnotated + 1
                        f.write(transSQ + values[2][index:] + "\n")  # write line

                elif semicoding and not values[0] == "." and not values[1] == ".":
                    bannot = False
                    # funcion match annot to its our CDSexons and match to CDSexonsSQ
                    bannot = checkFeatureInCDS(
                        dc_SQcoding,
                        dc_GFF3coding,
                        transSQ,
                        transGFF3,
                        int(values[0]),
                        int(values[1]),
                        strand,
                    )
                    if bannot:
                        index = values[2].find("\t")
                        if values[2].endswith("\n"):
                            featuresAnnotated = featuresAnnotated + 1
                            f.write(transSQ + values[2][index:])  # write line
                        else:
                            featuresAnnotated = featuresAnnotated + 1
                            f.write(transSQ + values[2][index:] + "\n")  # write line

                elif semicoding and values[0] == "." and values[1] == ".":
                    index = values[2].find("\t")
                    if values[2].endswith("\n"):
                        featuresAnnotated = featuresAnnotated + 1
                        f.write(transSQ + values[2][index:])  # write line
                    else:
                        featuresAnnotated = featuresAnnotated + 1
                        f.write(transSQ + values[2][index:] + "\n")  # write line

            #######################
            # TRANSCRIPT ANNOTATION#
            #######################

            if (
                not values[0] == "."
                and not values[1] == "."
                and text[-1].endswith("T\n")
            ):
                bannot = False
                bannot = checkFeatureInTranscript(
                    dc_SQexons,
                    dc_GFF3transExons,
                    transSQ,
                    transGFF3,
                    int(values[0]),
                    int(values[1]),
                    strand,
                )

                if bannot:
                    index = values[2].find("\t")
                    if values[2].endswith("\n"):
                        featuresAnnotated = featuresAnnotated + 1
                        f.write(transSQ + values[2][index:])  # write line
                    else:
                        featuresAnnotated = featuresAnnotated + 1
                        f.write(transSQ + values[2][index:] + "\n")  # write line
        transcriptsAnnotated = transcriptsAnnotated + 1
    f.close()

    print(
        "\n\n\t·Annoted a total of "
        + str(featuresAnnotated)
        + " annotation features from reference GFF3 file."
    )
    perct = featuresAnnotated / totalAnotations * 100
    print(
        "\t·Annoted a total of "
        + "%.2f" % perct
        + " % of the reference GFF3 file annotations.\n\n"
    )


# UPDATE GFF3 - new columns information
def addPosType(res, line, posType):
    if line.endswith(";"):
        res.write(line + " PosType=" + posType + "\n")
    else:
        res.write(line[:-1] + "; PosType=" + posType + "\n")


def updateGTF(filename, filenameMod):
    # open new file
    res = open(filenameMod, "w")
    # open annotation file and process all data
    with open(filename, "r") as f:
        # process all entries - no header line in file
        for line in f:
            if len(line) == 0:
                break
            else:
                if line and line[0] != "#":
                    fields = line.split("\t")
                    if len(fields) == 9:

                        text = fields[8].split(" ")
                        if text[-1].startswith("PosType"):
                            res.write(line)

                        elif fields[1] == "tappAS":
                            if fields[2] == "transcript":
                                addPosType(res, line, "T")
                            elif fields[2] == "gene":
                                addPosType(res, line, "T")
                            elif fields[2] == "CDS":
                                addPosType(res, line, "T")
                            elif fields[2] == "genomic":
                                addPosType(res, line, "G")
                            elif fields[2] == "exon":
                                addPosType(res, line, "G")
                            elif fields[2] == "splice_junction":
                                addPosType(res, line, "G")
                            elif fields[2] == "protein":
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "COILS":
                            if fields[2] == "COILED":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "GeneOntology":
                            if fields[2] in ("C", "cellular_component"):
                                addPosType(res, line, "N")
                            elif fields[2] in ("F", "molecular_function"):
                                addPosType(res, line, "N")
                            elif fields[2] in ("P", "biological_process"):
                                addPosType(res, line, "N")
                            elif fields[2] in ("eco"):
                                addPosType(res, line, "N")  # Fran tomato annot
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                ##break

                        elif fields[1] == "MOBIDB_LITE":
                            if fields[2] == "DISORDER":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "NMD":
                            if fields[2] == "NMD":
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] in ("PAR-CLIP", "PAR-clip"):
                            if fields[2] in (
                                "RNA_binding",
                                "RNA_Binding_Protein",
                                "RBP_Binding",
                            ) or fields[2].startswith("RNA_binding_"):
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] == "PFAM":
                            if fields[2] == "DOMAIN":
                                addPosType(res, line, "P")
                            elif fields[2] in ("CLAN", "clan"):
                                addPosType(res, line, "N")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                # break

                        elif fields[1] == "Provean":
                            if fields[2] == "FunctionalImpact":
                                addPosType(res, line, "N")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                # break

                        elif fields[1] in ("REACTOME", "Reactome"):
                            if fields[2] in ("PATHWAY", "pathway", "Pathway"):
                                addPosType(res, line, "N")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                # break

                        elif fields[1] == "RepeatMasker":
                            if fields[2] == "repeat":
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] == "SIGNALP_EUK":
                            if fields[2] == "SIGNAL":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "TMHMM":
                            if fields[2] == "TRANSMEM":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "TranscriptAttributes":
                            addPosType(res, line, "T")

                        elif fields[1] == "UTRsite":
                            if fields[2] == "uORF":
                                addPosType(res, line, "T")
                            elif fields[2] == "5UTRmotif":
                                addPosType(res, line, "T")
                            elif fields[2] == "PAS":
                                addPosType(res, line, "T")
                            elif fields[2] == "3UTRmotif":
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] in (
                            "UniProtKB/Swiss-Prot_Phosphosite",
                            "Swissprot_Phosphosite",
                        ):
                            if fields[2] == "ACT_SITE":
                                addPosType(res, line, "P")
                            elif fields[2] == "BINDING":
                                addPosType(res, line, "P")
                            elif fields[2] == "PTM":
                                addPosType(res, line, "P")
                            elif fields[2] == "MOTIF":
                                addPosType(res, line, "P")
                            elif fields[2] == "COILED":
                                addPosType(res, line, "P")
                            elif fields[2] == "TRANSMEM":
                                addPosType(res, line, "P")
                            elif fields[2] == "COMPBIAS":
                                addPosType(res, line, "P")
                            elif fields[2] == "INTRAMEM":
                                addPosType(res, line, "P")
                            elif fields[2] == "NON_STD":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] in ("cNLS_mapper", "NLS_mapper"):
                            if fields[2] == "MOTIF":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] in ("miRWalk", "mirWalk"):
                            if fields[2] in ("miRNA", "miRNA_Binding"):
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] == "scanForMotifs":
                            if fields[2] == "PAS":
                                addPosType(res, line, "T")
                            elif fields[2] in ("3UTRmotif", "3'UTRmotif"):
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] == "MetaCyc":
                            if fields[2] == "pathway":
                                addPosType(res, line, "N")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                # break

                        elif fields[1] == "KEGG":
                            if fields[2] in ("pathway", "Pathway"):
                                addPosType(res, line, "N")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                # break

                        elif fields[1] == "SUPERFAMILY":
                            if fields[2] == "DOMAIN":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "SMART":
                            if fields[2] == "DOMAIN":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "TIGRFAM":
                            if fields[2] == "DOMAIN":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "psRNATarget":
                            if fields[2] == "miRNA":
                                addPosType(res, line, "T")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using T type to annotate."
                                )
                                addPosType(res, line, "T")
                                # break

                        elif fields[1] == "CORUM":
                            if fields[2] == "Complex":
                                addPosType(res, line, "P")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using P type to annotate."
                                )
                                addPosType(res, line, "P")
                                # break

                        elif fields[1] == "Orthologues":
                            if fields[2] == "S.tuberosum":
                                addPosType(res, line, "N")
                            elif fields[2] in ("A.thaliana"):
                                addPosType(res, line, "N")
                            else:
                                print(
                                    "IsoAnnotLite can not identify the feature "
                                    + str(fields[2])
                                    + " in source "
                                    + str(fields[1])
                                    + ", using N type to annotate."
                                )
                                addPosType(res, line, "N")
                                # break

                        else:
                            print(
                                "IsoAnnotLite can not identify the source "
                                + str(fields[1])
                                + ", in line:\n"
                                + line
                                + "\nUSing N type to annotate."
                            )
                            addPosType(res, line, "N")
                            # break

                    else:
                        print("Error in line (has not 9 fields):\n" + line)
                        break

        res.close()


def readGFFandGetData(filenameMod):
    # open annotation file and process all data
    dcTrans = {}
    dcExon = {}
    dcTransFeatures = {}
    dcGenomic = {}
    dcSpliceJunctions = {}
    dcProt = {}
    dcProtFeatures = {}
    dcTranscriptAttributes = {}

    dcTransID = {}

    with open(filenameMod, "r") as f:
        # process all entries - no header line in file
        for line in f:
            if len(line) == 0:
                break
            else:
                if line and line[0] != "#":
                    fields = line.split("\t")
                    if len(fields) == 9:

                        transcript = fields[0]
                        text = fields[8].split(" ")
                        # transcriptID = text[0]
                        # transcriptID = transcriptID[3:-1]

                        if fields[1] == "tappAS":
                            if fields[2] in ["transcript", "gene", "CDS"]:
                                if not dcTrans.get(str(transcript)):
                                    dcTrans.update({str(transcript): [line]})
                                else:
                                    dcTrans.update(
                                        {
                                            str(transcript): dcTrans.get(
                                                str(transcript)
                                            )
                                            + [line]
                                        }
                                    )
                                # extra dcTransID
                                # if not dcTransID.get(str(transcriptID)):
                                #    dcTransID.update({str(transcriptID) : [line]})
                                # else:
                                #    dcTransID.update({str(transcriptID) : dcTransID.get(str(transcriptID)) + [line]})
                            elif fields[2] in ["exon"]:
                                if not dcExon.get(str(transcript)):
                                    dcExon.update({str(transcript): [line]})
                                else:
                                    dcExon.update(
                                        {
                                            str(transcript): dcExon.get(str(transcript))
                                            + [line]
                                        }
                                    )
                            elif fields[2] in ["genomic"]:
                                if not dcGenomic.get(str(transcript)):
                                    dcGenomic.update({str(transcript): [line]})
                                else:
                                    dcGenomic.update(
                                        {
                                            str(transcript): dcGenomic.get(
                                                str(transcript)
                                            )
                                            + [line]
                                        }
                                    )
                            elif fields[2] in ["splice_junction"]:
                                if not dcSpliceJunctions.get(str(transcript)):
                                    dcSpliceJunctions.update({str(transcript): [line]})
                                else:
                                    dcSpliceJunctions.update(
                                        {
                                            str(transcript): dcSpliceJunctions.get(
                                                str(transcript)
                                            )
                                            + [line]
                                        }
                                    )
                            elif fields[2] in ["protein"]:
                                if not dcProt.get(str(transcript)):
                                    dcProt.update({str(transcript): [line]})
                                else:
                                    dcProt.update(
                                        {
                                            str(transcript): dcProt.get(str(transcript))
                                            + [line]
                                        }
                                    )
                        # Transcript Information
                        elif fields[1] == "TranscriptAttributes":
                            if not dcTranscriptAttributes.get(str(transcript)):
                                dcTranscriptAttributes.update({str(transcript): [line]})
                            else:
                                dcTranscriptAttributes.update(
                                    {
                                        str(transcript): dcTranscriptAttributes.get(
                                            str(transcript)
                                        )
                                        + [line]
                                    }
                                )
                        # Feature information
                        else:
                            if text[-1].endswith("T\n"):
                                if not dcTransFeatures.get(str(transcript)):
                                    dcTransFeatures.update({str(transcript): [line]})
                                else:
                                    dcTransFeatures.update(
                                        {
                                            str(transcript): dcTransFeatures.get(
                                                str(transcript)
                                            )
                                            + [line]
                                        }
                                    )
                            elif (
                                text[-1].endswith("P\n")
                                or text[-1].endswith("G\n")
                                or text[-1].endswith("N\n")
                            ):
                                if not dcProtFeatures.get(str(transcript)):
                                    dcProtFeatures.update({str(transcript): [line]})
                                else:
                                    dcProtFeatures.update(
                                        {
                                            str(transcript): dcProtFeatures.get(
                                                str(transcript)
                                            )
                                            + [line]
                                        }
                                    )

    return (
        dcTrans,
        dcExon,
        dcTransFeatures,
        dcGenomic,
        dcSpliceJunctions,
        dcProt,
        dcProtFeatures,
        dcTranscriptAttributes,
    )


def generateFinalGFF3(
    dcTrans,
    dcExon,
    dcTransFeatures,
    dcGenomic,
    dcSpliceJunctions,
    dcProt,
    dcProtFeatures,
    dcTranscriptAttributes,
    filename,
):
    # open new file
    res = open(filename, "w")
    strand = ""
    for SQtrans in dcTrans.keys():
        t = dcTrans.get(SQtrans)
        strand = t[0].split("\t")
        strand = strand[6]
        if t:
            for line in t:
                res.write(line)

        tf = dcTransFeatures.get(SQtrans)
        if tf:
            for line in tf:
                res.write(line)

        g = dcGenomic.get(SQtrans)
        if g:
            for line in g:
                res.write(line)

        e = dcExon.get(SQtrans)
        if e:
            if strand == "+":
                for line in e:
                    res.write(line)
            else:
                for i in range(len(e) - 1, -1, -1):
                    res.write(e[i])

        sj = dcSpliceJunctions.get(SQtrans)
        if sj:
            for line in sj:
                res.write(line)

        p = dcProt.get(SQtrans)
        if p:
            for line in p:
                res.write(line)

        pf = dcProtFeatures.get(SQtrans)
        if pf:
            for line in pf:
                res.write(line)

        ta = dcTranscriptAttributes.get(SQtrans)
        if ta:
            for line in ta:
                res.write(line)

    res.close()


############
# Parámetros
############

# -GTF de SQANTI3
# -Classification de SQANTI3
# -Junctions de SQANTI3
# -GFF3 de referencia


def main():
    global USE_GFF3
    global version
    # arguments
    parser = argparse.ArgumentParser(
        description="IsoAnnotLite "
        + str(version)
        + ": Transform SQANTI 3 output files to generate GFF3 to tappAS."
    )
    parser.add_argument(
        "corrected", help="\t\t*_corrected.gtf file from SQANTI 3 output."
    )
    parser.add_argument(
        "classification", help="\t\t*_classification.txt file from SQANTI 3 output."
    )
    parser.add_argument(
        "junction", help="\t\t*_junctions.txt file from SQANTI 3 output."
    )
    parser.add_argument(
        "-gff3",
        help="\t\ttappAS GFF3 file to map its annotation to your SQANTI 3 data (only if you use the same reference genome in SQANTI 3).",
        required=False,
    )

    args = parser.parse_args()

    # path and prefix for output files
    args.corrected = os.path.abspath(args.corrected)
    if not os.path.isfile(args.corrected):
        sys.stderr.write("ERROR: '%s' doesn't exist\n" % (args.corrected))
        sys.exit()

    args.classification = os.path.abspath(args.classification)
    if not os.path.isfile(args.classification):
        sys.stderr.write("ERROR: '%s' doesn't exist\n" % (args.classification))
        sys.exit()

    args.junction = os.path.abspath(args.junction)
    if not os.path.isfile(args.junction):
        sys.stderr.write("ERROR: '%s' doesn't exist\n" % (args.junction))
        sys.exit()

    if args.gff3:
        USE_GFF3 = True
        args.gff3 = os.path.abspath(args.gff3)
        if not os.path.isfile(args.gff3):
            sys.stderr.write("ERROR: '%s' doesn't exist\n" % (args.gff3))
            sys.exit()

    # Running functionality
    sys.stdout.write("\n\nRunning IsoAnnot Lite " + str(version) + "...\n")
    run(args)


def run(args):
    import time

    global USE_GFF3

    t1 = time.time()
    # corrected = input("Enter your file name for \"corrected.gtf\" file from SQANTI 3 (with extension): ")
    gtf = args.corrected
    # classification = input("Enter your file name for \"classification.txt\" file from SQANTI 3 (with extension): ")
    classification = args.classification
    # junctions = input("Enter your file name for \"junctions.txt\" file from SQANTI 3 (with extension): ")
    junctions = args.junction
    # GFF3 download from tappAS.org/downloads

    ########################
    # MAPPING SQANTI FILES #
    ########################

    if USE_GFF3:
        gff3 = args.gff3

        # File names
        filename = "tappAS_annot_from_SQANTI3.gff3"
        filenameMod = filename[:-5] + "_mod" + filename[-5:]

        #################
        # START PROCESS #
        #################
        print("\nReading SQANTI 3 Files and creating an auxiliar GFF...")

        # dc_SQexons = {trans : [[start,end], [start,end]...]}
        # dc_SQcoding = {trans : [CDSstart, CDSend, orf]}
        # dc_SQtransGene = {trans : [gene, category, transAssociated]}
        dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQstrand = createGTFFromSqanti(
            gtf, classification, junctions, filename
        )

        print("Reading reference annotation file and creating data variables...")
        # dc_GFF3 = {trans : [[start,end,line], [start,end,line], ...]}
        # dc_GFF3exonsTrans = {start : [trans, trans, ...]}
        # dc_GFF3transExons = {trans : [[start,end], [start,end]...]}
        # dc_GFF3coding = {trans : [CDSstart, CDSend]}
        dc_GFF3, dc_GFF3exonsTrans, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand = readGFF(
            gff3
        )  # dc_GFF3exons is sorted

        print("Transforming CDS local positions to genomic position...")
        # Transformar características a posiciones genómicas //revisar
        dc_SQcoding = transformCDStoGenomic(dc_SQcoding, dc_SQexons, dc_SQstrand)
        dc_GFF3coding = transformCDStoGenomic(
            dc_GFF3coding, dc_GFF3transExons, dc_GFF3strand
        )

        print("Transforming feature local positions to genomic position in GFF3...")
        # Transformar características a posiciones genómicas //revisar
        dc_GFF3_Genomic = transformTransFeaturesToGenomic(
            dc_GFF3, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand
        )

        print("Mapping transcript features betweeen GFFs...")
        mappingFeatures(
            dc_SQexons,
            dc_SQcoding,
            dc_SQtransGene,
            dc_GFF3exonsTrans,
            dc_GFF3transExons,
            dc_GFF3_Genomic,
            dc_GFF3coding,
            filename,
        )  # edit tappAS_annotation_from_Sqanti file

        print("Adding extra information to GFF3 columns...")
        updateGTF(filename, filenameMod)

        print("Reading GFF3 to sort it correctly...")
        dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes = readGFFandGetData(
            filenameMod
        )

        # Remove old files
        os.remove(filename)
        os.remove(filenameMod)

        dcTransFeatures = transformTransFeaturesToLocale(dcTransFeatures, dc_SQexons)

        print("Generating final GFF3...")
        generateFinalGFF3(
            dcTrans,
            dcExon,
            dcTransFeatures,
            dcGenomic,
            dcSpliceJunctions,
            dcProt,
            dcProtFeatures,
            dcTranscriptAttributes,
            filename,
        )

        t2 = time.time()
        time = t2 - t1
        print("Time used to generate new GFF3: " + "%.2f" % time + " seconds.\n")

        print("Exportation complete.\nYour GFF3 result is: '" + filename + "'\n")

    #####################
    # JUST SQANTI FILES #
    #####################

    else:
        # File names
        filename = "tappAS_annotation_from_SQANTI3.gff3"
        filenameMod = filename[:-5] + "_mod" + filename[-5:]

        #################
        # START PROCESS #
        #################
        print("\nReading SQANTI 3 Files and creating an auxiliar GFF...")

        # dc_SQexons = {trans : [[start,end], [start,end]...]}
        # dc_SQcoding = {trans : [CDSstart, CDSend, orf]}
        # dc_SQtransGene = {trans : [gene, category, transAssociated]}
        dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQstrand = createGTFFromSqanti(
            gtf, classification, junctions, filename
        )

        print("Adding extra information to relative columns...")
        updateGTF(filename, filenameMod)

        print("Reading GFF3 to sort it correctly...")
        dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes = readGFFandGetData(
            filenameMod
        )

        # Remove old files
        os.remove(filename)
        os.remove(filenameMod)

        print("Generating final GFF3...")
        generateFinalGFF3(
            dcTrans,
            dcExon,
            dcTransFeatures,
            dcGenomic,
            dcSpliceJunctions,
            dcProt,
            dcProtFeatures,
            dcTranscriptAttributes,
            filename,
        )

        t2 = time.time()
        time = t2 - t1
        print("Time used to generate new GFF3: " + "%.2f" % time + " seconds.\n")

        print("Exportation complete.\nYour GFF3 result is: '" + filename + "'\n")


if __name__ == "__main__":
    main()
