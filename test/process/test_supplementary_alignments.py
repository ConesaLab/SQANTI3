"""
Process tests for supplementary/secondary alignments in mapping mode (issue #216).
SQANTI3 must use only the primary alignment of each isoform, keep the GTF and corrected FASTA IDs
consistent, not count indels of split alignments, and report the non-primary alignments.
"""
import os
import re
import shutil
from collections import defaultdict

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from src.helpers import report_non_primary_alignments, sequence_correction
from src.utilities.cupcake.io.BioReaders import GMAPSAMReader
from src.utilities.cupcake.sequence.err_correct_w_genome import err_correct
from src.utilities.cupcake.sequence.sam_to_gff3 import convert_sam_to_gff3
from src.utilities.indels_annot import calc_indels_from_sam

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
GENOME = os.path.join(main_path, "test/test_data/genome/genome_test.fasta")
REF_GTF = os.path.join(main_path, "test/test_data/reference/test_reference.gtf")

# PB.1.1 is split: primary + supplementary (with an insertion) + secondary.
# PB.2.1 is a normal isoform with one deletion.
SAM = """@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr22\tLN:50818468
PB.1.1\t0\tchr22\t20000001\t60\t100M500N100M100S\t*\t0\t0\t*\t*
PB.1.1\t2048\tchr22\t40000001\t60\t200H50M2I48M\t*\t0\t0\t*\t*
PB.1.1\t256\tchr22\t45000001\t0\t100M200S\t*\t0\t0\t*\t*
PB.2.1\t16\tchr22\t30000001\t60\t70M1D80M\t*\t0\t0\t*\t*
"""

@pytest.fixture(scope="module")
def genome_dict():
    return {r.id: r for r in SeqIO.parse(GENOME, "fasta")}

@pytest.fixture
def sam_file(tmp_path):
    sam = tmp_path / "test_corrected.sam"
    sam.write_text(SAM)
    return str(sam)

def fasta_ids(fasta):
    return [r.id for r in SeqIO.parse(fasta, "fasta")]

def gff_transcript_ids(gff):
    return [re.search(r"ID=([^;]+)", l).group(1) for l in open(gff) if "\tgene\t" in l]

### Primary-only reading ###

def test_reader_default_returns_all_alignments(sam_file):
    assert [r.qID for r in GMAPSAMReader(sam_file, True)] == ["PB.1.1", "PB.1.1", "PB.1.1", "PB.2.1"]

def test_reader_skip_non_primary(sam_file):
    recs = list(GMAPSAMReader(sam_file, True, skip_non_primary=True))
    assert [r.qID for r in recs] == ["PB.1.1", "PB.2.1"]
    assert recs[0].sStart == 20000000

def test_corrected_fasta_and_gtf_ids_match(sam_file, genome_dict, tmp_path):
    fasta = str(tmp_path / "corrected.fasta")
    gff = str(tmp_path / "corrected.gff3")
    err_correct(GENOME, sam_file, fasta, genome_dict=genome_dict)
    convert_sam_to_gff3(sam_file, gff, source="test")

    assert fasta_ids(fasta) == ["PB.1.1", "PB.2.1"]
    assert gff_transcript_ids(gff) == ["PB.1.1", "PB.2.1"]

def test_indels_from_supplementary_not_counted(sam_file):
    _, indelsTotal = calc_indels_from_sam(sam_file)
    assert indelsTotal["PB.1.1"] == 0  # its only indel is in the supplementary alignment
    assert indelsTotal["PB.2.1"] == 1

### Report of non-primary alignments ###

def test_report_lists_all_alignments_of_split_isoforms(sam_file, tmp_path):
    out = tmp_path / "supplementary.txt"
    assert report_non_primary_alignments(sam_file, str(out)) == 1

    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    assert rows[0] == ["isoform", "alignment", "chrom", "strand", "start", "end", "MAPQ"]
    assert rows[1:] == [
        ["PB.1.1", "primary", "chr22", "+", "20000001", "20000700", "60"],
        ["PB.1.1", "supplementary", "chr22", "+", "40000001", "40000098", "60"],
        ["PB.1.1", "secondary", "chr22", "+", "45000001", "45000100", "0"],
    ]

def test_report_not_written_without_split_isoforms(tmp_path):
    sam = tmp_path / "primary_only.sam"
    sam.write_text("\n".join(l for l in SAM.splitlines() if not re.search(r"\t(2048|256)\t", l)) + "\n")
    out = tmp_path / "supplementary.txt"
    assert report_non_primary_alignments(str(sam), str(out)) == 0
    assert not out.exists()

### End-to-end with minimap2 ###

def spliced_sequence(genome_dict, exons, strand):
    seq = Seq("".join(str(genome_dict["chr22"].seq[s - 1:e]) for s, e in sorted(exons)))
    return seq if strand == "+" else seq.reverse_complement()

@pytest.mark.skipif(shutil.which("minimap2") is None or shutil.which("gffread") is None,
                    reason="minimap2 and gffread are required")
def test_sequence_correction_chimeric_isoform(genome_dict, tmp_path):
    exons, strands = defaultdict(list), {}
    for line in open(REF_GTF):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "exon":
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            exons[tid].append((int(f[3]), int(f[4])))
            strands[tid] = f[6]
    # Chimera of two transcripts ~39 Mb apart, plus a normal control
    chimera = (spliced_sequence(genome_dict, exons["ENST00000657645.1"], strands["ENST00000657645.1"]) +
               spliced_sequence(genome_dict, exons["ENST00000496652.5"], strands["ENST00000496652.5"]))
    control = spliced_sequence(genome_dict, exons["ENST00000397906.6"], strands["ENST00000397906.6"])
    isoforms = tmp_path / "isoforms.fasta"
    isoforms.write_text(f">PB.1.1\n{chimera}\n>PB.2.1\n{control}\n")

    outdir = tmp_path / "out"
    outdir.mkdir()
    sequence_correction(str(outdir), "test", cpus=2, chunks=1, fasta=True, genome_dict=genome_dict,
                        badstrandGTF=str(outdir / "unknown_strand.gtf"), genome=GENOME,
                        isoforms=str(isoforms), aligner_choice="minimap2")

    sam_flags = [int(l.split("\t")[1]) for l in open(outdir / "test_corrected.sam") if not l.startswith("@")]
    assert any(flag & 2048 for flag in sam_flags), "minimap2 did not split the chimera; test is not informative"

    assert fasta_ids(outdir / "test_corrected.fasta") == ["PB.1.1", "PB.2.1"]
    gtf_ids = {re.search(r'transcript_id "([^"]+)"', l).group(1) for l in open(outdir / "test_corrected.gtf")}
    assert gtf_ids == {"PB.1.1", "PB.2.1"}

    report = [l.split("\t") for l in open(outdir / "test_supplementary_alignments.txt")][1:]
    assert {r[0] for r in report} == {"PB.1.1"}
    assert sorted(r[1] for r in report) == ["primary", "supplementary"]
