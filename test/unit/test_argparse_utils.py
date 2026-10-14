import logging
import os, sys
import pytest
import pysam

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.argparse_utils import valid_unique_ids, index_sr_bams

@pytest.fixture
def tester_logger():
    return logging.getLogger("tester_logger")

def test_valid_unique_ids_ok(tmp_path, tester_logger):
    fasta = tmp_path / "isoforms.fasta"
    fasta.write_text(">PB.1.1|a\nACGT\n>PB.2.1|a\nACGT\n")
    assert valid_unique_ids(str(fasta), tester_logger) == str(fasta)

def test_valid_unique_ids_aborts_on_collision(tmp_path, tester_logger, caplog):
    fasta = tmp_path / "isoforms.fasta"
    fasta.write_text(">PB.1.1|copyA\nACGT\n>PB.1.1|copyB\nACGT\n>PB.2.1\nACGT\n")
    caplog.set_level(logging.ERROR)
    with pytest.raises(SystemExit):
        valid_unique_ids(str(fasta), tester_logger)
    assert "1 sequence IDs" in caplog.text
    assert "PB.1.1 <- PB.1.1|copyA, PB.1.1|copyB" in caplog.text

def write_sorted_bam(bam_path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        read = pysam.AlignedSegment(bam.header)
        read.query_name, read.reference_id, read.reference_start = "r1", 0, 100
        read.cigarstring, read.query_sequence = "50M", "A" * 50
        bam.write(read)

def test_index_sr_bams_directory(tmp_path, tester_logger):
    """Every BAM of an --SR_bam directory gets an index in index_dir, none next to the input (issue #611)."""
    bam_dir = tmp_path / "bams"
    bam_dir.mkdir()
    for name in ["rep1.bam", "rep2.bam"]:
        write_sorted_bam(bam_dir / name)
    index_dir = tmp_path / "bam_index"
    index_sr_bams(str(bam_dir), str(index_dir), 1, tester_logger)
    assert len(os.listdir(index_dir)) == 2
    assert sorted(os.listdir(bam_dir)) == ["rep1.bam", "rep2.bam"]

def test_index_sr_bams_missing_bam_in_fofn(tmp_path, tester_logger, caplog):
    fofn = tmp_path / "bams.fofn"
    fofn.write_text(f"{tmp_path}/missing.bam\n")
    caplog.set_level(logging.ERROR)
    with pytest.raises(SystemExit):
        index_sr_bams(str(fofn), str(tmp_path / "bam_index"), 1, tester_logger)
    assert "missing.bam not found" in caplog.text
