import gzip
import logging
import os, sys
from argparse import Namespace
from unittest import mock
import pytest
import pysam

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.argparse_utils import valid_unique_ids, valid_fasta, contains_gene_name, valid_gtf, filter_args_validation, \
    index_sr_bams

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

def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)
    return str(path)

@pytest.mark.parametrize("name", ["isoforms.fasta.gz", "isoforms.fa.gz", "reads.fastq.gz", "reads.fq.gz"])
def test_valid_fasta_gzipped_allowed(tmp_path, tester_logger, name):
    path = write_gz(tmp_path / name, ">PB.1.1\nACGT\n")
    assert valid_fasta(path, tester_logger, allow_gz=True) == path

def test_valid_fasta_gzipped_not_allowed(tmp_path, tester_logger, caplog):
    """Inputs read with plain open() (genome, ORF input, filter/rescue FASTAs) must be uncompressed."""
    path = write_gz(tmp_path / "genome.fa.gz", ">chr1\nACGT\n")
    caplog.set_level(logging.ERROR)
    with pytest.raises(SystemExit):
        valid_fasta(path, tester_logger)
    assert "is gzipped" in caplog.text

def test_valid_fasta_gzipped_wrong_extension(tmp_path, tester_logger, caplog):
    path = write_gz(tmp_path / "isoforms.txt.gz", ">PB.1.1\nACGT\n")
    caplog.set_level(logging.ERROR)
    with pytest.raises(SystemExit):
        valid_fasta(path, tester_logger, allow_gz=True)
    assert "is not a FASTA file" in caplog.text

def test_contains_gene_name_gzipped(tmp_path):
    assert contains_gene_name(write_gz(tmp_path / "with.fasta.gz", ">PB.1.1 gene_name=X\nACGT\n"))
    assert not contains_gene_name(write_gz(tmp_path / "without.fasta.gz", ">PB.1.1\nACGT\n"))

GTF_LINES = (
    'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
)

def fake_gffread(content=GTF_LINES):
    """Stands in for run_command: writes content to the -o path of the gffread command."""
    def run(cmd, *args, **kwargs):
        with open(cmd.split()[-1], "w") as out:
            out.write(content)
    return run

class TestValidGtf:
    """Test suite for valid_gtf."""

    @pytest.mark.parametrize("name", ["reference.gff", "reference.gff3"])
    def test_gff_is_converted_into_out_dir(self, tmp_path, tester_logger, name):
        """The converted GTF goes to out_dir only: the input folder may be read-only (issue #611)."""
        in_dir, out_dir = tmp_path / "input", tmp_path / "out"
        in_dir.mkdir(); out_dir.mkdir()
        gff = in_dir / name
        gff.write_text("##gff-version 3\n")
        with mock.patch("src.argparse_utils.run_command", side_effect=fake_gffread()) as gffread:
            result = valid_gtf(str(gff), tester_logger, str(out_dir))
        assert gffread.call_count == 1
        assert result == str(out_dir / "reference.converted.gtf")
        assert sorted(os.listdir(in_dir)) == [name]

    def test_gtf_is_used_directly(self, tmp_path, tester_logger):
        gtf = tmp_path / "reference.gtf"
        gtf.write_text(GTF_LINES)
        with mock.patch("src.argparse_utils.run_command") as gffread:
            result = valid_gtf(str(gtf), tester_logger, str(tmp_path / "out"))
        gffread.assert_not_called()
        assert result == str(gtf)

    def test_unknown_extension_aborts(self, tmp_path, tester_logger, caplog):
        bad = tmp_path / "reference.txt"
        bad.write_text(GTF_LINES)
        caplog.set_level(logging.ERROR)
        with pytest.raises(SystemExit):
            valid_gtf(str(bad), tester_logger, str(tmp_path))
        assert "is not a GTF file" in caplog.text

    def test_empty_conversion_aborts(self, tmp_path, tester_logger, caplog):
        """gffread exits without error on a file that is not GFF, so an empty result must stop the run."""
        gff = tmp_path / "reference.gff3"
        gff.write_text("not a gff\n")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        caplog.set_level(logging.ERROR)
        with mock.patch("src.argparse_utils.run_command", side_effect=fake_gffread("")):
            with pytest.raises(SystemExit):
                valid_gtf(str(gff), tester_logger, str(out_dir))
        assert "Please check that the input is a valid GFF file" in caplog.text

def test_filter_args_validation_uses_converted_gtf(tmp_path):
    """The caller must keep the path valid_gtf returns, or the run goes on with the GFF."""
    classification = tmp_path / "classification.txt"
    classification.write_text("isoform\n")
    json_filter = tmp_path / "filter.json"
    json_filter.write_text("{}")
    gff = tmp_path / "isoforms.gff3"
    gff.write_text("##gff-version 3\n")
    out_dir = tmp_path / "out"
    args = Namespace(sqanti_class=str(classification), filter_isoforms=None, filter_gtf=str(gff),
                     filter_sam=None, filter_faa=None, dir=str(out_dir), subcommand="rules",
                     json_filter=str(json_filter))
    with mock.patch("src.argparse_utils.run_command", side_effect=fake_gffread()):
        filter_args_validation(args)
    assert args.filter_gtf == str(out_dir / "isoforms.converted.gtf")
    assert args.filter_gtf_input == str(gff)

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
