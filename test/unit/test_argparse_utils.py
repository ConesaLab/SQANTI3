import gzip
import logging
import os, sys
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.argparse_utils import valid_unique_ids, valid_fasta, contains_gene_name

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
