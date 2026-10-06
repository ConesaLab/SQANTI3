import logging
import os, sys
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.argparse_utils import valid_unique_ids

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
