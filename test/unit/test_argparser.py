import os, sys, logging
from unittest import mock

import pytest

main_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../..")
sys.path.append(main_path)

from src.argparse_utils import valid_gtf

GTF_LINE = (
    "chr1\ttest\tgene\t1\t100\t.\t+\t.\t"
    'gene_id "G1"; gene_name "G1";\n'
    "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\t"
    'gene_id "G1"; transcript_id "T1"; gene_name "G1";\n'
)


@pytest.fixture
def tester_logger():
    logger = logging.getLogger("tester_logger")
    logger.setLevel(logging.INFO)
    return logger


def _write(path, content=GTF_LINE):
    with open(path, "w") as handle:
        handle.write(content)
    return path


def test_valid_gtf_accepts_gff3(tmp_path, tester_logger):
    gff3 = _write(str(tmp_path / "reference.gff3"))
    converted = str(tmp_path / "reference.gtf3")

    def fake_call(cmd):
        # gffread writes the converted file in real usage
        _write(converted)
        return 0

    with mock.patch("src.argparse_utils.subprocess.call", side_effect=fake_call):
        result = valid_gtf(gff3, tester_logger)

    assert result == converted


def test_valid_gtf_accepts_gff(tmp_path, tester_logger):
    gff = _write(str(tmp_path / "reference.gff"))
    converted = str(tmp_path / "reference.gtf")

    def fake_call(cmd):
        _write(converted)
        return 0

    with mock.patch("src.argparse_utils.subprocess.call", side_effect=fake_call):
        result = valid_gtf(gff, tester_logger)

    assert result == converted


def test_valid_gtf_accepts_gtf_directly(tmp_path, tester_logger):
    gtf = _write(str(tmp_path / "reference.gtf"))
    with mock.patch("src.argparse_utils.subprocess.call") as mocked:
        result = valid_gtf(gtf, tester_logger)
    mocked.assert_not_called()
    assert result == gtf


def test_valid_gtf_rejects_unknown_extension(tmp_path, tester_logger):
    bad = _write(str(tmp_path / "reference.txt"))
    with pytest.raises(SystemExit):
        valid_gtf(bad, tester_logger)
