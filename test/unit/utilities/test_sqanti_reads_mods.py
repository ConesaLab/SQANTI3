import array
import glob
import logging
import os
import sys
from types import SimpleNamespace

import pandas as pd
import pysam
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities import sqanti_reads_mods as mods
from src.utilities.sqanti_reads_mods import (
    check_uBAM, gen_mod_tsv, mod_files_check, process_read, ujc_mod_sites, validate_MODS
)

logger = logging.getLogger("test_mods")

SEQ = "ACGT" * 25  # 100 nt, an A every 4 bases starting at read position 0


def make_read(header, name, mm="A+a,0;", ml=(255,), start=None):
    """A read with MM/ML tags. Unmapped unless start (0-based) is given."""
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = SEQ
    read.query_qualities = pysam.qualitystring_to_array("I" * len(SEQ))
    if start is None:
        read.flag = 4
    else:
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigartuples = [(0, len(SEQ))]
    read.set_tag("MM", mm, "Z")
    read.set_tag("ML", array.array("B", ml))
    return read


def write_ubam(path, reads_kwargs):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"}})
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for kwargs in reads_kwargs:
            bam.write(make_read(header, **kwargs))


def write_mods_bam(out_dir, file_acc, sampleID, reads_kwargs):
    """Sorted and indexed <sampleID>_mods.bam, as align_uBAM leaves it."""
    sample_dir = out_dir / file_acc
    sample_dir.mkdir(parents=True, exist_ok=True)
    path = str(sample_dir / f"{sampleID}_mods.bam")
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chr1", "LN": 1000000}]})
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for kwargs in sorted(reads_kwargs, key=lambda k: k["start"]):
            bam.write(make_read(header, **kwargs))
    pysam.index(path)
    return path


def write_classification(out_dir, file_acc, sampleID, rows):
    path = out_dir / file_acc / f"{sampleID}_reads_classification.txt"
    pd.DataFrame(rows, columns=["isoform", "structural_category", "jxn_string", "jxnHash"]).to_csv(
        path, sep="\t", index=False)
    return path


def mods_args(out_dir, modtypes=("m6A",), prob=0.9):
    return SimpleNamespace(OUTPUT=str(out_dir), MODTYPES=list(modtypes),
                           MODTOTHRESH={mod: prob for mod in modtypes})


class TestModFilesCheck:
    """Test suite for mod_files_check."""

    def test_sample_name_that_is_prefix_of_another(self, tmp_path, monkeypatch):
        """sample1 must not be run with the uBAM of sample10, whatever order glob returns."""
        for name in ("sample1.bam", "sample10.bam", "genome.fasta", "ref.gtf"):
            (tmp_path / name).write_text("")
        real_glob = glob.glob
        monkeypatch.setattr(mods.glob, "glob", lambda pattern: sorted(real_glob(pattern), reverse=True))
        monkeypatch.setattr(mods, "check_uBAM", lambda path, reads_logger: {"m6A"})
        args = SimpleNamespace(input_dir=str(tmp_path), refFasta=str(tmp_path / "genome.fasta"),
                               refGTF=str(tmp_path / "ref.gtf"), verbose=False, MODTYPES=["m6A"])
        _, uBAM = mod_files_check(args, "sample1", logger)
        assert os.path.basename(uBAM) == "sample1.bam"

    def test_file_with_suffix(self, tmp_path, monkeypatch):
        for name in ("sample1_mods.bam", "genome.fasta", "ref.gtf"):
            (tmp_path / name).write_text("")
        monkeypatch.setattr(mods, "check_uBAM", lambda path, reads_logger: {"m6A"})
        args = SimpleNamespace(input_dir=str(tmp_path), refFasta=str(tmp_path / "genome.fasta"),
                               refGTF=str(tmp_path / "ref.gtf"), verbose=False, MODTYPES=["m6A"])
        _, uBAM = mod_files_check(args, "sample1", logger)
        assert os.path.basename(uBAM) == "sample1_mods.bam"

    def test_several_matches_and_none_exact(self, tmp_path, caplog):
        """The uBAM is never guessed: the error lists the candidates."""
        for name in ("sample1_a.bam", "sample10.bam"):
            (tmp_path / name).write_text("")
        args = SimpleNamespace(input_dir=str(tmp_path))
        with pytest.raises(SystemExit):
            mod_files_check(args, "sample1", logger)
        assert "sample10.bam, sample1_a.bam" in caplog.text
