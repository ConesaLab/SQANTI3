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


class TestValidateMODS:
    """Test suite for validate_MODS."""

    def base_args(self, **kwargs):
        args = SimpleNamespace(MODPROB=[0.9], MODTYPES=["m6A"], sqanti_dirs=None, SKIPHASH=False,
                               aligner_choice="minimap2", chunks=1)
        args.__dict__.update(kwargs)
        return args

    def test_other_aligner_is_rejected(self, caplog):
        """Mods mode always aligns with minimap2, so asking for uLTRA stops the run."""
        args = self.base_args(aligner_choice="uLTRA")
        with pytest.raises(SystemExit):
            validate_MODS(args, logger)
        assert "--aligner_choice uLTRA" in caplog.text

    def test_chunks_are_reported_as_ignored(self, caplog):
        """The alignment is done once before QC, so --chunks only gives a warning."""
        args = self.base_args(chunks=4)
        with caplog.at_level(logging.WARNING):
            validate_MODS(args, logger)
        assert "--chunks is ignored" in caplog.text


class TestEmptySamples:
    """A sample without usable reads or without modifications must not stop a multi-sample run."""

    def test_sample_without_mapped_reads(self, tmp_path):
        write_mods_bam(tmp_path, "s1", "s1", [])
        write_classification(tmp_path, "s1", "s1", [])
        args = mods_args(tmp_path)
        gen_mod_tsv(args, "s1", "s1", logger)
        out = ujc_mod_sites(args, "s1", "s1", logger)
        assert len(pd.read_csv(out, sep="\t")) == 0
        assert len(pd.read_csv(tmp_path / "s1" / "s1_ujc_mod_sites.tsv", sep="\t")) == 0

    def test_sample_without_modifications_above_threshold(self, tmp_path):
        write_mods_bam(tmp_path, "s1", "s1", [{"name": "r1", "ml": (10,), "start": 999}])
        write_classification(tmp_path, "s1", "s1", [["r1", "full-splice_match", "chr1_+_monoexon_T1", "h1"]])
        args = mods_args(tmp_path)
        gen_mod_tsv(args, "s1", "s1", logger)
        out = ujc_mod_sites(args, "s1", "s1", logger)
        classif = pd.read_csv(out, sep="\t")
        assert classif["isoform"].tolist() == ["r1"]
        assert classif["m6A_n_mod"].tolist() == [0]
        assert len(pd.read_csv(tmp_path / "s1" / "s1_ujc_mod_sites.tsv", sep="\t")) == 0


class TestGenModTsv:
    """Test suite for gen_mod_tsv."""

    def test_columns_of_the_per_read_table(self, tmp_path):
        write_mods_bam(tmp_path, "s1", "s1", [{"name": "r1", "start": 999}])
        table = pd.read_csv(gen_mod_tsv(mods_args(tmp_path, modtypes=("m6A", "m5C")), "s1", "s1", logger), sep="\t")
        assert table.columns.tolist() == ["isoform", "read_length", "read_start", "read_end", "chrom", "strand",
                                          "m6A_positions", "m6A_probs", "m6A_n_mod",
                                          "m5C_positions", "m5C_probs", "m5C_n_mod"]
        assert table[["isoform", "read_start", "read_end", "m6A_n_mod", "m5C_n_mod"]].values.tolist() == [["r1", 1000, 1099, 1, 0]]

    def test_sample_without_reads_keeps_the_header(self, tmp_path, caplog):
        write_mods_bam(tmp_path, "s1", "s1", [])
        table = pd.read_csv(gen_mod_tsv(mods_args(tmp_path), "s1", "s1", logger), sep="\t")
        assert len(table) == 0
        assert "m6A_n_mod" in table.columns
        assert "no mapped reads with modification tags" in caplog.text
