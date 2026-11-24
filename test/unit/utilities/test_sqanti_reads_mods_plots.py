import os
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities import sqanti_reads_mods_plots as mods_plots
from src.utilities.sqanti_reads_mods_plots import add_ujc_positions, run_mods_plots


def reads_table(rows):
    return pd.DataFrame(rows, columns=["sampleID", "isoform", "read_start", "read_end", "strand",
                                       "jxn_string", "jxnHash"])


def sites_table(rows):
    return pd.DataFrame(rows, columns=["sampleID", "jxnHash", "mod_type", "mod_pos"])


class TestAddUjcPositions:
    """Test suite for add_ujc_positions."""

    def test_multiexon_site_position(self):
        """Exons 1000-1099 and 2000-2099 on the + strand: position 2009 is 110 nt into the UJC."""
        reads = reads_table([["s1", "r1", 1000, 2099, "+", "chr1_+_1100_1999", "h1"]])
        sites = sites_table([["s1", "h1", "m6A", 2009]])
        placed, _, _ = add_ujc_positions(sites, reads)
        assert placed["ujc_length"].tolist() == [200]
        assert placed["tx_pos"].tolist() == [109]
        assert placed["dist_3p"].tolist() == [90]
        assert placed["dist_junction"].tolist() == [10]

    @pytest.mark.xfail(reason="PR #646 issue 1: monoexon reads of different exons are joined")
    def test_monoexon_reads_of_different_exons_share_a_hash(self):
        """Monoexon reads assigned to a reference transcript are named after the transcript, so
        fragments in different exons share a UJC. A site must not be placed on an exon that spans
        the intron between them."""
        jxn = "chr1_+_monoexon_ENST0001_incomplete-splice_match"
        reads = reads_table([["s1", "r1", 1000, 1999, "+", jxn, "h1"],
                             ["s1", "r2", 201000, 201999, "+", jxn, "h1"]])
        sites = sites_table([["s1", "h1", "m6A", 1500]])
        placed, _, _ = add_ujc_positions(sites, reads)
        # 1000 nt is the longest read; a larger model means the two fragments were joined
        assert placed["ujc_length"].iloc[0] <= 1000
        assert placed["dist_3p"].iloc[0] <= 499


class TestRunModsPlots:
    """Test suite for run_mods_plots."""

    def test_no_sample_with_reads(self, tmp_path, caplog):
        """With every sample empty there is nothing to draw: one warning and no report."""
        classif = tmp_path / "s1_reads_mods_classification.txt"
        sites = tmp_path / "s1_ujc_mod_sites.tsv"
        classif.write_text("\t".join(mods_plots.READ_COLS + ["m6A_n_mod"]) + "\n")
        sites.write_text("jxnHash\tchrom\tstrand\tmod_type\tmod_pos\tn_mod\tmean_prob\tn_reads_ujc\tn_cov\toccupancy\n")
        design = tmp_path / "design.csv"
        design.write_text(f"sampleID,mod_classification_file,ujc_mod_sites_file\ns1,{classif},{sites}\n")
        mods_plots.reads_logger.addHandler(caplog.handler)
        try:
            result = run_mods_plots(str(design), str(tmp_path), "test", SimpleNamespace(MODTYPES=["m6A"]))
        finally:
            mods_plots.reads_logger.removeHandler(caplog.handler)
        assert result is None
        assert "the modification plots are skipped" in caplog.text
        assert not (tmp_path / "test_mods_plots.pdf").exists()
