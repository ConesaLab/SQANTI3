import os
import sys
from types import SimpleNamespace

import pandas as pd

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities import sqanti_reads_mods_plots as mods_plots
from src.utilities.sqanti_reads_mods_plots import run_mods_plots


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
