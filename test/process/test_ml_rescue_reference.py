"""
Process-level tests for run_randomforest_on_reference.R (ML rescue).
"""
import os
import sys
import shutil
import subprocess
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.commands import RESCUE_RANDOM_FOREST

REFERENCE_CLASSIFICATION = os.path.join(main_path, "test/test_data/isoforms/test_isoforms_classification.tsv")

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def train_forest(rds_file, columns):
    """Train a minimal random forest on random data with the given columns and save it like the ML filter does."""
    cols = ", ".join(f"`{c}` = runif(20)" for c in columns)
    code = (
        "set.seed(1); "
        f"x <- data.frame({cols}, check.names = FALSE); "
        "y <- factor(rep(c('POS', 'NEG'), 10)); "
        "rf <- caret::train(x = x, y = y, method = 'rf', "
        "trControl = caret::trainControl(method = 'none', classProbs = TRUE), tuneGrid = data.frame(mtry = 1)); "
        f"saveRDS(rf, '{rds_file}')"
    )
    subprocess.run(["Rscript", "-e", code], check=True, capture_output=True)


def run_reference_forest(rds_file, out_dir):
    # First 10 isoforms only: the full test classification has rows read.table cannot parse
    ref_class = out_dir / "reference_classification.txt"
    with open(REFERENCE_CLASSIFICATION) as h:
        ref_class.write_text("".join(h.readlines()[:11]))
    return subprocess.run(
        ["Rscript", RESCUE_RANDOM_FOREST, "-c", str(ref_class),
         "-o", "test", "-d", str(out_dir), "-r", str(rds_file)],
        capture_output=True, text=True,
    )


def test_reference_forest_runs_with_reference_columns(tmp_path):
    rds = tmp_path / "test_randomforest.RData"
    train_forest(rds, ["length", "exons", "sd_cov"])
    result = run_reference_forest(rds, tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "test_reference_isoform_predict.tsv").exists()


def test_reference_forest_reports_multisample_columns(tmp_path):
    """A forest trained on multi-sample counts uses columns that reference transcripts do not have.
    The script must say which ones and what to do, not fail with 'undefined columns selected'."""
    rds = tmp_path / "test_randomforest.RData"
    train_forest(rds, ["length", "exons", "FL.s1", "prevalence"])
    result = run_reference_forest(rds, tmp_path)
    assert result.returncode != 0
    assert "not in the reference classification: FL.s1, prevalence" in result.stderr
    assert "--remove_columns" in result.stderr
    assert "undefined columns selected" not in result.stderr
    assert not (tmp_path / "test_reference_isoform_predict.tsv").exists()
