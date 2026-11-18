"""Unit tests for compute_qc_pca — the sample-similarity PCA built over the same
per-sample metric matrix the cohort-outlier scorecard scores (plus transcript
divergency), with constant/absent features dropped.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities.sqanti_reads_plots import compute_qc_pca


def _sc(n=4):
    rows = []
    for i, s in enumerate([f"S{k+1}" for k in range(n)]):
        rows.append({"sampleID": s,
                     "median_length": 1000 + i * 300, "perc_ISM": 10 + i * 5,
                     "composition_drift": 0.1 * i, "jxn_per_1k_reads": 50 - i * 3,
                     "perc_reads_RTS": 2.0,          # constant -> dropped
                     "perc_within_cage": None})      # absent  -> dropped
    return rows


def test_uses_scorecard_metrics_and_drops_constant_absent():
    fmap = {"S1": "a", "S2": "a", "S3": "b", "S4": "b"}
    td = {"S1": 0.2, "S2": 0.25, "S3": 0.3, "S4": 0.35}
    pca_DF, loadings, var = compute_qc_pca(_sc(), fmap, "grp",
                                           extra={"transcript_divergency": td})
    # PC columns present, factor mapped
    assert 0 in pca_DF.columns and 1 in pca_DF.columns
    assert pca_DF["grp"].tolist() == ["a", "a", "b", "b"]
    feats = set(map(str, loadings.index))
    # named QC metrics feed the PCA (interpretable loadings), TD folded in
    assert {"median_length", "perc_ISM", "composition_drift",
            "jxn_per_1k_reads", "transcript_divergency"} <= feats
    # constant + all-absent metrics dropped
    assert "perc_reads_RTS" not in feats and "perc_within_cage" not in feats
    assert len(var) >= 2


def test_empty_and_single_sample_safe():
    # no metrics -> trivial frames
    pca_DF, loadings, var = compute_qc_pca([], {}, "temp_factor")
    assert list(pca_DF.columns) == ["sampleID", "temp_factor"] and loadings.empty and len(var) == 0
    # a single sample -> no PC1 (the renderers' PCA guard then skips the plot)
    one = [{"sampleID": "S1", "median_length": 1000, "perc_ISM": 10}]
    pca_DF, loadings, var = compute_qc_pca(one, {"S1": 0}, "temp_factor")
    assert 1 not in pca_DF.columns
