"""Unit tests for B11/B12/B13 — orthogonal-support views (CAGE 5', polyA 3',
short-read). Validated on synthetic data because the chr22 fixture has all
support columns 100% NaN (the views correctly no-op there).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_pdf import PdfPages

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities.sqanti_reads_plots import (
    _summarize_support,
    _column_is_populated,
    load_sqanti_file,
    compute_cage_metrics,
    compute_polya_metrics,
    compute_sr_metrics,
    plot_support_pages,
    assemble_scorecard_metrics,
)
from src.utilities.sqanti_reads_report import _metric_cohort_figure


def _class_df():
    return pd.DataFrame({
        "within_CAGE_peak": ["TRUE", "FALSE", "TRUE", "TRUE"],       # 3/4 = 75%
        "dist_to_CAGE_peak": [10, -20, 0, 5],                        # |.| median 7.5
        "within_polyA_site": ["TRUE", "TRUE", "FALSE", "FALSE"],     # 50%
        "polyA_motif_found": ["TRUE", "FALSE", "FALSE", "FALSE"],    # 25%
        "dist_to_polyA_site": [1, -3, 2, 4],                         # |.| median 2.5
        "ratio_TSS": [1.0, 2.0, 3.0, 4.0],                          # median 2.5
    })


def _jxn_df():
    return pd.DataFrame({"total_coverage_unique": [0, 5, 10, 0]})    # >0 -> 50%


def test_summarize_support_populated():
    row = _summarize_support(_class_df(), _jxn_df(), "S", "temp_factor", 0).iloc[0]
    assert row["perc_within_cage"] == pytest.approx(75.0)
    assert row["median_dist_cage"] == pytest.approx(7.5)
    assert row["perc_within_polya"] == pytest.approx(50.0)
    assert row["perc_polya_motif"] == pytest.approx(25.0)
    assert row["median_dist_polya"] == pytest.approx(2.5)
    assert row["median_ratio_TSS"] == pytest.approx(2.5)
    assert row["perc_jxn_SR_supported"] == pytest.approx(50.0)


def test_summarize_support_absent_columns_all_nan():
    # No support columns at all -> every scalar NaN (the fixture's situation).
    empty_class = pd.DataFrame({"structural_category": ["FSM", "ISM"]})
    empty_jxn = pd.DataFrame({"canonical": ["canonical", "canonical"]})
    row = _summarize_support(empty_class, empty_jxn, "S", "temp_factor", 0).iloc[0]
    for c in ("perc_within_cage", "median_dist_cage", "perc_within_polya",
              "perc_polya_motif", "median_dist_polya", "median_ratio_TSS",
              "perc_jxn_SR_supported"):
        assert pd.isna(row[c])


def test_compute_metrics_populated_and_empty():
    pop = _summarize_support(_class_df(), _jxn_df(), "S", "temp_factor", 0)
    assert compute_cage_metrics(pop)["samples"] == ["S"]
    assert compute_polya_metrics(pop)["samples"] == ["S"]
    assert compute_sr_metrics(pop)["samples"] == ["S"]
    # all-NaN support -> every view no-ops
    nan_row = pd.DataFrame([{"sampleID": "S", "temp_factor": 0,
                             **{c: np.nan for c in
                                ("perc_within_cage", "median_dist_cage", "perc_within_polya",
                                 "perc_polya_motif", "median_dist_polya", "median_ratio_TSS",
                                 "perc_jxn_SR_supported")}}])
    assert compute_cage_metrics(nan_row)["samples"] == []
    assert compute_polya_metrics(nan_row)["samples"] == []
    assert compute_sr_metrics(nan_row)["samples"] == []
    assert compute_cage_metrics(None)["samples"] == []


def test_pages_and_figures(tmp_path):
    # Two samples so the cohort bar has peers.
    s1 = _summarize_support(_class_df(), _jxn_df(), "A", "temp_factor", 0)
    s2 = _summarize_support(_class_df(), _jxn_df(), "B", "temp_factor", 0)
    support = pd.concat([s1, s2], ignore_index=True)
    cage, polya, sr = (compute_cage_metrics(support),
                       compute_polya_metrics(support), compute_sr_metrics(support))
    assert _metric_cohort_figure(cage["per_sample"], "perc_within_cage", "t", "y") is not None
    assert _metric_cohort_figure(polya["per_sample"], "perc_within_polya", "t", "y") is not None
    assert _metric_cohort_figure(sr["per_sample"], "perc_jxn_SR_supported", "t", "y") is not None
    p = tmp_path / "b.pdf"
    with PdfPages(str(p)) as pdf:
        plot_support_pages(pdf, cage, polya, sr)
    pypdf = pytest.importorskip("pypdf")
    assert len(pypdf.PdfReader(str(p)).pages) == 4      # CAGE + polyA + 2 SR
    # empty metrics -> no pages
    p2 = tmp_path / "empty.pdf"
    empty = {"samples": [], "per_sample": pd.DataFrame(columns=["sampleID"])}
    with PdfPages(str(p2)) as pdf:
        plot_support_pages(pdf, empty, empty, empty)
    assert not os.path.exists(str(p2)) or len(pypdf.PdfReader(str(p2)).pages) == 0


def test_scorecard_feed_includes_support():
    support = pd.concat([
        _summarize_support(_class_df(), _jxn_df(), "A", "temp_factor", 0),
        _summarize_support(_class_df(), _jxn_df(), "B", "temp_factor", 0),
    ], ignore_index=True)
    rows = assemble_scorecard_metrics(["A", "B"], support_DF=support)
    d = {r["sampleID"]: r for r in rows}
    assert d["A"]["perc_within_cage"] == pytest.approx(75.0)
    assert d["A"]["perc_within_polya"] == pytest.approx(50.0)
    assert d["A"]["perc_jxn_SR_supported"] == pytest.approx(50.0)
    # absent support -> keys present but None (scorecard drops them automatically)
    rows2 = assemble_scorecard_metrics(["A"], support_DF=None)
    assert rows2[0]["perc_within_cage"] is None


def test_tolerant_load_skips_missing_columns(tmp_path):
    f = tmp_path / "class.txt"
    # File has isoform + within_CAGE_peak but NOT within_polyA_site.
    f.write_text("isoform\twithin_CAGE_peak\nr1\tTRUE\nr2\tFALSE\n")
    cols = ["isoform", "within_CAGE_peak", "within_polyA_site"]
    dtypes = {"isoform": "string", "within_CAGE_peak": "string",
              "within_polyA_site": "string"}
    # Non-tolerant would raise on the missing column; tolerant loads what exists.
    df = load_sqanti_file(str(f), cols, dtypes, tolerant=True)
    assert "within_CAGE_peak" in df.columns and "within_polyA_site" not in df.columns
    assert _column_is_populated(df, "within_CAGE_peak")
    assert not _column_is_populated(df, "within_polyA_site")
