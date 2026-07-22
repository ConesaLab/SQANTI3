"""Unit tests for A4 — depth-normalised junction yield (junctions per 1000 reads).

Covers the arithmetic, empty-safety, the PDF page renderer, and the scorecard
feed (a cohort where one sample has a low depth-normalised yield is flagged).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import pytest
from matplotlib.backends.backend_pdf import PdfPages

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities.sqanti_reads_plots import (
    compute_jxn_yield_vs_depth,
    plot_jxn_yield_page,
    assemble_scorecard_metrics,
    compute_sample_scorecard,
)

CFG = {"sample_scorecard": {"min_samples": 4, "z_warn": 2.5, "z_fail": 3.5,
                            "min_metrics_warn": 1, "min_metrics_fail": 1,
                            "metrics": {"jxn_per_1k_reads": "low"}}}


def test_arithmetic():
    m = compute_jxn_yield_vs_depth({"A": 8000, "B": 3000}, {"A": 1000, "B": 1000})
    ps = m["per_sample"].set_index("sampleID")
    assert ps.loc["A", "jxn_per_1k_reads"] == pytest.approx(8000.0)
    assert ps.loc["B", "jxn_per_1k_reads"] == pytest.approx(3000.0)
    assert set(m["samples"]) == {"A", "B"}


def test_empty_safe():
    for m in (compute_jxn_yield_vs_depth({}, {}),
              compute_jxn_yield_vs_depth({"A": 10}, {}),
              compute_jxn_yield_vs_depth({"A": 10}, {"A": 0})):  # zero reads dropped
        assert m["samples"] == []
        assert list(m["per_sample"].columns) == ["sampleID", "n_jxn", "n_reads",
                                                 "jxn_per_1k_reads"]


def test_page_renders_and_noop(tmp_path):
    pypdf = pytest.importorskip("pypdf")
    m = compute_jxn_yield_vs_depth({"A": 8000, "B": 3000}, {"A": 1000, "B": 1000})
    p = tmp_path / "y.pdf"
    with PdfPages(str(p)) as pdf:
        plot_jxn_yield_page(pdf, m)                                    # 1 page
        plot_jxn_yield_page(pdf, compute_jxn_yield_vs_depth({}, {}))   # no-op, adds nothing
    assert len(pypdf.PdfReader(str(p)).pages) == 1


def test_scorecard_feed_flags_low_yield():
    # 5 samples, one (LOW) with a much lower depth-normalised yield -> flagged.
    jxn = {"S1": 8000, "S2": 8200, "S3": 7900, "S4": 8100, "LOW": 2000}
    reads = {s: 1000 for s in jxn}
    ym = compute_jxn_yield_vs_depth(jxn, reads)
    rows = assemble_scorecard_metrics(list(jxn), yield_metrics=ym)
    assert all("jxn_per_1k_reads" in r for r in rows)
    sc = compute_sample_scorecard(rows, CFG)
    assert sc["enabled"] is True
    assert sc["overall"]["LOW"] in ("warn", "fail")
    assert sc["overall"]["S1"] == "pass"
