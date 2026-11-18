"""Unit tests for the two derived read-quality scalars and their cohort-context
pages: the ISM truncation-fragment fraction and the novel non-canonical junction
burden (B8 / B10). Covers the metric helpers, the PDF page renderer, and the
scorecard feed.
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
    _ism_fragment_pct,
    _novel_noncanonical_pct,
    plot_quality_metric_pages,
    assemble_scorecard_metrics,
)

QC = {
    "perc_ISM_fragments": {"warn": 60.0, "fail": 80.0, "worse": "high"},
    "perc_novel_noncanonical_jxn": {"warn": 2.0, "fail": 5.0, "worse": "high"},
}


def _ism_df():
    # S0: mostly fragments (90/100); S1: mostly mono-exon (10 frag/100)
    return pd.DataFrame({
        "sampleID": ["S0", "S1"],
        "3prime_fragment": [60.0, 5.0],
        "5prime_fragment": [20.0, 3.0],
        "internal_fragment": [10.0, 2.0],
        "mono-exon": [5.0, 80.0],
        "intron_retention": [5.0, 10.0],
    })


def _nov_df():
    return pd.DataFrame({
        "sampleID": ["S0", "S1"],
        "known_canonical": [900.0, 950.0],
        "known_non_canonical": [10.0, 5.0],
        "novel_canonical": [80.0, 40.0],
        "novel_non_canonical": [10.0, 5.0],
    })


# --------------------------------------------------------------------------- #
# ISM fragment fraction
# --------------------------------------------------------------------------- #
def test_ism_fragment_pct_values():
    d = _ism_fragment_pct(_ism_df())
    assert d["S0"] == pytest.approx(90.0)   # (60+20+10)/100
    assert d["S1"] == pytest.approx(10.0)   # (5+3+2)/100


def test_ism_fragment_pct_empty_and_missing():
    assert _ism_fragment_pct(None) == {}
    assert _ism_fragment_pct(pd.DataFrame({"other": [1]})) == {}
    # only the sampleID column -> no numeric subcols -> empty
    assert _ism_fragment_pct(pd.DataFrame({"sampleID": ["S0"]})) == {}


def test_ism_fragment_pct_missing_fragment_cols_is_zero():
    # a sample with no fragment columns at all -> 0% fragments, not a crash
    df = pd.DataFrame({"sampleID": ["S0"], "mono-exon": [50.0], "intron_retention": [50.0]})
    assert _ism_fragment_pct(df)["S0"] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# novel non-canonical burden
# --------------------------------------------------------------------------- #
def test_novel_noncanonical_pct_values():
    d = _novel_noncanonical_pct(_nov_df())
    assert d["S0"] == pytest.approx(10.0 / 1000.0 * 100)   # 1.0%
    assert d["S1"] == pytest.approx(5.0 / 1000.0 * 100)    # 0.5%


def test_novel_noncanonical_pct_missing_cols():
    assert _novel_noncanonical_pct(None) == {}
    assert _novel_noncanonical_pct(pd.DataFrame({"sampleID": ["S0"],
                                                 "known_canonical": [1]})) == {}


# --------------------------------------------------------------------------- #
# scorecard feed + page rendering
# --------------------------------------------------------------------------- #
def test_ism_fragments_has_no_absolute_threshold_but_is_scorecard_metric():
    # ISM fragment fraction is a shape signal, not a quality verdict: it must NOT
    # carry an absolute warn/fail cut-off, but it must still feed the cohort-
    # relative scorecard. (Novel non-canonical, being an artefact class, keeps its
    # absolute threshold.)
    from src.utilities.sqanti_reads_config import load_config
    cfg = load_config()
    assert "perc_ISM_fragments" not in cfg["qc_flags"]
    assert cfg["sample_scorecard"]["metrics"].get("perc_ISM_fragments") == "high"
    assert "perc_novel_noncanonical_jxn" in cfg["qc_flags"]


def test_scorecard_metrics_include_derived():
    rows = assemble_scorecard_metrics(["S0", "S1"], ISM_DF=_ism_df(), nov_can_DF=_nov_df())
    by = {r["sampleID"]: r for r in rows}
    assert by["S0"]["perc_ISM_fragments"] == pytest.approx(90.0)
    assert by["S0"]["perc_novel_noncanonical_jxn"] == pytest.approx(1.0)


def test_quality_pages_render_two_pages(tmp_path):
    out = tmp_path / "q.pdf"
    with PdfPages(str(out)) as pdf:
        plot_quality_metric_pages(pdf, _ism_df(), _nov_df(), QC, scorecard=None)
    import pypdf
    assert len(pypdf.PdfReader(str(out)).pages) == 2


def test_quality_pages_noop_on_missing(tmp_path):
    out = tmp_path / "q.pdf"
    with PdfPages(str(out)) as pdf:
        plot_quality_metric_pages(pdf, None, None, QC)
    assert not out.exists()
