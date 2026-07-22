"""Unit tests for the per-sample cohort-context artefact pages (RT-switching and
intra-priming) that accompany the outlier scorecard: the PDF page renderer and
the HTML figure builder, including their no-op guards on missing data.
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
    _plot_metric_cohort_page,
    plot_artefact_metric_pages,
)
from src.utilities.sqanti_reads_report import _metric_cohort_figure

QC = {
    "perc_reads_RTS": {"warn": 15.0, "fail": 30.0, "worse": "high"},
    "perc_reads_intrapriming": {"warn": 20.0, "fail": 40.0, "worse": "high"},
}


def _err_df():
    return pd.DataFrame({
        "sampleID": ["S0", "S1", "S2", "S3", "S4"],
        "perc_reads_RTS": [16.0, 16.2, 16.4, 16.6, 41.0],
        "perc_reads_intrapriming": [20.0, 20.1, 20.2, 20.3, 50.0],
    })


def _scorecard():
    # S4 flagged fail on both artefact metrics; peers pass.
    cell = {s: {"perc_reads_RTS": "pass", "perc_reads_intrapriming": "pass"}
            for s in ["S0", "S1", "S2", "S3"]}
    cell["S4"] = {"perc_reads_RTS": "fail", "perc_reads_intrapriming": "fail"}
    return {"enabled": True,
            "metrics": ["perc_reads_RTS", "perc_reads_intrapriming"],
            "cell_flags": cell}


def test_artefact_pages_render_two_pages(tmp_path):
    out = tmp_path / "art.pdf"
    with PdfPages(str(out)) as pdf:
        plot_artefact_metric_pages(pdf, _err_df(), QC, scorecard=_scorecard())
    import PyPDF2
    assert len(PyPDF2.PdfReader(str(out)).pages) == 2


def test_artefact_pages_no_scorecard_still_renders(tmp_path):
    out = tmp_path / "art.pdf"
    with PdfPages(str(out)) as pdf:
        plot_artefact_metric_pages(pdf, _err_df(), QC, scorecard=None)
    import PyPDF2
    assert len(PyPDF2.PdfReader(str(out)).pages) == 2


def test_artefact_pages_noop_on_missing(tmp_path):
    out = tmp_path / "art.pdf"
    with PdfPages(str(out)) as pdf:
        plot_artefact_metric_pages(pdf, None, QC)                       # no df
        plot_artefact_metric_pages(pdf, pd.DataFrame({"sampleID": ["x"]}), QC)  # no metric cols
        _plot_metric_cohort_page(pdf, _err_df(), "does_not_exist",
                                 title="t", ylabel="y")                 # missing metric
    # No page was drawn, so PdfPages never flushed a file to disk. The guards
    # holding (no exception, no output) is the assertion.
    assert not out.exists()


def test_html_metric_figure_builds_and_flags():
    fig = _metric_cohort_figure(_err_df(), "perc_reads_RTS", "RTS", "%",
                                scorecard=_scorecard(), threshold=QC["perc_reads_RTS"])
    assert fig is not None
    # one bar trace with 5 samples
    bar = fig.data[0]
    assert list(bar.x) == ["S0", "S1", "S2", "S3", "S4"]
    # S4 bar is flag-colored (fail red), peers are the base color
    colors = list(bar.marker.color)
    assert colors[4] == "#F44336"
    assert colors[0] != "#F44336"


def test_html_metric_figure_missing_returns_none():
    assert _metric_cohort_figure(_err_df(), "nope", "x", "y") is None
    assert _metric_cohort_figure(None, "perc_reads_RTS", "x", "y") is None
