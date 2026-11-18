"""Unit tests for A6 — pairwise UJC overlap index (|A∩B|/|A|).

Validated on synthetic presence matrices (the real fixture's three samples share
almost no exact junction chains, so its overlap is ~0 — correct but degenerate).
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

from src.utilities.sqanti_reads_plots import compute_ujc_overlap, plot_ujc_overlap_page
from src.utilities.sqanti_reads_report import _ujc_overlap_figure


def _upset(cols, index):
    """Build a boolean presence matrix (jxnHash rows × sample cols) + SC column."""
    data = {s: vals for s, vals in cols.items()}
    df = pd.DataFrame(data, index=index).astype(bool)
    df["structural_category"] = ["novel_in_catalog"] * len(index)
    return df


def test_asymmetry_and_ratio():
    # A={h1,h2,h3}, B={h2,h3,h4,h5}: shared=2, |A|=3, |B|=4.
    upset = _upset({"A": [1, 1, 1, 0, 0], "B": [0, 1, 1, 1, 1]},
                   ["h1", "h2", "h3", "h4", "h5"])
    m = compute_ujc_overlap({"samples": ["A", "B"], "upset": upset})["matrix"]
    assert m.loc["A", "B"] == pytest.approx(2 / 3)   # 2 of A's 3 are in B
    assert m.loc["B", "A"] == pytest.approx(2 / 4)   # 2 of B's 4 are in A
    assert m.loc["A", "A"] == 1.0 and m.loc["B", "B"] == 1.0  # diagonal


def test_identical_and_disjoint():
    upset = _upset({"X": [1, 1, 0], "Y": [1, 1, 0], "Z": [0, 0, 1]},
                   ["1", "2", "3"])
    m = compute_ujc_overlap({"samples": ["X", "Y", "Z"], "upset": upset})["matrix"]
    assert m.loc["X", "Y"] == 1.0 and m.loc["Y", "X"] == 1.0  # identical repertoires
    assert m.loc["X", "Z"] == 0.0 and m.loc["Z", "X"] == 0.0  # disjoint


def test_empty_is_noop():
    upset = _upset({"X": [1, 0]}, ["1", "2"])
    ov = compute_ujc_overlap({"samples": ["X"], "upset": upset})
    assert ov["matrix"] is None and ov["samples"] == []
    assert compute_ujc_overlap(None)["matrix"] is None
    assert compute_ujc_overlap({})["matrix"] is None
    assert _ujc_overlap_figure(ov) is None


def test_page_and_figure_build(tmp_path):
    upset = _upset({"A": [1, 1, 1, 0], "B": [0, 1, 1, 1]}, ["h1", "h2", "h3", "h4"])
    ov = compute_ujc_overlap({"samples": ["A", "B"], "upset": upset})
    assert _ujc_overlap_figure(ov) is not None
    p = tmp_path / "ov.pdf"
    with PdfPages(str(p)) as pdf:
        plot_ujc_overlap_page(pdf, ov)
    pypdf = pytest.importorskip("pypdf")
    assert len(pypdf.PdfReader(str(p)).pages) == 1
    # empty overlap draws no page
    p2 = tmp_path / "empty.pdf"
    with PdfPages(str(p2)) as pdf:
        plot_ujc_overlap_page(pdf, {"samples": [], "matrix": None})
    assert not os.path.exists(str(p2)) or len(pypdf.PdfReader(str(p2)).pages) == 0
