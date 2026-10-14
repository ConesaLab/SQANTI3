"""Unit tests for A5 — multi-axis replicate concordance (within-group agreement
on structural composition, length profile, and splice-site imprecision).
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
    compute_replicate_concordance,
    plot_replicate_concordance_page,
)
from src.utilities.sqanti_reads_report import _replicate_concordance_figure


def _inputs(comp_rows, len_rows, impr_rows):
    pivot = pd.DataFrame([{"sampleID": s, **c} for s, c in comp_rows.items()])
    length = pd.DataFrame([{"sampleID": s, "q25_length": q, "median_length": m,
                            "q75_length": q3} for s, (q, m, q3) in len_rows.items()])
    jxn = {"samples": list(impr_rows),
           "per_sample": pd.DataFrame([{"sampleID": s, "perc_imprecise": v}
                                       for s, v in impr_rows.items()])}
    return pivot, length, jxn


def test_consistent_group_scores_high():
    comp = {s: {"FSM": 60, "ISM": 40} for s in ("A", "B", "C")}
    ln = {s: (1000, 2000, 3000) for s in ("A", "B", "C")}
    impr = {s: 5.0 for s in ("A", "B", "C")}
    piv, length, jxn = _inputs(comp, ln, impr)
    fmap = {"A": "g1", "B": "g1", "C": "g1"}
    m = compute_replicate_concordance(piv, length, jxn, fmap)
    assert m["enabled"] is True
    ps = m["per_sample"].set_index("sampleID")
    for ax in m["axes"]:
        assert ps.loc["A", ax] == pytest.approx(1.0, abs=1e-6)


def test_shifted_rep_diverges_on_composition():
    comp = {"A": {"FSM": 60, "ISM": 40}, "B": {"FSM": 60, "ISM": 40},
            "ODD": {"FSM": 10, "ISM": 90}}
    ln = {s: (1000, 2000, 3000) for s in ("A", "B", "ODD")}
    impr = {s: 5.0 for s in ("A", "B", "ODD")}
    piv, length, jxn = _inputs(comp, ln, impr)
    fmap = {"A": "g", "B": "g", "ODD": "g"}
    ps = compute_replicate_concordance(piv, length, jxn, fmap)["per_sample"].set_index("sampleID")
    # ODD agrees less on composition than the two consistent reps do
    assert ps.loc["ODD", "category_composition"] < ps.loc["A", "category_composition"]
    # ...but still matches on length + imprecision
    assert ps.loc["ODD", "length_profile"] == pytest.approx(1.0, abs=1e-6)
    assert ps.loc["ODD", "imprecision"] == pytest.approx(1.0, abs=1e-6)


def test_disabled_reasons(tmp_path):
    comp = {s: {"FSM": 50, "ISM": 50} for s in ("A", "B")}
    ln = {s: (1, 2, 3) for s in ("A", "B")}
    impr = {s: 1.0 for s in ("A", "B")}
    piv, length, jxn = _inputs(comp, ln, impr)
    # no factor
    m0 = compute_replicate_concordance(piv, length, jxn, None)
    assert m0["enabled"] is False and "factor" in m0["reason"]
    # all singletons
    m1 = compute_replicate_concordance(piv, length, jxn, {"A": "g1", "B": "g2"})
    assert m1["enabled"] is False and "replicate" in m1["reason"]
    # disabled -> no page, no figure
    assert _replicate_concordance_figure(m1) is None
    p = tmp_path / "r.pdf"
    with PdfPages(str(p)) as pdf:
        plot_replicate_concordance_page(pdf, m1)
    pypdf = pytest.importorskip("pypdf")
    assert not os.path.exists(str(p)) or len(pypdf.PdfReader(str(p)).pages) == 0


def test_mad_zero_for_identical_group():
    # A8: identical replicates -> zero within-group spread on every axis.
    comp = {s: {"FSM": 60, "ISM": 40} for s in ("A", "B", "C")}
    ln = {s: (1000, 2000, 3000) for s in ("A", "B", "C")}
    impr = {s: 5.0 for s in ("A", "B", "C")}
    piv, length, jxn = _inputs(comp, ln, impr)
    m = compute_replicate_concordance(piv, length, jxn,
                                      {"A": "g", "B": "g", "C": "g"})
    ps = m["per_sample"].set_index("sampleID")
    for col in m["mad_axes"]:
        assert ps.loc["A", col] == pytest.approx(0.0, abs=1e-9)


def test_mad_tracks_group_spread():
    # A8: a group spread across length + imprecision (but sharing one category mix)
    # -> non-zero MAD on those axes, zero on composition. MAD is robust: with three
    # differing values the median absolute deviation is the middle deviation.
    comp = {s: {"FSM": 60, "ISM": 40} for s in ("A", "B", "C")}
    ln = {"A": (1000, 2000, 3000), "B": (1000, 2200, 3000), "C": (1000, 2400, 3000)}
    impr = {"A": 4.0, "B": 6.0, "C": 8.0}
    piv, length, jxn = _inputs(comp, ln, impr)
    m = compute_replicate_concordance(piv, length, jxn,
                                      {"A": "g", "B": "g", "C": "g"})
    ps = m["per_sample"].set_index("sampleID")
    assert ps.loc["A", "mad_category_composition"] == pytest.approx(0.0, abs=1e-9)
    assert ps.loc["A", "mad_length_profile"] == pytest.approx(200.0)   # median|{200,0,200}|
    assert ps.loc["A", "mad_imprecision"] == pytest.approx(2.0)        # median|{2,0,2}|
    # MAD is a group-level quantity: identical for every member of the group.
    assert ps.loc["A", "mad_length_profile"] == pytest.approx(ps.loc["C", "mad_length_profile"])


def test_mad_nonzero_for_two_replicate_group():
    # A8 value-add: MAD stays informative for exactly two replicates, where the
    # pairwise agreement bars are identical by construction.
    comp = {"A": {"FSM": 60, "ISM": 40}, "B": {"FSM": 60, "ISM": 40}}
    ln = {"A": (1000, 2000, 3000), "B": (1000, 2400, 3000)}
    impr = {"A": 4.0, "B": 8.0}
    piv, length, jxn = _inputs(comp, ln, impr)
    m = compute_replicate_concordance(piv, length, jxn, {"A": "g", "B": "g"})
    ps = m["per_sample"].set_index("sampleID")
    assert ps.loc["A", "mad_length_profile"] == pytest.approx(200.0)   # |2000-2400|/2
    assert ps.loc["A", "mad_imprecision"] == pytest.approx(2.0)        # |4-8|/2


def test_fixture_like_groupB_singleton_excluded():
    comp = {"SQ_R1": {"FSM": 50, "ISM": 50}, "SQ_R2": {"FSM": 55, "ISM": 45},
            "SQ_R3": {"FSM": 40, "ISM": 60}}
    ln = {s: (1000, 2000, 3000) for s in comp}
    impr = {s: 5.0 for s in comp}
    piv, length, jxn = _inputs(comp, ln, impr)
    fmap = {"SQ_R1": "groupA", "SQ_R2": "groupA", "SQ_R3": "groupB"}
    m = compute_replicate_concordance(piv, length, jxn, fmap)
    assert m["enabled"] is True
    scored = set(m["per_sample"]["sampleID"])
    assert scored == {"SQ_R1", "SQ_R2"}          # groupB singleton excluded
    assert _replicate_concordance_figure(m) is not None
