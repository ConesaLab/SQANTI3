"""Unit tests for the splice-site fuzziness (junction-offset) analysis added to
``sqanti_reads_plots``: strand-aware site assignment, the collapsed offset
summary, the derived metrics (spectrum / precision profile / per-sample scalar /
canonical split), and empty-input safety of both the metric and the PDF page
renderer.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities.sqanti_reads_plots import (
    _summarize_jxn_offsets,
    compute_jxn_offset_metrics,
    plot_jxn_offset_pages,
)

OFFSET_COLS = ["sampleID", "group", "site", "offset", "canonical", "count"]


def _jxn(strand, d_start, d_end, canonical="canonical"):
    """One junction row with the columns _summarize_jxn_offsets reads."""
    return {
        "strand": strand,
        "diff_to_Ref_start_site": d_start,
        "diff_to_Ref_end_site": d_end,
        "canonical": canonical,
    }


def test_strand_aware_site_assignment():
    """On + the start diff is the donor and the end diff the acceptor; on - it is
    reversed. Offsets within the window are kept and labelled accordingly."""
    df = pd.DataFrame([
        _jxn("+", 3, -5),   # donor +3, acceptor -5
        _jxn("-", 7, -2),   # acceptor +7, donor -2  (reversed)
    ])
    out = _summarize_jxn_offsets(df, "S1", "group", "g", window=15)
    got = {(r.site, r.offset) for r in out.itertuples()}
    assert ("donor", 3) in got
    assert ("acceptor", -5) in got
    assert ("acceptor", 7) in got
    assert ("donor", -2) in got
    assert int(out["count"].sum()) == 4  # two junctions -> four site observations


def test_window_excludes_far_offsets():
    """Offsets beyond +/- window are treated as novel sites and dropped here."""
    df = pd.DataFrame([_jxn("+", 3, 100)])   # acceptor 100 bp away -> excluded
    out = _summarize_jxn_offsets(df, "S1", "group", "g", window=15)
    assert set(zip(out["site"], out["offset"])) == {("donor", 3)}


def test_nan_diffs_are_dropped():
    df = pd.DataFrame([_jxn("+", np.nan, np.nan), _jxn("+", 0, 0)])
    out = _summarize_jxn_offsets(df, "S1", "group", "g", window=15)
    assert int(out["count"].sum()) == 2  # only the second junction's two sites


def test_empty_input_returns_typed_frame():
    out = _summarize_jxn_offsets(pd.DataFrame(), "S1", "group", "g")
    assert list(out.columns) == OFFSET_COLS
    assert out.empty


def _toy_offsets():
    """Two samples with a known offset composition.

    S1: 8 exact + 2 imprecise (offsets +4, +4)  -> 20% imprecise, all downstream
    S2: 9 exact + 1 imprecise (offset -3)        -> 10% imprecise, upstream
    """
    rows = [
        {"sampleID": "S1", "group": "g", "site": "donor", "offset": 0, "canonical": "canonical", "count": 8},
        {"sampleID": "S1", "group": "g", "site": "donor", "offset": 4, "canonical": "non_canonical", "count": 2},
        {"sampleID": "S2", "group": "g", "site": "donor", "offset": 0, "canonical": "canonical", "count": 9},
        {"sampleID": "S2", "group": "g", "site": "donor", "offset": -3, "canonical": "canonical", "count": 1},
    ]
    return pd.DataFrame(rows)


def test_per_sample_scalars():
    m = compute_jxn_offset_metrics(_toy_offsets(), window=15)
    ps = m["per_sample"].set_index("sampleID")
    assert ps.at["S1", "perc_imprecise"] == pytest.approx(20.0)
    assert ps.at["S1", "perc_exact"] == pytest.approx(80.0)
    assert ps.at["S1", "median_abs_offset"] == pytest.approx(4.0)
    assert ps.at["S1", "perc_downstream"] == pytest.approx(100.0)
    assert ps.at["S2", "perc_imprecise"] == pytest.approx(10.0)
    assert ps.at["S2", "perc_downstream"] == pytest.approx(0.0)  # the one fuzzy site is upstream


def test_precision_profile_monotone_and_bounds():
    m = compute_jxn_offset_metrics(_toy_offsets(), window=15)
    prof = m["profile"]
    for s in m["samples"]:
        d = prof[prof["sampleID"] == s].sort_values("k")
        vals = d["perc_within"].to_numpy()
        assert np.all(np.diff(vals) >= -1e-9)          # non-decreasing in k
        assert d[d["k"] == 15]["perc_within"].iloc[0] == pytest.approx(100.0)
    # k=0 equals the exact-match rate
    s1_k0 = prof[(prof["sampleID"] == "S1") & (prof["k"] == 0)]["perc_within"].iloc[0]
    assert s1_k0 == pytest.approx(80.0)


def test_by_class_split():
    m = compute_jxn_offset_metrics(_toy_offsets(), window=15)
    bc = m["by_class"].set_index(["sampleID", "canonical"])
    # S1 non_canonical bucket is entirely the two fuzzy observations
    assert bc.at[("S1", "non_canonical"), "perc_imprecise"] == pytest.approx(100.0)
    assert bc.at[("S1", "canonical"), "perc_imprecise"] == pytest.approx(0.0)


def test_empty_metrics_and_render_noop(tmp_path):
    from matplotlib.backends.backend_pdf import PdfPages
    empty = compute_jxn_offset_metrics(
        pd.DataFrame(columns=OFFSET_COLS), window=15)
    assert empty["samples"] == []
    # renderer must no-op cleanly on empty metrics and on None
    out = tmp_path / "empty.pdf"
    with PdfPages(str(out)) as pdf:
        plot_jxn_offset_pages(pdf, empty)
        plot_jxn_offset_pages(pdf, None)
