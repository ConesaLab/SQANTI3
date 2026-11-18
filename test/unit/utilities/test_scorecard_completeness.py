"""Unit tests for the cohort-relative sample-outlier scorecard and the 5'/3'
read-end completeness analysis added to ``sqanti_reads_plots``.

The scorecard tests exercise the three regimes the tool must survive on *any*
dataset: all-samples-clean (no false alarms), one degraded outlier (flagged on
concordant metrics), and a whole-cohort shift that is internally consistent
(stays quiet — the case a fixed absolute threshold would get wrong). The
completeness tests cover the strand/end digest, the ECDF monotonicity, the
within-window scalars, and empty-input safety.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities.sqanti_reads_plots import (
    _robust_z,
    _summarize_completeness,
    compute_completeness_metrics,
    compute_sample_scorecard,
    plot_scorecard_page,
    plot_completeness_pages,
)
from src.utilities.sqanti_reads_config import load_config

COMP_COLS = ["sampleID", "group", "end", "abs_dist", "count"]


# --------------------------------------------------------------------------- #
# robust z
# --------------------------------------------------------------------------- #
def test_robust_z_flags_single_outlier():
    # a realistic cohort: small spread among peers, one clear outlier. With a
    # non-zero MAD the outlier lands well beyond the fail threshold.
    z = _robust_z([10, 11, 9, 10, 40])
    assert np.all(np.abs(z[:4]) < 2)       # peers stay small
    assert z[4] > 3.5                       # the outlier is far out


def test_robust_z_mad_zero_falls_back_to_std():
    # majority identical -> MAD=0 -> std fallback (outlier inflates std, so its z
    # is attenuated but still the largest). Must not divide-by-zero or NaN.
    z = _robust_z([10, 10, 10, 10, 40])
    assert np.isfinite(z).all()
    assert z[4] == max(z)


def test_robust_z_all_equal_is_zero():
    assert np.allclose(_robust_z([5, 5, 5, 5]), 0.0)


# --------------------------------------------------------------------------- #
# completeness summarizer + metrics
# --------------------------------------------------------------------------- #
def _class_df(tss, tts):
    return pd.DataFrame({"diff_to_gene_TSS": tss, "diff_to_gene_TTS": tts})


def test_completeness_summarize_shapes_and_clip():
    df = _class_df([0, 10, -20, 100000], [0, -5, 3, 2])
    out = _summarize_completeness(df, "S1", "group", "g", window=50)
    assert list(out.columns) == COMP_COLS
    assert set(out["end"]) == {"5prime", "3prime"}
    # the 100000 bp read is clipped to the 10*window ceiling (500)
    assert out["abs_dist"].max() == 500


def test_completeness_missing_columns_returns_empty():
    out = _summarize_completeness(pd.DataFrame({"length": [1, 2]}), "S1", "group", "g")
    assert out.empty and list(out.columns) == COMP_COLS


def test_completeness_metrics_within_window_and_profile():
    # 6 reads within 50 bp at 5', 4 beyond -> 60% complete
    tss = [0, 10, 20, 30, 40, 50, 200, 300, 400, 500]
    tts = [0] * 10
    df = _class_df(tss, tts)
    comp = _summarize_completeness(df, "S1", "group", "g", window=50)
    m = compute_completeness_metrics(comp, window=50)
    ps = m["per_sample"].set_index("sampleID")
    assert ps.at["S1", "perc_5p_within_window"] == pytest.approx(60.0)
    assert ps.at["S1", "perc_3p_within_window"] == pytest.approx(100.0)
    # profile is non-decreasing in k for each end
    prof = m["profile"]
    for end in ("5prime", "3prime"):
        v = prof[prof["end"] == end].sort_values("k")["perc_within"].to_numpy()
        assert np.all(np.diff(v) >= -1e-9)


def test_completeness_metrics_empty():
    m = compute_completeness_metrics(pd.DataFrame(columns=COMP_COLS), window=50)
    assert m["samples"] == []
    assert m["per_sample"].empty


# --------------------------------------------------------------------------- #
# scorecard — the three cohort regimes
# --------------------------------------------------------------------------- #
def _cohort(n=5, outlier_idx=None):
    rows = []
    for i in range(n):
        r = {"sampleID": f"S{i}", "median_length": 2000 + 10 * i,
             "perc_ISM": 25 + 0.5 * i, "perc_novel_junctions": 6 + 0.1 * i,
             "perc_5p_within_window": 58 + i, "perc_3p_within_window": 70 + i,
             "perc_reads_RTS": 16 + 0.2 * i, "perc_reads_intrapriming": 20 + 0.1 * i,
             "perc_sites_imprecise": 2.5 + 0.05 * i}
        rows.append(r)
    if outlier_idx is not None:
        rows[outlier_idx].update(median_length=1100, perc_ISM=38,
                                 perc_novel_junctions=9.5, perc_5p_within_window=46)
    return rows


@pytest.fixture
def cfg():
    return load_config()


def test_scorecard_disabled_below_min_samples(cfg):
    sc = compute_sample_scorecard(_cohort(3), cfg)
    assert sc["enabled"] is False
    assert "3" in sc["reason"]


def test_scorecard_clean_cohort_is_quiet(cfg):
    sc = compute_sample_scorecard(_cohort(5), cfg)
    assert sc["enabled"] is True
    assert set(sc["overall"].values()) == {"pass"}


def test_scorecard_flags_single_outlier(cfg):
    sc = compute_sample_scorecard(_cohort(5, outlier_idx=4), cfg)
    assert sc["overall"]["S4"] == "fail"
    assert all(sc["overall"][f"S{i}"] == "pass" for i in range(4))
    assert sc["n_flagged"]["S4"] >= 2


def test_scorecard_consistent_shift_stays_quiet(cfg):
    # whole cohort short but internally consistent -> no flags (relative scoring)
    rows = _cohort(5)
    for r in rows:
        r["median_length"] = 300 + 5 * int(r["sampleID"][1])
    sc = compute_sample_scorecard(rows, cfg)
    assert set(sc["overall"].values()) == {"pass"}


def test_scorecard_good_side_never_flagged(cfg):
    # one sample much BETTER (longer, higher completeness) must not be flagged
    rows = _cohort(5)
    rows[2].update(median_length=6000, perc_5p_within_window=95,
                   perc_ISM=10, perc_novel_junctions=2)
    sc = compute_sample_scorecard(rows, cfg)
    assert sc["overall"]["S2"] == "pass"


def test_scorecard_drops_absent_metric(cfg):
    # imprecision absent for all -> silently dropped, others still scored
    rows = _cohort(5, outlier_idx=4)
    for r in rows:
        r["perc_sites_imprecise"] = None
    sc = compute_sample_scorecard(rows, cfg)
    assert "perc_sites_imprecise" not in sc["metrics"]
    assert sc["overall"]["S4"] == "fail"


def test_scorecard_render_noop_when_disabled(tmp_path):
    from matplotlib.backends.backend_pdf import PdfPages
    out = tmp_path / "sc.pdf"
    with PdfPages(str(out)) as pdf:
        plot_scorecard_page(pdf, {"enabled": False})
        plot_scorecard_page(pdf, None)
        plot_completeness_pages(pdf, {"samples": []})
