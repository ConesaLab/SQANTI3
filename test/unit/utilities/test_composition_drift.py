"""Unit tests for A3 — composition drift (pairwise Jensen–Shannon distance on
the structural-category composition vectors). Covers the JSD arithmetic, empty-
safety, the PDF page + HTML figure no-op, and the scorecard feed.
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
    compute_composition_drift,
    plot_composition_drift_page,
    assemble_scorecard_metrics,
    compute_sample_scorecard,
)
from src.utilities.sqanti_reads_report import _composition_drift_figure


def _pivot(rows):
    """rows: {sampleID: {cat: pct}} -> pivot DataFrame with sampleID + cat cols."""
    return pd.DataFrame([{"sampleID": s, **cats} for s, cats in rows.items()]).fillna(0)


def test_identical_rows_zero_distance():
    df = _pivot({"A": {"FSM": 50, "ISM": 50}, "B": {"FSM": 50, "ISM": 50}})
    m = compute_composition_drift(df)
    assert m["matrix"].loc["A", "B"] == pytest.approx(0.0, abs=1e-9)
    assert m["per_sample"].set_index("sampleID").loc["A", "mean_jsd"] == pytest.approx(0.0, abs=1e-9)


def test_disjoint_rows_max_distance():
    df = _pivot({"A": {"FSM": 100, "NNC": 0}, "B": {"FSM": 0, "NNC": 100}})
    m = compute_composition_drift(df)
    # JS distance (base 2) between disjoint distributions is 1.0.
    assert m["matrix"].loc["A", "B"] == pytest.approx(1.0, abs=1e-6)


def test_symmetric_diagonal_zero():
    df = _pivot({"A": {"FSM": 70, "ISM": 30}, "B": {"FSM": 40, "ISM": 60},
                 "C": {"FSM": 55, "ISM": 45}})
    m = compute_composition_drift(df)
    mat = m["matrix"]
    assert list(np.diag(mat.values)) == [0.0, 0.0, 0.0]
    assert np.allclose(mat.values, mat.values.T)
    assert len(m["per_sample"]) == 3


def test_empty_safe(tmp_path):
    for df in (None,
               _pivot({"A": {"FSM": 100}}),                 # single sample
               pd.DataFrame({"sampleID": ["A", "B"]})):      # no category cols
        m = compute_composition_drift(df)
        assert m["samples"] == []
        assert _composition_drift_figure(m) is None
        p = tmp_path / "d.pdf"
        with PdfPages(str(p)) as pdf:
            plot_composition_drift_page(pdf, m)
        pypdf = pytest.importorskip("pypdf")
        # no page written on empty
        assert not os.path.exists(str(p)) or len(pypdf.PdfReader(str(p)).pages) == 0


def test_scorecard_feed():
    # 5 samples, one (ODD) with an atypical composition -> high mean_jsd -> flagged.
    rows = {f"S{i}": {"FSM": 60, "ISM": 40} for i in range(4)}
    rows["ODD"] = {"FSM": 10, "ISM": 90}
    dm = compute_composition_drift(_pivot(rows))
    metrics = assemble_scorecard_metrics(list(rows), drift_metrics=dm)
    assert all("composition_drift" in r for r in metrics)
    cfg = {"sample_scorecard": {"min_samples": 4, "z_warn": 2.0, "z_fail": 3.5,
                                "min_metrics_warn": 1, "min_metrics_fail": 1,
                                "metrics": {"composition_drift": "high"}}}
    sc = compute_sample_scorecard(metrics, cfg)
    assert sc["enabled"] is True
    assert sc["overall"]["ODD"] in ("warn", "fail")
