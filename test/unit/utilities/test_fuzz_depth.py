"""Unit tests for F4 — splice-site imprecision vs junction depth. The chr22
fixture has ``total_coverage_unique`` 100% NaN, so F4 no-ops there; the real
behaviour is validated on a synthetic jxn_DF with coverage populated.
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
    _summarize_fuzz_by_depth,
    compute_fuzz_depth_metrics,
    plot_fuzz_depth_page,
)
from src.utilities.sqanti_reads_report import _fuzz_depth_figure


def _jxn(cov, dstart, dend, n, strand="+"):
    """n junctions, all with coverage `cov` and start/end offsets dstart/dend."""
    return pd.DataFrame({
        "chrom": ["chr1"] * n, "strand": [strand] * n,
        "genomic_start_coord": list(range(1000, 1000 + n)),
        "genomic_end_coord": list(range(5000, 5000 + n)),
        "diff_to_Ref_start_site": [dstart] * n,
        "diff_to_Ref_end_site": [dend] * n,
        "total_coverage_unique": [cov] * n,
    })


def test_imprecision_decreases_with_depth():
    # low-depth junctions imprecise (offset 3/2), high-depth precise (offset 0/0)
    jxn = pd.concat([_jxn(2, 3, 2, 10), _jxn(100, 0, 0, 10)], ignore_index=True)
    fd = _summarize_fuzz_by_depth(jxn, "S", "temp_factor", 0)
    m = compute_fuzz_depth_metrics(fd)
    prof = m["profile"].set_index("depth_bin")["perc_imprecise"]
    assert prof.loc["1-4"] == pytest.approx(100.0)   # all low-depth obs imprecise
    assert prof.loc["100+"] == pytest.approx(0.0)     # all high-depth obs precise


def test_all_nan_coverage_is_noop(tmp_path):
    jxn = _jxn(2, 3, 2, 5)
    jxn["total_coverage_unique"] = pd.array([pd.NA] * len(jxn), dtype="Int64")
    fd = _summarize_fuzz_by_depth(jxn, "S", "temp_factor", 0)
    assert fd.empty and list(fd.columns) == ["sampleID", "temp_factor", "depth_bin",
                                             "n_sites", "n_imprecise"]
    m = compute_fuzz_depth_metrics(fd)
    assert m["samples"] == []
    assert _fuzz_depth_figure(m) is None
    p = tmp_path / "fd.pdf"
    with PdfPages(str(p)) as pdf:
        plot_fuzz_depth_page(pdf, m)
    pypdf = pytest.importorskip("pypdf")
    assert not os.path.exists(str(p)) or len(pypdf.PdfReader(str(p)).pages) == 0


def test_missing_column_is_noop():
    jxn = _jxn(2, 3, 2, 5).drop(columns=["total_coverage_unique"])
    fd = _summarize_fuzz_by_depth(jxn, "S", "temp_factor", 0)
    assert fd.empty


def test_page_and_figure_build(tmp_path):
    jxn = pd.concat([_jxn(2, 3, 2, 10), _jxn(30, 1, 0, 10), _jxn(100, 0, 0, 10)],
                    ignore_index=True)
    m = compute_fuzz_depth_metrics(_summarize_fuzz_by_depth(jxn, "S", "temp_factor", 0))
    assert _fuzz_depth_figure(m) is not None
    p = tmp_path / "fd.pdf"
    with PdfPages(str(p)) as pdf:
        plot_fuzz_depth_page(pdf, m)
    pypdf = pytest.importorskip("pypdf")
    assert len(pypdf.PdfReader(str(p)).pages) == 1
