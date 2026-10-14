"""Unit tests for F6 — tandem splice-site (NAGNAG) detector. The real fixture has
almost no tandem signal, so the arithmetic is validated on a synthetic offset
spectrum (a spike at ±3 bp). Covers compute, PDF page + HTML no-op.
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

from src.utilities.sqanti_reads_plots import compute_tandem_sites, plot_tandem_sites_page
from src.utilities.sqanti_reads_report import _tandem_sites_figure


def _metrics(spectrum_rows):
    """spectrum_rows: list of (sampleID, offset, count) -> jxn_offset_metrics-like."""
    sp = pd.DataFrame([{"sampleID": s, "site": "donor", "offset": o, "count": c}
                       for s, o, c in spectrum_rows])
    return {"samples": sorted(sp["sampleID"].unique().tolist()), "spectrum": sp}


def test_spike_at_3bp_high_fraction():
    # 30 obs at +3, 10 spread elsewhere -> 30/40 = 75% tandem at 3bp.
    rows = [("A", 3, 30), ("A", 1, 5), ("A", 2, 3), ("A", 7, 2), ("A", 0, 100)]
    m = compute_tandem_sites(_metrics(rows))
    ps = m["per_sample"].set_index("sampleID")
    assert ps.loc["A", "perc_tandem_3bp"] == pytest.approx(75.0)
    # +/-4 and +/-6 add nothing here, so perc_tandem == perc_tandem_3bp
    assert ps.loc["A", "perc_tandem"] == pytest.approx(75.0)


def test_flat_spectrum_low_fraction():
    rows = [("A", o, 5) for o in range(-10, 11) if o != 0]  # uniform, no tandem excess
    m = compute_tandem_sites(_metrics(rows))
    ps = m["per_sample"].set_index("sampleID")
    # 2 of 20 nonzero offsets are ±3 -> 10%
    assert ps.loc["A", "perc_tandem_3bp"] == pytest.approx(10.0)


def test_empty_safe(tmp_path):
    for m_in in (None, {"samples": [], "spectrum": pd.DataFrame()},
                 _metrics([("A", 0, 50)])):  # only exact (offset 0) -> no imprecise obs
        m = compute_tandem_sites(m_in)
        assert m["samples"] == []
        assert _tandem_sites_figure(m) is None
        p = tmp_path / "t.pdf"
        with PdfPages(str(p)) as pdf:
            plot_tandem_sites_page(pdf, m)
        pypdf = pytest.importorskip("pypdf")
        assert not os.path.exists(str(p)) or len(pypdf.PdfReader(str(p)).pages) == 0


def test_page_and_figure_build(tmp_path):
    m = compute_tandem_sites(_metrics([("A", 3, 20), ("A", 4, 5), ("A", 1, 10),
                                       ("B", -3, 8), ("B", 2, 12)]))
    assert _tandem_sites_figure(m) is not None
    p = tmp_path / "t.pdf"
    with PdfPages(str(p)) as pdf:
        plot_tandem_sites_page(pdf, m)
    pypdf = pytest.importorskip("pypdf")
    assert len(pypdf.PdfReader(str(p)).pages) == 1
