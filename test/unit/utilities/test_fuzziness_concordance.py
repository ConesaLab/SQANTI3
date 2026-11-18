"""Unit tests for F7 — replicate concordance of splice-site precision, and the
new _summarize_site_offsets summarizer (whose ref_key must match calc_jxn_cv).
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
    _summarize_site_offsets,
    compute_fuzziness_concordance,
    plot_fuzziness_concordance_page,
)
from src.utilities.sqanti_reads_report import _fuzziness_concordance_figure


def _site_df(rows):
    """rows: (sampleID, group, ref_key, site_type, median_offset)."""
    return pd.DataFrame([{"sampleID": s, "group": g, "ref_key": rk,
                          "site_type": st, "median_offset": mo, "n": 5}
                         for s, g, rk, st, mo in rows])


# ---- _summarize_site_offsets: ref_key must match calc_jxn_cv's reconstruction ----

def test_ref_key_matches_calc_jxn_cv():
    # + strand: start->donor ref = start+diff; end->acceptor ref = end+diff.
    jxn = pd.DataFrame({
        "isoform": ["i1"], "chrom": ["chr1"], "strand": ["+"],
        "genomic_start_coord": [1000], "genomic_end_coord": [2000],
        "diff_to_Ref_start_site": [2], "diff_to_Ref_end_site": [-1],
        "junction_number": ["J1"], "junction_category": ["known"], "canonical": ["canonical"],
    })
    out = _summarize_site_offsets(jxn, "S", "temp_factor", 0).set_index("site_type")
    # calc_jxn_cv: ref_junction_start = genomic_start_coord + diff_to_Ref_start_site
    assert out.loc["donor", "ref_key"] == "chr1:+:1002"     # 1000 + 2
    assert out.loc["acceptor", "ref_key"] == "chr1:+:1999"  # 2000 + (-1)
    assert out.loc["donor", "median_offset"] == 2
    assert out.loc["acceptor", "median_offset"] == -1


def test_summarizer_empty_safe():
    cols = ["sampleID", "temp_factor", "ref_key", "site_type", "median_offset", "n"]
    assert list(_summarize_site_offsets(None, "S", "temp_factor", 0).columns) == cols
    assert _summarize_site_offsets(pd.DataFrame(), "S", "temp_factor", 0).empty


# ---- compute_fuzziness_concordance ----

def test_identical_replicates_high_agreement():
    sites = [("K1", "donor", 0), ("K2", "acceptor", 1), ("K3", "donor", 2)]
    rows = [("A", "g", rk, st, mo) for rk, st, mo in sites] + \
           [("B", "g", rk, st, mo) for rk, st, mo in sites]
    m = compute_fuzziness_concordance(_site_df(rows), {"A": "g", "B": "g"})
    assert m["enabled"] is True
    ps = m["per_sample"].set_index("sampleID")
    assert ps.loc["A", "site_precision_agreement"] == pytest.approx(1.0, abs=1e-6)


def test_disagreeing_replicate_low_agreement():
    # A and B agree; C shifts every shared site by >tol (3bp).
    base = [("K1", "donor", 0), ("K2", "acceptor", 0), ("K3", "donor", 0)]
    rows = ([("A", "g", rk, st, mo) for rk, st, mo in base]
            + [("B", "g", rk, st, mo) for rk, st, mo in base]
            + [("C", "g", rk, st, 8) for rk, st, _ in base])
    ps = compute_fuzziness_concordance(_site_df(rows), {"A": "g", "B": "g", "C": "g"}
                                       )["per_sample"].set_index("sampleID")
    assert ps.loc["A", "site_precision_agreement"] > ps.loc["C", "site_precision_agreement"]
    assert ps.loc["C", "site_precision_agreement"] == pytest.approx(0.0, abs=1e-6)


def test_disabled_and_empty(tmp_path):
    df = _site_df([("A", "g", "K1", "donor", 0), ("B", "g", "K1", "donor", 0)])
    assert compute_fuzziness_concordance(df, None)["enabled"] is False
    assert compute_fuzziness_concordance(pd.DataFrame(), {"A": "g", "B": "g"})["enabled"] is False
    m_single = compute_fuzziness_concordance(df, {"A": "g1", "B": "g2"})
    assert m_single["enabled"] is False and "replicate" in m_single["reason"]
    assert _fuzziness_concordance_figure(m_single) is None
    p = tmp_path / "f.pdf"
    with PdfPages(str(p)) as pdf:
        plot_fuzziness_concordance_page(pdf, m_single)
    pypdf = pytest.importorskip("pypdf")
    assert not os.path.exists(str(p)) or len(pypdf.PdfReader(str(p)).pages) == 0
