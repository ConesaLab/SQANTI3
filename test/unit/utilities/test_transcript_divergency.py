"""Unit tests for transcript divergency (TD; Monzó et al. 2025, Genome Research).

Per annotated gene, TD = (NIC+NNC)/(FSM+NIC+NNC), over genes with FSM+NIC+NNC >=
min_reads; per-sample mean; and pairwise between-group per-gene TD differences
when faceted (reads pooled within a group per gene).
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
    compute_transcript_divergency,
    plot_transcript_divergency_page,
)
from src.utilities.sqanti_reads_report import _td_difference_figure


def _gc(rows, exp_factor="temp_factor"):
    """rows: list of (gene, sample, factorval, FSM, NIC, NNC)."""
    df = pd.DataFrame(rows, columns=["associated_gene", "sampleID", exp_factor,
                                     "full-splice_match", "novel_in_catalog",
                                     "novel_not_in_catalog"])
    df["flag_annotated_gene"] = 1
    return df


def test_td_formula_and_filter():
    # G1: FSM8/NIC1/NNC1 -> denom10, TD=0.2; G2: FSM20 -> TD=0; G3: denom5 -> dropped.
    df = _gc([("G1", "S", 0, 8, 1, 1), ("G2", "S", 0, 20, 0, 0), ("G3", "S", 0, 3, 1, 1)])
    m = compute_transcript_divergency(df, "temp_factor", min_reads=10)
    pg = m["per_gene"].set_index("associated_gene")["TD"]
    assert set(pg.index) == {"G1", "G2"}                      # G3 filtered
    assert pg["G1"] == pytest.approx(0.2)
    assert pg["G2"] == pytest.approx(0.0)
    assert m["per_sample"].set_index("sampleID")["td_mean"]["S"] == pytest.approx(0.1)
    assert m["group_pairs"] == []                             # not faceted


def test_only_annotated_genes():
    df = _gc([("G1", "S", 0, 8, 1, 1), ("NOVELG", "S", 0, 0, 5, 5)])
    df.loc[df["associated_gene"] == "NOVELG", "flag_annotated_gene"] = 0
    m = compute_transcript_divergency(df, "temp_factor", min_reads=10)
    assert set(m["per_gene"]["associated_gene"]) == {"G1"}


def test_faceted_pairwise_two_groups():
    df = _gc([("GX", "A", "g1", 8, 2, 0),     # group g1 pooled TD = 2/10 = 0.2
              ("GX", "B", "g2", 18, 1, 1)],   # group g2 pooled TD = 2/20 = 0.1
             exp_factor="grp")
    m = compute_transcript_divergency(df, "grp", min_reads=10)
    assert len(m["group_pairs"]) == 1
    p = m["group_pairs"][0]
    assert (p["a"], p["b"]) == ("g1", "g2")
    assert p["diffs"][0] == pytest.approx(0.1)   # TD[g1]-TD[g2]
    assert p["n_a_higher"] == 1 and p["n_b_higher"] == 0 and p["n_genes"] == 1


def test_faceted_three_groups_all_pairwise():
    rows = []
    for g, fsm in [("g1", 10), ("g2", 10), ("g3", 10)]:
        rows.append(("GX", f"s_{g}", g, fsm, 2, 0))
    m = compute_transcript_divergency(_gc(rows, exp_factor="grp"), "grp", min_reads=5)
    assert len(m["group_pairs"]) == 3            # 3 choose 2
    assert {(p["a"], p["b"]) for p in m["group_pairs"]} == {("g1", "g2"), ("g1", "g3"), ("g2", "g3")}


def test_empty_is_noop(tmp_path):
    # all genes below threshold -> empty
    df = _gc([("G1", "S", 0, 2, 1, 0)])
    m = compute_transcript_divergency(df, "temp_factor", min_reads=10)
    assert m["samples"] == [] and m["per_gene"].empty
    assert compute_transcript_divergency(None, "temp_factor")["samples"] == []
    # PDF page draws nothing on empty
    p = tmp_path / "td.pdf"
    with PdfPages(str(p)) as pdf:
        plot_transcript_divergency_page(pdf, m)
    pypdf = pytest.importorskip("pypdf")
    assert not os.path.exists(str(p)) or len(pypdf.PdfReader(str(p)).pages) == 0
    # HTML difference figure is None on empty
    assert _td_difference_figure({"a": "x", "b": "y", "diffs": np.array([]), "n_genes": 0}) is None


def test_pages_and_figure_build(tmp_path):
    df = _gc([("GX", "A", "g1", 8, 2, 0), ("GX", "B", "g2", 18, 1, 1),
              ("GY", "A", "g1", 30, 0, 0), ("GY", "B", "g2", 12, 3, 5)], exp_factor="grp")
    m = compute_transcript_divergency(df, "grp", min_reads=10)
    # PDF: per-sample bar + per-gene box + one difference histogram = 3 pages
    p = tmp_path / "td.pdf"
    with PdfPages(str(p)) as pdf:
        plot_transcript_divergency_page(pdf, m)
    pypdf = pytest.importorskip("pypdf")
    assert len(pypdf.PdfReader(str(p)).pages) == 3
    # HTML difference figure builds
    assert _td_difference_figure(m["group_pairs"][0]) is not None
