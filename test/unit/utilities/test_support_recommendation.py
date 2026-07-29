"""Unit tests for Part C — the closing 'optional orthogonal support' note listing
support types that were not supplied. Renders only when something is missing, and
is deliberately non-judgemental (optional, not required).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import pytest
from matplotlib.backends.backend_pdf import PdfPages

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

from src.utilities.sqanti_reads_plots import (
    missing_support_types,
    plot_optional_support_note_page,
)
from src.utilities.sqanti_reads_report import _optional_support_footer

_ABSENT = {"samples": [], "per_sample": None}
_PRESENT = {"samples": ["A", "B"], "per_sample": None}


def test_all_absent_names_all_three():
    missing = missing_support_types(_ABSENT, _ABSENT, _ABSENT)
    labels = " ".join(l for l, _ in missing).lower()
    flags = " ".join(f for _, f in missing)
    assert len(missing) == 3
    assert "cage" in labels and "polya" in labels and "short-read" in labels
    assert "--CAGE_peak" in flags and "--polyA_peak" in flags and "--short_reads" in flags


def test_fully_supported_is_empty():
    assert missing_support_types(_PRESENT, _PRESENT, _PRESENT) == []


def test_partial_missing():
    # CAGE present, polyA + short-read absent -> two entries.
    missing = missing_support_types(_PRESENT, _ABSENT, _ABSENT)
    labels = " ".join(l for l, _ in missing).lower()
    assert len(missing) == 2 and "cage" not in labels


def test_pdf_page_only_when_missing(tmp_path):
    pypdf = pytest.importorskip("pypdf")
    # missing -> one page
    p = tmp_path / "note.pdf"
    with PdfPages(str(p)) as pdf:
        plot_optional_support_note_page(pdf, _ABSENT, _ABSENT, _ABSENT)
    reader = pypdf.PdfReader(str(p))
    assert len(reader.pages) == 1
    text = (reader.pages[0].extract_text() or "").lower()
    assert "optional" in text and "not" in text  # "...not required..."
    # fully supported -> no page
    p2 = tmp_path / "none.pdf"
    with PdfPages(str(p2)) as pdf:
        plot_optional_support_note_page(pdf, _PRESENT, _PRESENT, _PRESENT)
    assert not os.path.exists(str(p2)) or len(pypdf.PdfReader(str(p2)).pages) == 0


def test_html_footer_only_when_missing():
    html = _optional_support_footer(_ABSENT, _ABSENT, _ABSENT)
    assert "Optional orthogonal support" in html
    assert "optional" in html.lower() and "required" in html.lower()
    assert "--CAGE_peak" in html
    # fully supported -> empty string (no section)
    assert _optional_support_footer(_PRESENT, _PRESENT, _PRESENT) == ""
