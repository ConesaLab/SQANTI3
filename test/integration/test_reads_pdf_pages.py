"""Page-image regression harness for the SQANTI-reads PDF report.

This is the safety net for the plot-renderer unification: the CSV/page-count/
title tests in ``test_sqanti_reads.py`` do not notice a plot that renders with the
wrong bars, colours or layout. Here every page of the fixture report is
rasterized with PyMuPDF and compared, pixel-wise, against a committed reference
image, so any visual drift across the ~30 plots fails a test.

References live in ``test/test_data/sqanti_reads/baseline/report_pages/`` as one
PNG per page plus a ``page_hashes.json`` manifest (sha256 per page). Regenerate
them deliberately after an intended visual change::

    SQANTI_READS_UPDATE_REFS=1 pytest test/integration/test_reads_pdf_pages.py

Assertion layers (see the roadmap §2 note on cross-version pixel jitter):
  * page count always checked;
  * each page must be non-blank and match the reference dimensions;
  * pixels compared with a small mean-absolute-difference tolerance (tolerates
    freetype/matplotlib anti-aliasing jitter, catches real regressions);
  * the exact sha256 manifest match is additionally enforced only when
    ``SQANTI_READS_STRICT_PIXELS=1`` (for a pinned CI image).
"""
import hashlib
import json
import os

import numpy as np
import pytest

from .conftest import BASELINE_DIR, requires_external_tools, run_reads

fitz = pytest.importorskip("fitz")  # PyMuPDF; self-contained, no system poppler

pytestmark = requires_external_tools

# Rasterization + comparison knobs.
RENDER_DPI = 100
# Mean absolute pixel difference (0-255 scale) allowed before a page is a
# regression. Small enough to catch a wrong bar/colour, loose enough to absorb
# anti-aliasing jitter across matplotlib/freetype builds.
MAD_TOLERANCE = float(os.environ.get("SQANTI_READS_MAD_TOLERANCE", "2.0"))
UPDATE_REFS = os.environ.get("SQANTI_READS_UPDATE_REFS") == "1"
STRICT_PIXELS = os.environ.get("SQANTI_READS_STRICT_PIXELS") == "1"

REF_ROOT = os.path.join(BASELINE_DIR, "report_pages")

# (design file, --factor arg, reference subdir).
DESIGNS = [
    ("design.csv", None, "design"),
    ("design_with_factor.csv", "group", "design_with_factor"),
]


def _page_array(page):
    """Render a PDF page to an (H, W, 3) uint8 RGB array at RENDER_DPI."""
    pix = page.get_pixmap(dpi=RENDER_DPI, alpha=False)
    arr = np.frombuffer(pix.samples, dtype=np.uint8)
    return arr.reshape(pix.height, pix.width, pix.n)


def _sha256(arr):
    return hashlib.sha256(arr.tobytes()).hexdigest()


def _render_report(tmp_path, design, factor):
    out_dir = str(tmp_path / "out")
    result = run_reads(out_dir, design, factor=factor)
    assert result.returncode == 0, f"non-zero exit\nSTDERR:\n{result.stderr}"
    pdf = os.path.join(out_dir, "sqantiReads_report.pdf")
    assert os.path.isfile(pdf), "report PDF was not produced"
    return pdf


def _write_references(subdir, pages):
    ref_dir = os.path.join(REF_ROOT, subdir)
    os.makedirs(ref_dir, exist_ok=True)
    # Drop any stale page images from a previous (higher) page count.
    for fname in os.listdir(ref_dir):
        if fname.startswith("page_") and fname.endswith(".png"):
            os.remove(os.path.join(ref_dir, fname))
    manifest = {}
    for i, arr in enumerate(pages):
        name = f"page_{i:02d}.png"
        _save_png(os.path.join(ref_dir, name), arr)
        manifest[name] = {
            "sha256": _sha256(arr),
            "height": int(arr.shape[0]),
            "width": int(arr.shape[1]),
        }
    with open(os.path.join(ref_dir, "page_hashes.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)


def _save_png(path, arr):
    height, width = arr.shape[0], arr.shape[1]
    pix = fitz.Pixmap(fitz.csRGB, width, height, arr.tobytes(), False)
    pix.save(path)


def _load_png(path):
    pix = fitz.Pixmap(path)
    arr = np.frombuffer(pix.samples, dtype=np.uint8)
    return arr.reshape(pix.height, pix.width, pix.n)


@pytest.mark.parametrize("design,factor,subdir", DESIGNS, ids=[d[2] for d in DESIGNS])
def test_report_pages_match_reference(tmp_path, design, factor, subdir):
    pdf = _render_report(tmp_path, design, factor)
    with fitz.open(pdf) as doc:
        pages = [_page_array(doc[i]) for i in range(doc.page_count)]

    ref_dir = os.path.join(REF_ROOT, subdir)
    manifest_path = os.path.join(ref_dir, "page_hashes.json")

    if UPDATE_REFS:
        _write_references(subdir, pages)
        pytest.skip(f"regenerated {len(pages)} reference pages for {subdir}")

    if not os.path.isfile(manifest_path):
        pytest.skip(
            f"no reference images for {subdir}; generate with "
            "SQANTI_READS_UPDATE_REFS=1"
        )

    with open(manifest_path) as fh:
        manifest = json.load(fh)

    assert len(pages) == len(manifest), (
        f"{subdir}: page count changed ({len(pages)} vs {len(manifest)} "
        "reference pages) — regenerate references if intentional"
    )

    failures = []
    for i, arr in enumerate(pages):
        name = f"page_{i:02d}.png"
        spec = manifest[name]
        # Structural checks (always).
        assert arr.shape[0] == spec["height"] and arr.shape[1] == spec["width"], (
            f"{subdir}/{name}: dimensions {arr.shape[:2]} != "
            f"({spec['height']}, {spec['width']})"
        )
        assert not np.all(arr == 255), f"{subdir}/{name}: page is blank"

        if _sha256(arr) == spec["sha256"]:
            continue  # exact match, nothing more to check
        if STRICT_PIXELS:
            failures.append(f"{name}: sha256 mismatch (strict mode)")
            continue
        # Tolerance path: compare against the committed reference PNG.
        ref = _load_png(os.path.join(ref_dir, name))
        mad = float(np.mean(np.abs(arr.astype(np.int16) - ref.astype(np.int16))))
        if mad > MAD_TOLERANCE:
            failures.append(f"{name}: mean abs diff {mad:.3f} > {MAD_TOLERANCE}")

    assert not failures, (
        f"{subdir}: {len(failures)} page(s) drifted from reference — inspect the "
        f"committed PNGs, and regenerate with SQANTI_READS_UPDATE_REFS=1 if the "
        f"change is intentional:\n  " + "\n  ".join(failures)
    )
