"""Shared fixtures/helpers for the SQANTI-reads integration tests.

``test_sqanti_reads.py`` (CSV / page-count / title / HTML-section regression)
drives the tool end-to-end on the vendored chr22 fixture in
``test/test_data/sqanti_reads`` via :func:`run_reads`.
"""
import os
import shutil
import subprocess
import sys

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "test", "test_data", "sqanti_reads")
SQANTI_DIRS = os.path.join(DATA_DIR, "sqanti_dirs")
REF_GTF = os.path.join(DATA_DIR, "reference", "gencode.v38.basic_chr22.gtf")
BASELINE_DIR = os.path.join(DATA_DIR, "baseline")

# Read-count cutoff that lets the full pipeline run on this small fixture
# (the fixture's genes top out at ~53 reads, so the default -ge 100 yields none).
GENE_EXPRESSION = "10"

# External tools used by the UJC-hashing step; the reads suite is skipped if any
# are missing so the wider test run stays portable.
REQUIRED_TOOLS = ["gffread", "gtftools", "bedtools"]
_missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
requires_external_tools = pytest.mark.skipif(
    bool(_missing),
    reason="missing external tools required by SQANTI-reads: " + ", ".join(_missing),
)


def run_reads(out_dir, design, factor=None, report="pdf", config=None, jobs=None):
    """Run ``sqanti3_reads.py`` in fast mode on the fixture; return CompletedProcess."""
    cmd = [
        sys.executable,
        os.path.join(REPO_ROOT, "sqanti3_reads.py"),
        "--design", os.path.join(DATA_DIR, design),
        "--sqanti_dirs", SQANTI_DIRS,
        "--refGTF", REF_GTF,
        "--refFasta", REF_GTF,  # fast mode does not read it; satisfies argparse today
        "--output", out_dir,
        "--prefix", "sqantiReads",
        "-ge", GENE_EXPRESSION,
        "--report", report,
    ]
    if factor:
        cmd += ["--factor", factor]
    if config:
        cmd += ["--config", config]
    if jobs:
        cmd += ["--jobs", str(jobs)]
    return subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
