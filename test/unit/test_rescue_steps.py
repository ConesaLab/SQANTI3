"""
Unit tests for run_candidate_mapping in rescue_steps.py.
minimap2 is mocked; the tests check the targets FASTA deduplication and
when an existing SAM file is reused.
"""
import os
import sys
import shutil
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio import SeqIO

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.rescue_steps import run_candidate_mapping

TEST_SAM = Path(main_path) / "test" / "test_data" / "rescue" / "test_mapping.sam"

REF_FASTA = ">REF1\nACGTACGT\n>REF2\nACGTACGA\n>REF3\nACGTACGC\n>REF1\nACGTACGT\n"
# REF1 is also reported as a long read transcript (Bambu-like)
LR_FASTA = ">PB.1.1\nACGTAAAA\n>PB.2.1\nACGTCCCC\n>REF1\nACGTACGT\n>PB.9.1\nACGTGGGG\n"

TARGETS = ["REF1", "REF2", "REF3"]
CANDIDATES = ["PB.1.1", "PB.2.1"]


@pytest.fixture
def workdir():
    with tempfile.TemporaryDirectory() as tmpdir:
        (Path(tmpdir) / "logs" / "rescue").mkdir(parents=True)
        (Path(tmpdir) / "ref.fasta").write_text(REF_FASTA)
        (Path(tmpdir) / "corrected.fasta").write_text(LR_FASTA)
        yield tmpdir


def fake_minimap2(cmd, *args, **kwargs):
    """Write a valid SAM to the redirection target of the minimap2 command."""
    out_file = cmd.split(">")[-1].strip()
    shutil.copy(TEST_SAM, out_file)


def run(workdir, targets=TARGETS, candidates=CANDIDATES):
    return run_candidate_mapping(
        f"{workdir}/ref.fasta", targets, candidates,
        f"{workdir}/corrected.fasta", workdir, "isoform"
    )


@patch('src.rescue_steps.run_command', side_effect=fake_minimap2)
def test_targets_fasta_has_unique_ids(mock_run, workdir):
    run(workdir)

    ids = [r.id for r in SeqIO.parse(f"{workdir}/isoform_rescue_targets.fasta", "fasta")]
    assert sorted(ids) == ["REF1", "REF2", "REF3"]


@patch('src.rescue_steps.run_command', side_effect=fake_minimap2)
def test_first_run_maps_and_writes_fingerprint(mock_run, workdir):
    hits = run(workdir)

    sam = Path(workdir) / "isoform_mapped_rescue.sam"
    assert mock_run.call_count == 1
    assert sam.is_file()
    assert Path(f"{sam}.sha256").is_file()
    assert not Path(f"{sam}.tmp").exists()
    assert not hits.empty


@patch('src.rescue_steps.run_command', side_effect=fake_minimap2)
def test_same_inputs_reuse_sam(mock_run, workdir):
    run(workdir)
    run(workdir)

    assert mock_run.call_count == 1


@patch('src.rescue_steps.run_command', side_effect=fake_minimap2)
def test_sam_without_fingerprint_is_remapped(mock_run, workdir):
    """SAM files from older SQANTI3 versions have no fingerprint and may be broken."""
    shutil.copy(TEST_SAM, Path(workdir) / "isoform_mapped_rescue.sam")

    run(workdir)

    assert mock_run.call_count == 1
    assert Path(workdir, "isoform_mapped_rescue.sam.sha256").is_file()


@patch('src.rescue_steps.run_command', side_effect=fake_minimap2)
def test_changed_inputs_are_remapped(mock_run, workdir):
    run(workdir)
    run(workdir, targets=["REF1", "REF2"])
    run(workdir, candidates=["PB.1.1"])

    assert mock_run.call_count == 3


@patch('src.rescue_steps.run_command', side_effect=SystemExit(1))
def test_failed_mapping_leaves_no_sam(mock_run, workdir):
    with pytest.raises(SystemExit):
        run(workdir)

    assert not Path(workdir, "isoform_mapped_rescue.sam").exists()
    assert not Path(workdir, "isoform_mapped_rescue.sam.sha256").exists()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
