import os
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

import sqanti3_reads


class TestGetMethodRunSQANTI3:
    """Test suite for get_method_runSQANTI3."""

    def test_mods_mode_stops_when_qc_fails(self, tmp_path, monkeypatch):
        """align_uBAM has already replaced the SAM, so going on after a failed QC would join the new
        alignment to the classification of a previous run."""
        monkeypatch.setattr(sqanti3_reads, "sqantiqcPath", main_path, raising=False)
        monkeypatch.setattr(sqanti3_reads, "get_mods_info", lambda args, df, logger: {"s1": "s1.bam"})
        monkeypatch.setattr(sqanti3_reads, "align_uBAM", lambda *a, **k: str(tmp_path / "s1.fastq"))
        # Any way of launching QC reports a failure
        monkeypatch.setattr(sqanti3_reads.subprocess, "call", lambda *a, **k: 1)

        def failed(*a, **k):
            raise sqanti3_reads.subprocess.CalledProcessError(1, "sqanti3_qc.py")
        monkeypatch.setattr(sqanti3_reads.subprocess, "check_call", failed)
        monkeypatch.setattr(sqanti3_reads.subprocess, "run", failed)

        args = SimpleNamespace(mods_mode=True, sqanti_dirs=None, input_dir=str(tmp_path), OUTPUT=str(tmp_path),
                               refGTF="ref.gtf", refFasta="genome.fasta", min_ref_len=0, cpus=1,
                               sites="ATAC,GCAG,GTAG", aligner_choice="minimap2", chunks=1, jobs=1,
                               verbose=False)
        df = pd.DataFrame({"file_acc": ["s1"], "sampleID": ["s1"]})
        with pytest.raises((SystemExit, sqanti3_reads.subprocess.CalledProcessError, RuntimeError)):
            sqanti3_reads.get_method_runSQANTI3(args, df)
