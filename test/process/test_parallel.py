import pytest
import subprocess
import os
import tempfile
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))

class TestParallelPipeline:
    def test_parallel_run_basic(self):
        """Test the SQANTI3 pipeline using the --chunks parameter to ensure parallel processing works."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "python", os.path.join(main_path, "sqanti3_qc.py"),
                "--isoforms", os.path.join(main_path, "test", "test_data", "isoforms", "test_isoforms.fasta"),
                "--refGTF", os.path.join(main_path, "test", "test_data", "reference", "test_reference.gtf"),
                "--refFasta", os.path.join(main_path, "test", "test_data", "genome", "genome_test.fasta"),
                "--dir", tmpdir,
                "--output", "test_parallel",
                "--chunks", "2",
                "--report", "skip"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            assert result.returncode == 0, f"QC parallel run failed!\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            
            # Check for expected outputs in parallel mode
            assert os.path.exists(os.path.join(tmpdir, "test_parallel_classification.txt"))
            assert os.path.exists(os.path.join(tmpdir, "test_parallel_junctions.txt"))

