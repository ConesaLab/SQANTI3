"""
Integration tests for SQANTI3 TUSCO functionality
"""
import pytest
import subprocess
import os
import tempfile
import sys

# Add the parent directory to the path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))


class TestTUSCOIntegration:
    """Integration tests for TUSCO benchmarking functionality"""

    @pytest.mark.integration
    def test_tusco_human_chr22(self):
        """
        Test TUSCO with human gene set on chromosome 22 test data.
        This test runs the full SQANTI3 QC pipeline with TUSCO benchmarking enabled.
        """
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as tmpdir:
            # Construct the command
            cmd = [
                "python", "sqanti3_qc.py",
                "--isoforms", "test/test_data/test_isoforms.gtf",
                "--refGTF", "data/reference/gencode.v38.basic_chr22.gtf",
                "--refFasta", "data/reference/GRCh38.p13_chr22.fasta",
                "--tusco", "human",
                "-o", "TUSCO_human_test",
                "-d", tmpdir,
                ""  # Skip ORF prediction to speed up the test
            ]

            # Run the command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            # Check that the command succeeded
            assert result.returncode == 0, f"SQANTI3 QC with TUSCO failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"

            # Check that key output files were created
            expected_files = [
                "TUSCO_human_test_classification.txt",
                "TUSCO_human_test_junctions.txt",
                "TUSCO_human_test_TUSCO_results.tsv",
                "TUSCO_human_test_TUSCO_report.html"
            ]

            for filename in expected_files:
                filepath = os.path.join(tmpdir, filename)
                assert os.path.exists(filepath), f"Expected output file not found: {filename}"

                # Check that files are not empty
                assert os.path.getsize(filepath) > 0, f"Output file is empty: {filename}"

            # Optionally, check TUSCO results file has expected content
            tusco_results_path = os.path.join(tmpdir, "TUSCO_human_test_TUSCO_results.tsv")
            with open(tusco_results_path, 'r') as f:
                lines = f.readlines()
                # Should have header + at least some data
                assert len(lines) > 1, "TUSCO results file appears to have no data"

                # Check header contains expected columns
                header = lines[0].strip()
                expected_columns = ["transcript_id", "associated_gene", "TUSCO_category"]
                for col in expected_columns:
                    assert col in header, f"Missing expected column in TUSCO results: {col}"

    @pytest.mark.integration
    @pytest.mark.skipif(not os.path.exists("data/mouse_test_data.gtf"),
                        reason="Mouse test data not available")
    def test_tusco_mouse(self):
        """
        Test TUSCO with mouse gene set (if test data is available).
        This test is skipped if mouse test data is not present.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Similar test structure for mouse data
            # This is a placeholder for when mouse test data is available
            pass

    @pytest.mark.integration
    def test_tusco_invalid_species(self):
        """
        Test that TUSCO properly handles invalid species specification.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "python", "sqanti3_qc.py",
                "--isoforms", "data/UHR_chr22.gtf",
                "--refGTF", "data/reference/gencode.v38.basic_chr22.gtf",
                "--refFasta", "data/reference/GRCh38.p13_chr22.fasta",
                "--tusco", "invalid_species",
                "-o", "TUSCO_invalid_test",
                "-d", tmpdir,
                "--skipORF"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30
            )

            # Should fail with non-zero exit code
            assert result.returncode != 0, "SQANTI3 should fail with invalid TUSCO species"

            # Check error message mentions invalid species
            error_output = result.stderr.lower()
            assert "invalid" in error_output or "species" in error_output or "tusco" in error_output, \
                "Error message should mention the invalid TUSCO species"


def test_tusco_import():
    """Test that TUSCO-related modules can be imported"""
    try:
        from src.qc_pipeline import run
        # If we get here, import was successful
        assert True
    except ImportError as e:
        # Check if it's a missing dependency that would be in conda environment
        if "Bio" in str(e) or "bx" in str(e) or "BCBio" in str(e):
            pytest.skip(f"Skipping test - requires conda environment: {e}")
        else:
            pytest.fail(f"Failed to import TUSCO-related modules: {e}")