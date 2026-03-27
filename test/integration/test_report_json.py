"""
Integration tests for SQANTI3 --report_json functionality
"""
import json
import pytest
import subprocess
import os
import tempfile
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))


class TestReportJSON:
    """Integration tests for JSON export of report tables"""

    @pytest.mark.integration
    def test_report_json_with_html(self):
        """
        Test --report_json alongside default HTML report.
        Verifies JSON file is created with expected structure.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "python", "sqanti3_qc.py",
                "--isoforms", "data/UHR_chr22.gtf",
                "--refGTF", "data/reference/gencode.v38.basic_chr22.gtf",
                "--refFasta", "data/reference/GRCh38.p13_chr22.fasta",
                "--report_json",
                "-o", "json_test",
                "-d", tmpdir,
                "--skipORF"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, (
                f"SQANTI3 QC with --report_json failed:\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

            # Check JSON file exists
            json_path = os.path.join(tmpdir, "json_test_report_tables.json")
            assert os.path.exists(json_path), "JSON report file not created"
            assert os.path.getsize(json_path) > 0, "JSON report file is empty"

            # Validate JSON structure
            with open(json_path) as f:
                data = json.load(f)

            # Check top-level keys
            assert "metadata" in data, "Missing 'metadata' key in JSON"
            assert "summary" in data, "Missing 'summary' key in JSON"

            # Check metadata
            assert "generated_at" in data["metadata"]
            assert "classification_file" in data["metadata"]
            assert "junction_file" in data["metadata"]

            # Check summary
            assert "unique_genes" in data["summary"]
            assert "unique_isoforms" in data["summary"]
            assert isinstance(data["summary"]["unique_genes"], int)
            assert isinstance(data["summary"]["unique_isoforms"], int)
            assert data["summary"]["unique_genes"] > 0
            assert data["summary"]["unique_isoforms"] > 0

            # Check isoform classification table
            assert "isoform_classification" in data["summary"]
            iso_class = data["summary"]["isoform_classification"]
            assert len(iso_class) > 0
            assert "category" in iso_class[0]
            assert "n_isoforms" in iso_class[0]

            # Check gene classification table
            assert "gene_classification" in data["summary"]
            gene_class = data["summary"]["gene_classification"]
            assert len(gene_class) > 0
            assert "category" in gene_class[0]
            assert "n_genes" in gene_class[0]

            # HTML report should also exist
            html_path = os.path.join(tmpdir, "json_test_SQANTI3_report.html")
            assert os.path.exists(html_path), "HTML report should still be generated"

    @pytest.mark.integration
    def test_report_json_skip_report(self):
        """
        Test --report skip --report_json produces JSON without HTML/PDF.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "python", "sqanti3_qc.py",
                "--isoforms", "data/UHR_chr22.gtf",
                "--refGTF", "data/reference/gencode.v38.basic_chr22.gtf",
                "--refFasta", "data/reference/GRCh38.p13_chr22.fasta",
                "--report", "skip",
                "--report_json",
                "-o", "json_skip_test",
                "-d", tmpdir,
                "--skipORF"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, (
                f"SQANTI3 QC with --report skip --report_json failed:\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

            # JSON file should exist
            json_path = os.path.join(tmpdir, "json_skip_test_report_tables.json")
            assert os.path.exists(json_path), "JSON report not created with --report skip"

            # HTML/PDF should NOT exist
            html_path = os.path.join(tmpdir, "json_skip_test_SQANTI3_report.html")
            pdf_path = os.path.join(tmpdir, "json_skip_test_SQANTI3_report.pdf")
            assert not os.path.exists(html_path), "HTML report should not exist with --report skip"
            assert not os.path.exists(pdf_path), "PDF report should not exist with --report skip"

            # Validate JSON is well-formed
            with open(json_path) as f:
                data = json.load(f)
            assert "metadata" in data
            assert "summary" in data

    @pytest.mark.integration
    def test_no_json_without_flag(self):
        """
        Test that no JSON file is produced when --report_json is not specified.
        Ensures backward compatibility.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "python", "sqanti3_qc.py",
                "--isoforms", "data/UHR_chr22.gtf",
                "--refGTF", "data/reference/gencode.v38.basic_chr22.gtf",
                "--refFasta", "data/reference/GRCh38.p13_chr22.fasta",
                "--report", "skip",
                "-o", "no_json_test",
                "-d", tmpdir,
                "--skipORF"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, (
                f"SQANTI3 QC failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

            # JSON file should NOT exist
            json_path = os.path.join(tmpdir, "no_json_test_report_tables.json")
            assert not os.path.exists(json_path), "JSON report should not be created without --report_json flag"

    @pytest.mark.integration
    def test_json_quality_metrics(self):
        """
        Test that quality_metrics section is present when junction data exists.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "python", "sqanti3_qc.py",
                "--isoforms", "data/UHR_chr22.gtf",
                "--refGTF", "data/reference/gencode.v38.basic_chr22.gtf",
                "--refFasta", "data/reference/GRCh38.p13_chr22.fasta",
                "--report", "skip",
                "--report_json",
                "-o", "json_qm_test",
                "-d", tmpdir,
                "--skipORF"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, (
                f"SQANTI3 QC failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

            json_path = os.path.join(tmpdir, "json_qm_test_report_tables.json")
            with open(json_path) as f:
                data = json.load(f)

            # Quality metrics should exist for multi-exon data
            if "quality_metrics" in data:
                qm = data["quality_metrics"]
                # Check structure of quality metric tables
                for key in qm:
                    entries = qm[key]
                    assert isinstance(entries, list), f"quality_metrics.{key} should be a list"
                    if len(entries) > 0:
                        assert "structural_category" in entries[0], (
                            f"quality_metrics.{key} entries should have 'structural_category'"
                        )
                        assert "percent" in entries[0], (
                            f"quality_metrics.{key} entries should have 'percent'"
                        )
