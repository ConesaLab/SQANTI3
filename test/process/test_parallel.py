import pytest
import subprocess
import os
import tempfile
import sys

from Bio import SeqIO

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

    def test_parallel_run_fastq(self):
        """Test that --chunks splits a FASTQ input and gives the same classification as its FASTA."""
        fasta = os.path.join(main_path, "test", "test_data", "isoforms", "test_isoforms.fasta")
        with tempfile.TemporaryDirectory() as tmpdir:
            # Quality 29 is written as '>', so a FASTQ read as FASTA would give wrong records
            fastq = os.path.join(tmpdir, "test_isoforms.fastq")
            records = []
            for r in SeqIO.parse(fasta, "fasta"):
                r.letter_annotations["phred_quality"] = [29] * len(r)
                records.append(r)
            SeqIO.write(records, fastq, "fastq")

            classifications = []
            for isoforms in (fasta, fastq):
                outdir = os.path.join(tmpdir, os.path.basename(isoforms).split('.')[-1])
                cmd = [
                    "python", os.path.join(main_path, "sqanti3_qc.py"),
                    "--isoforms", isoforms,
                    "--refGTF", os.path.join(main_path, "test", "test_data", "reference", "test_reference.gtf"),
                    "--refFasta", os.path.join(main_path, "test", "test_data", "genome", "genome_test.fasta"),
                    "--dir", outdir,
                    "--output", "test_parallel",
                    "--chunks", "2",
                    "--report", "skip"
                ]

                result = subprocess.run(cmd, capture_output=True, text=True)

                assert result.returncode == 0, f"QC parallel run failed on {isoforms}!\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
                with open(os.path.join(outdir, "test_parallel_classification.txt")) as h:
                    classifications.append(h.read())

            assert classifications[0] == classifications[1]

    def test_parallel_run_aligns_once(self):
        """Test that --chunks aligns a FASTA input once instead of once per chunk."""
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
            with open(os.path.join(tmpdir, "logs", "sqanti3_qc.log")) as h:
                assert h.read().count("Aligning reads with Minimap2") == 1
