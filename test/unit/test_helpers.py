
import logging
from typing import Dict
import pytest
from Bio import SeqIO
from pathlib import Path
import sys, os


main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.helpers import (
    process_gtf_line, rename_isoform_seqids, get_corr_filenames, get_isoform_hits_name, 
    get_class_junc_filenames, get_omitted_name
)


class TestRenameIsoformSeqids:
    """Test suite for rename_isoform_seqids function."""
    
    @pytest.fixture
    def mock_fasta_file(self):
        """Fixture to return the path to the mock FASTA file."""
        return Path(main_path, "test", "test_data", "isoforms", "isoform_mock.fasta")

    @pytest.fixture
    def mock_fastq_file(self):
        """Fixture to return the path to the mock FASTQ file."""
        return Path(main_path, "test", "test_data", "isoforms", "isoform_mock.fastq")

    @pytest.fixture
    def mock_gz_file(self):
        """Fixture to return the path to the mock gzipped FASTA file."""
        return Path(main_path, "test", "test_data", "isoforms", "isoform_mock.fasta.gz")

    def test_rename_isoform_seqids_fasta(self, mock_fasta_file):
        """Test renaming IDs in a standard FASTA file."""
        # Arrange
        output_fasta = str(mock_fasta_file.with_name("isoform_mock.renamed.fasta"))
        expected_ids = ["PB.1.1", "PB.2.1", "PBfusion.3.1"]

        # Act
        result_file = rename_isoform_seqids(str(mock_fasta_file))
        
        # Assert
        assert result_file == output_fasta, "Output file name does not match expected."
        assert os.path.exists(output_fasta), "Output file was not created."

        # Verify the content of the output FASTA file
        with open(output_fasta, "r") as output:
            renamed_records = list(SeqIO.parse(output, "fasta"))
            assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

            for record, expected_id in zip(renamed_records, expected_ids):
                assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

        # Cleanup: Remove the output file after test
        os.remove(output_fasta)

    def test_rename_isoform_seqids_fastq(self, mock_fastq_file):
        """Test renaming IDs in a FASTQ file (output is FASTA)."""
        # Arrange
        output_fastq = str(mock_fastq_file.with_name("isoform_mock.renamed.fasta"))
        expected_ids = ["PB.1.1", "PB.2.1", "PBfusion.3.1"]

        # Act
        result_file = rename_isoform_seqids(str(mock_fastq_file))

        # Assert
        assert result_file == output_fastq, "Output file name does not match expected."
        assert os.path.exists(output_fastq), "Output file was not created."

        # Verify the content of the output FASTA file
        with open(output_fastq, "r") as output:
            renamed_records = list(SeqIO.parse(output, "fasta"))
            assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

            for record, expected_id in zip(renamed_records, expected_ids):
                assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

        # Cleanup: Remove the output file after test
        os.remove(output_fastq)

    def test_rename_isoform_seqids_gz(self, mock_gz_file):
        """Test renaming IDs in a gzipped FASTA file."""
        # Arrange
        output_fastq = str(mock_gz_file.with_name("isoform_mock.renamed.fasta"))
        expected_ids = ["PB.1.1", "PB.2.1", "PBfusion.3.1"]

        # Act
        result_file = rename_isoform_seqids(str(mock_gz_file))

        # Assert
        assert result_file == output_fastq, "Output file name does not match expected."
        assert os.path.exists(output_fastq), "Output file was not created."

        # Verify the content of the output FASTA file
        with open(output_fastq, "r") as output:
            renamed_records = list(SeqIO.parse(output, "fasta"))
            assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

            for record, expected_id in zip(renamed_records, expected_ids):
                assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

        # Cleanup: Remove the output file after test
        os.remove(output_fastq)

    def test_rename_isoform_seqids_arbitrary_format(self, mock_fasta_file):
        """Test that any ID format is accepted and properly cleaned."""
        # Arrange
        arbitrary_fasta = str(mock_fasta_file.with_name("arbitrary_format.fasta"))
        with open(arbitrary_fasta, "w") as f:
            f.write(">ENST00000456328|additional_info some description\nACGTACGT\n")
            f.write(">transcript_123 gene=ABC\nGCTAGCTA\n")
            f.write(">simple_id\nATATATAT\n")

        output_fasta = str(mock_fasta_file.with_name("arbitrary_format.renamed.fasta"))
        expected_ids = ["ENST00000456328", "transcript_123", "simple_id"]

        # Act
        result_file = rename_isoform_seqids(arbitrary_fasta)

        # Assert
        assert result_file == output_fasta, "Output file was not created."
        
        with open(output_fasta, "r") as output:
            renamed_records = list(SeqIO.parse(output, "fasta"))
            assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

            for record, expected_id in zip(renamed_records, expected_ids):
                assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

        # Cleanup
        os.remove(arbitrary_fasta)
        os.remove(output_fasta)


class TestFilenameGetters:
    """Test suite for filename getter functions."""

    @pytest.fixture
    def sample_paths(self):
        return {
            'outdir': '/tmp',
            'prefix': 'test'
        }

    def test_get_corr_filenames(self, sample_paths):
        """Test get_corr_filenames returns correct paths."""
        expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
        corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF = get_corr_filenames(**sample_paths)
        
        assert corrGTF == expected_prefix + "_corrected.gtf"
        assert corrSAM == expected_prefix + "_corrected.sam"
        assert corrFASTA == expected_prefix + "_corrected.fasta"
        assert corrORF == expected_prefix + "_corrected.faa"
        assert corrCDS_GTF_GFF == expected_prefix + "_corrected.cds.gff3"

    def test_get_isoform_hits_name(self, sample_paths):
        """Test get_isoform_hits_name returns correct path."""
        expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
        isoform_hits_name = get_isoform_hits_name(**sample_paths)
        
        assert isoform_hits_name == expected_prefix + "_isoform_hits.txt"

    def test_get_class_junc_filenames(self, sample_paths):
        """Test get_class_junc_filenames returns correct paths."""
        expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
        outputClassPath, outputJuncPath = get_class_junc_filenames(**sample_paths)
        
        assert outputClassPath == expected_prefix + "_classification.txt"
        assert outputJuncPath == expected_prefix + "_junctions.txt"

    def test_get_omitted_name(self, sample_paths):
        """Test get_omitted_name returns correct path."""
        expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
        omitted_name = get_omitted_name(**sample_paths)
        
        assert omitted_name == expected_prefix + "_omitted_due_to_min_ref_len.txt"

    @pytest.mark.parametrize("outdir,prefix", [
        ("/var/log", "app"),
        ("/home/user/documents", "report"),
        (".", "test"),
        ("/tmp/my folder", "test file"),
    ])
    def test_all_functions_parametrized(self, outdir, prefix):
        """Test all filename getters with various path combinations."""
        expected_prefix = os.path.abspath(os.path.join(outdir, prefix))
        
        # Test get_corr_filenames
        corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF = get_corr_filenames(outdir, prefix)
        assert corrGTF == expected_prefix + "_corrected.gtf"
        assert corrSAM == expected_prefix + "_corrected.sam"
        assert corrFASTA == expected_prefix + "_corrected.fasta"
        assert corrORF == expected_prefix + "_corrected.faa"
        assert corrCDS_GTF_GFF == expected_prefix + "_corrected.cds.gff3"
        
        # Test get_isoform_hits_name
        isoform_hits_name = get_isoform_hits_name(outdir, prefix)
        assert isoform_hits_name == expected_prefix + "_isoform_hits.txt"
        
        # Test get_class_junc_filenames
        outputClassPath, outputJuncPath = get_class_junc_filenames(outdir, prefix)
        assert outputClassPath == expected_prefix + "_classification.txt"
        assert outputJuncPath == expected_prefix + "_junctions.txt"
        
        # Test get_omitted_name
        omitted_name = get_omitted_name(outdir, prefix)
        assert omitted_name == expected_prefix + "_omitted_due_to_min_ref_len.txt"

    def test_empty_prefix(self):
        """Test filename getters with empty prefix."""
        outdir = "/tmp"
        prefix = ""
        expected_prefix = os.path.abspath(os.path.join(outdir, prefix))
        
        # Test all functions with empty prefix
        corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF = get_corr_filenames(outdir, prefix)
        isoform_hits_name = get_isoform_hits_name(outdir, prefix)
        outputClassPath, outputJuncPath = get_class_junc_filenames(outdir, prefix)
        omitted_name = get_omitted_name(outdir, prefix)
        
        assert all(path.startswith(expected_prefix) for path in [corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF, isoform_hits_name, outputClassPath, outputJuncPath, omitted_name])


class TestProcessGtfLine:
    """Test suite for process_gtf_line function."""

    @pytest.fixture
    def genome_dict(self) -> Dict[str, str]:
        return {"chr1": "ATCG", "chr2": "GCTA"}

    @pytest.fixture
    def mock_corrGTF_out(self):
        return main_path + "/test/test_data/corrected.gtf"

    @pytest.fixture
    def mock_discard_gtf(self):
        return main_path + "/test/test_data/discard.gtf"

    @pytest.fixture
    def tester_logger(self):
        """Create a logger for testing."""
        logger = logging.getLogger("tester_logger")
        logger.setLevel(logging.INFO)
        return logger

    def test_comment_line(self, genome_dict, mock_corrGTF_out, mock_discard_gtf, tester_logger):
        """Test that comment lines are ignored."""
        try:
            os.remove(mock_corrGTF_out)
            os.remove(mock_discard_gtf)
        except FileNotFoundError:
             # Files may not exist from previous runs; it's safe to ignore if they are absent.
            pass
        line = "# This is a comment\n"
        result = process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf, logger=tester_logger)
        assert result is None, "The result was not empty"
        assert not os.path.exists(mock_corrGTF_out), "The corrected GTF file was created"
        assert not os.path.exists(mock_discard_gtf), "The discarded GTF file was created"

    def test_malformed_line(self, genome_dict, mock_corrGTF_out, mock_discard_gtf, caplog, tester_logger):
        """Test that malformed lines are skipped with warning."""
        line = "chr1\tgene\n"
        process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf, logger=tester_logger)
        assert "Skipping malformed GTF line" in caplog.text
        assert not os.path.exists(mock_corrGTF_out)
        assert not os.path.exists(mock_discard_gtf)

    def test_chromosome_not_in_genome(self, genome_dict, mock_corrGTF_out, mock_discard_gtf, caplog, tester_logger):
        """Test that lines with unknown chromosomes raise ValueError."""
        caplog.set_level(logging.ERROR)
        line = "chr3\tEnsembl\texon\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
        with pytest.raises(ValueError):
            process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf, tester_logger)
        
        assert "GTF chromosome chr3 not found in genome reference file." in caplog.text

    def test_valid_transcript_line(self, genome_dict, mock_corrGTF_out, mock_discard_gtf):
        """Test that valid transcript lines are written to corrected GTF."""
        line = "chr1\tEnsembl\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
        
        with open(mock_corrGTF_out, "w") as f:
            process_gtf_line(line, genome_dict, f, mock_discard_gtf)
        f.close()
        with open(mock_corrGTF_out, "r") as f:
            assert f.read() == line
        assert not os.path.exists(mock_discard_gtf)
        os.remove(mock_corrGTF_out)

    def test_valid_exon_line(self, genome_dict, mock_corrGTF_out, mock_discard_gtf):
        """Test that valid exon lines are written to corrected GTF."""
        line = "chr2\tEnsembl\texon\t1\t1000\t.\t-\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
        with open(mock_corrGTF_out, "w") as f:
            process_gtf_line(line, genome_dict, f, mock_discard_gtf)
        f.close()
        with open(mock_corrGTF_out, "r") as f:
            assert f.read() == line
        assert not os.path.exists(mock_discard_gtf)
        os.remove(mock_corrGTF_out)

    def test_unknown_strand(self, genome_dict, mock_corrGTF_out, mock_discard_gtf, caplog, tester_logger):
        """Test that lines with unknown strand are written to discard GTF."""
        line = "chr1\tEnsembl\texon\t1\t1000\t.\t.\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
        with open(mock_discard_gtf, "w") as f:
            process_gtf_line(line, genome_dict, mock_corrGTF_out, f, logger=tester_logger)
        f.close()
        assert "Discarding unknown strand transcript" in caplog.text
        assert not os.path.exists(mock_corrGTF_out)
        with open(mock_discard_gtf, "r") as f:
            assert f.read() == line
        os.remove(mock_discard_gtf)

    def test_non_transcript_exon_line(self, genome_dict, mock_corrGTF_out, mock_discard_gtf):
        """Test that non-transcript/exon lines are ignored."""
        line = "chr1\tEnsembl\tgene\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\";\n"
        process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)
        assert not os.path.exists(mock_corrGTF_out)
        assert not os.path.exists(mock_discard_gtf)

    def test_adding_lines_to_file(self, genome_dict, mock_corrGTF_out, mock_discard_gtf):
        """Test that multiple valid lines are correctly written."""
        line1 = "chr1\tEnsembl\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
        line2 = "chr2\tEnsembl\texon\t1\t1000\t.\t-\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
        with open(mock_corrGTF_out, "w") as f:
            process_gtf_line(line1, genome_dict, f, mock_discard_gtf)
            process_gtf_line(line2, genome_dict, f, mock_discard_gtf)
        f.close()
        with open(mock_corrGTF_out, "r") as f:
            assert f.read() == line1 + line2
        assert not os.path.exists(mock_discard_gtf)
        os.remove(mock_corrGTF_out)
