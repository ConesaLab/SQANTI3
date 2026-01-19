"""
Unit tests for candidate_mapping_helpers.py functions.
Tests helper functions for FASTA manipulation and SAM file processing.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os
import tempfile
from unittest.mock import patch, MagicMock
from Bio import SeqIO

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.candidate_mapping_helpers import (
    filter_transcriptome,
    process_sam_file,
    prepare_fasta_transcriptome
)


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(main_path) / "test" / "test_data" / "rescue"


@pytest.fixture
def test_fasta_file(test_data_dir):
    """Return path to test FASTA file."""
    return test_data_dir / "test_transcriptome.fasta"


@pytest.fixture
def test_sam_file(test_data_dir):
    """Return path to test SAM file."""
    return test_data_dir / "test_mapping.sam"


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create logs subdirectory
        logs_dir = Path(tmpdir) / "logs"
        logs_dir.mkdir(exist_ok=True)
        yield tmpdir


class TestFilterTranscriptome:
    """Test suite for filter_transcriptome function."""
    
    def test_filter_single_target(self, test_fasta_file):
        """Filter for a single target ID."""
        target_ids = ['REF1']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 1
        assert result[0].id == 'REF1'
        assert len(result[0].seq) > 0
    
    def test_filter_multiple_targets(self, test_fasta_file):
        """Filter for multiple target IDs."""
        target_ids = ['REF1', 'REF3', 'REF5']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 3
        result_ids = [rec.id for rec in result]
        assert 'REF1' in result_ids
        assert 'REF3' in result_ids
        assert 'REF5' in result_ids
    
    def test_filter_no_matches(self, test_fasta_file):
        """Filter with target IDs that don't exist in FASTA."""
        target_ids = ['REF999', 'REF888']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 0
    
    def test_filter_empty_target_list(self, test_fasta_file):
        """Filter with empty target list."""
        target_ids = []
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 0
    
    def test_filter_all_targets(self, test_fasta_file):
        """Filter requesting all sequences in the file."""
        target_ids = ['REF1', 'REF2', 'REF3', 'REF4', 'REF5']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 5
    
    def test_filter_preserves_sequences(self, test_fasta_file):
        """Verify that sequences are preserved correctly."""
        target_ids = ['REF2']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 1
        assert result[0].id == 'REF2'
        # The sequence should match what's in the file
        assert str(result[0].seq) == 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'
    
    def test_filter_with_duplicates_in_targets(self, test_fasta_file):
        """Filter with duplicate IDs in target list."""
        target_ids = ['REF1', 'REF1', 'REF2']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        # Should only return unique sequences even if targets have duplicates
        result_ids = [rec.id for rec in result]
        assert 'REF1' in result_ids
        assert 'REF2' in result_ids
    
    def test_filter_partial_match(self, test_fasta_file):
        """Filter with mix of existing and non-existing IDs."""
        target_ids = ['REF1', 'REF999', 'REF3']
        
        result = filter_transcriptome(str(test_fasta_file), target_ids)
        
        assert len(result) == 2
        result_ids = [rec.id for rec in result]
        assert 'REF1' in result_ids
        assert 'REF3' in result_ids
        assert 'REF999' not in result_ids


class TestProcessSamFile:
    """Test suite for process_sam_file function."""
    
    def test_basic_sam_processing(self, test_sam_file):
        """Test basic SAM file processing."""
        result_df = process_sam_file(str(test_sam_file))
        
        assert isinstance(result_df, pd.DataFrame)
        assert result_df.shape == (7, 4)  
        assert list(result_df.columns) == ['rescue_candidate', 'mapping_hit', 'alignment_type', 'alignment_score']
    
    def test_sam_unmapped_read(self, test_sam_file):
        """Test that unmapped reads are included."""
        result_df = process_sam_file(str(test_sam_file))
        
        # PB.1.2 is unmapped (flag 4)
        unmapped = result_df[result_df['rescue_candidate'] == 'PB.1.2']
        assert len(unmapped) == 1
        assert unmapped.iloc[0]['alignment_type'] == 4
        # Unmapped reads have None as reference_name in pysam
        assert pd.isna(unmapped.iloc[0]['mapping_hit']) or unmapped.iloc[0]['mapping_hit'] == '*'
    
    def test_sam_mapped_reads(self, test_sam_file):
        """Test that mapped reads have correct information."""
        result_df = process_sam_file(str(test_sam_file))
        
        # PB.1.1 should be mapped to REF1 with AS score 200
        mapped = result_df[result_df['rescue_candidate'] == 'PB.1.1']
        assert len(mapped) == 1
        assert mapped.iloc[0]['mapping_hit'] == 'REF1'
        assert mapped.iloc[0]['alignment_type'] == 0
        assert mapped.iloc[0]['alignment_score'] == 200
    
    def test_sam_multiple_alignments(self, test_sam_file):
        """Test reads with multiple alignments."""
        result_df = process_sam_file(str(test_sam_file))
        
        # PB.3.1 has two alignments
        multi_aligned = result_df[result_df['rescue_candidate'] == 'PB.3.1']
        assert len(multi_aligned) == 2
        
        hits = set(multi_aligned['mapping_hit'])
        assert 'REF1' in hits
        assert 'REF2' in hits
    
    def test_sam_reverse_complement(self, test_sam_file):
        """Test reads mapped in reverse complement."""
        result_df = process_sam_file(str(test_sam_file))
        
        # PB.4.1 is reverse complemented (flag 16)
        reverse = result_df[result_df['rescue_candidate'] == 'PB.4.1']
        assert len(reverse) == 1
        assert reverse.iloc[0]['alignment_type'] == 16
        assert reverse.iloc[0]['mapping_hit'] == 'REF3'
    
    def test_sam_missing_alignment_score(self, test_sam_file):
        """Test handling of reads without AS tag."""
        result_df = process_sam_file(str(test_sam_file))
        
        # PB.4.1 doesn't have AS tag, should default to 0
        no_as_tag = result_df[result_df['rescue_candidate'] == 'PB.4.1']
        assert len(no_as_tag) == 1
        assert no_as_tag.iloc[0]['alignment_score'] == 0
    
    def test_sam_all_reads_processed(self, test_sam_file):
        """Test that all reads in SAM file are processed."""
        result_df = process_sam_file(str(test_sam_file))
        
        # Should have 7 rows (matching the 7 alignment lines in the SAM file)
        assert len(result_df) == 7
    
    def test_sam_dataframe_structure(self, test_sam_file):
        """Test DataFrame has correct structure."""
        result_df = process_sam_file(str(test_sam_file))
        
        # Check column names
        expected_columns = ['rescue_candidate', 'mapping_hit', 'alignment_type', 'alignment_score']
        assert list(result_df.columns) == expected_columns
        
        # Check data types are reasonable
        assert result_df['rescue_candidate'].dtype == object
        assert result_df['alignment_type'].dtype in [int, 'int64', 'int32']
        assert result_df['alignment_score'].dtype in [int, 'int64', 'int32']


class TestPrepareFastaTranscriptome:
    """Test suite for prepare_fasta_transcriptome function."""
    
    @patch('src.utilities.rescue.candidate_mapping_helpers.run_command')
    def test_prepare_fasta_success(self, mock_run_command, temp_output_dir, test_data_dir):
        """Test successful FASTA transcriptome preparation."""
        ref_gtf = str(test_data_dir / "reference_mini.gtf")
        ref_fasta = str(test_data_dir / "genome.fasta")
        
        # Create a fake output file for the test
        expected_output = Path(temp_output_dir) / "reference_mini.fasta"
        expected_output.write_text(">transcript1\nATGC\n")
        
        # Mock run_command to create the file
        def side_effect(*args, **kwargs):
            expected_output.touch()
        
        mock_run_command.side_effect = side_effect
        
        result = prepare_fasta_transcriptome(ref_gtf, ref_fasta, temp_output_dir)
        
        assert result == str(expected_output)
        assert mock_run_command.called
        
        # Verify the command was called with correct arguments
        call_args = mock_run_command.call_args[0][0]
        assert "gffread" in call_args
        assert ref_gtf in call_args
        assert ref_fasta in call_args
    
    @patch('src.utilities.rescue.candidate_mapping_helpers.run_command')
    def test_prepare_fasta_command_format(self, mock_run_command, temp_output_dir, test_data_dir):
        """Test that gffread command is formatted correctly."""
        ref_gtf = str(test_data_dir / "reference_mini.gtf")
        ref_fasta = str(test_data_dir / "genome.fasta")
        
        # Create expected output file
        expected_output = Path(temp_output_dir) / "reference_mini.fasta"
        expected_output.write_text(">transcript1\nATGC\n")
        
        prepare_fasta_transcriptome(ref_gtf, ref_fasta, temp_output_dir)
        
        # Get the command that was called
        call_args = mock_run_command.call_args[0][0]
        
        # Verify command structure
        assert call_args.startswith("gffread")
        assert "-w" in call_args
        assert "-g" in call_args
        assert str(expected_output) in call_args
    
    @patch('src.utilities.rescue.candidate_mapping_helpers.run_command')
    @patch('src.utilities.rescue.candidate_mapping_helpers.sys.exit')
    def test_prepare_fasta_file_not_created(self, mock_exit, mock_run_command, temp_output_dir, test_data_dir):
        """Test error handling when output file is not created."""
        ref_gtf = str(test_data_dir / "reference_mini.gtf")
        ref_fasta = str(test_data_dir / "genome.fasta")
        
        # Don't create the output file to simulate failure
        prepare_fasta_transcriptome(ref_gtf, ref_fasta, temp_output_dir)
        
        # Should call sys.exit(1) when file doesn't exist
        mock_exit.assert_called_once_with(1)
    
    @patch('src.utilities.rescue.candidate_mapping_helpers.run_command')
    def test_prepare_fasta_output_naming(self, mock_run_command, temp_output_dir):
        """Test that output file is named correctly based on GTF name."""
        ref_gtf = "/path/to/my_annotation.gtf"
        ref_fasta = "/path/to/genome.fasta"
        
        expected_output = Path(temp_output_dir) / "my_annotation.fasta"
        expected_output.write_text(">transcript1\nATGC\n")
        
        result = prepare_fasta_transcriptome(ref_gtf, ref_fasta, temp_output_dir)
        
        assert result == str(expected_output)
        assert "my_annotation.fasta" in result


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
