"""
Unit tests for rescue_output.py functions.
Tests GTF and FASTA output generation functions.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os
import tempfile
from Bio import SeqIO

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.rescue_output import (
    write_rescue_gtf,
    write_fasta_file,
    read_fasta_file,
    write_rescue_fasta
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
def test_gtf_file(test_data_dir):
    """Return path to test GTF file."""
    return test_data_dir / "reference_mini.gtf"


@pytest.fixture
def test_isoforms_gtf():
    """Return path to test isoforms GTF file."""
    return Path(main_path) / "test" / "test_data" / "isoforms" / "test_isoforms.gtf"


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_sequences():
    """Create sample Bio.SeqRecord sequences for testing."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    seq1 = SeqRecord(Seq("ATGCATGCATGC"), id="TEST1", description="Test sequence 1")
    seq2 = SeqRecord(Seq("GCTAGCTAGCTA"), id="TEST2", description="Test sequence 2")
    seq3 = SeqRecord(Seq("TTTTAAAACCCC"), id="TEST3", description="Test sequence 3")
    
    return [seq1, seq2, seq3]

class TestReadFastaFile:
    """Test suite for read_fasta_file function."""
    
    def test_read_all_sequences(self, test_fasta_file):
        """Test reading all sequences without filter."""
        records = read_fasta_file(str(test_fasta_file))
        
        assert len(records) == 5
        record_ids = [rec.id for rec in records]
        assert 'REF1' in record_ids
        assert 'REF2' in record_ids
        assert 'REF3' in record_ids
        assert 'REF4' in record_ids
        assert 'REF5' in record_ids
    
    def test_read_with_single_filter(self, test_fasta_file):
        """Test reading with single ID filter."""
        filter_ids = ['REF1']
        
        records = read_fasta_file(str(test_fasta_file), filter_ids)
        
        assert len(records) == 1
        assert records[0].id == 'REF1'
    
    def test_read_with_multiple_filters(self, test_fasta_file):
        """Test reading with multiple ID filters."""
        filter_ids = ['REF1', 'REF3', 'REF5']
        
        records = read_fasta_file(str(test_fasta_file), filter_ids)
        
        assert len(records) == 3
        record_ids = [rec.id for rec in records]
        assert 'REF1' in record_ids
        assert 'REF3' in record_ids
        assert 'REF5' in record_ids
        assert 'REF2' not in record_ids
        assert 'REF4' not in record_ids
    
    def test_read_with_nonexistent_filter(self, test_fasta_file):
        """Test reading with filter IDs that don't exist."""
        filter_ids = ['REF999', 'REF888']
        
        records = read_fasta_file(str(test_fasta_file), filter_ids)
        
        assert len(records) == 0
    
    def test_read_with_empty_filter(self, test_fasta_file):
        """Test reading with empty filter list."""
        filter_ids = []
        
        records = read_fasta_file(str(test_fasta_file), filter_ids)
        
        assert len(records) == 0
    
    def test_read_with_mixed_filters(self, test_fasta_file):
        """Test reading with mix of existing and non-existing IDs."""
        filter_ids = ['REF1', 'REF999', 'REF3']
        
        records = read_fasta_file(str(test_fasta_file), filter_ids)
        
        assert len(records) == 2
        record_ids = [rec.id for rec in records]
        assert 'REF1' in record_ids
        assert 'REF3' in record_ids
        assert 'REF999' not in record_ids
    
    def test_sequences_preserved(self, test_fasta_file):
        """Test that sequence content is preserved."""
        filter_ids = ['REF2']
        
        records = read_fasta_file(str(test_fasta_file), filter_ids)
        
        assert len(records) == 1
        assert str(records[0].seq) == 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'


class TestWriteRescueFasta:
    """Test suite for write_rescue_fasta function."""
    
    def test_merge_transcripts_and_references(self, test_fasta_file, temp_output_dir):
        """Test merging good transcripts with rescued references."""
        prefix = str(Path(temp_output_dir) / "test_output")
        
        # Good transcripts from LR data
        good_transcripts = ['REF1', 'REF2']
        # Rescued references to add
        inclusion_list = ['REF4', 'REF5']
        
        write_rescue_fasta(
            str(test_fasta_file),  # input (also serves as ref for this test)
            str(test_fasta_file),  # reference
            good_transcripts,
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.fasta")
        assert output_file.exists()
        
        # Read and verify
        records = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(records) == 4  # 2 good + 2 rescued
        
        record_ids = [rec.id for rec in records]
        assert 'REF1' in record_ids
        assert 'REF2' in record_ids
        assert 'REF4' in record_ids
        assert 'REF5' in record_ids
    
    def test_no_good_transcripts(self, test_fasta_file, temp_output_dir):
        """Test with no good transcripts, only rescued."""
        prefix = str(Path(temp_output_dir) / "test_output")
        
        good_transcripts = []
        inclusion_list = ['REF1', 'REF2']
        
        write_rescue_fasta(
            str(test_fasta_file),
            str(test_fasta_file),
            good_transcripts,
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.fasta")
        assert output_file.exists()
        
        records = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(records) == 2
        record_ids = [rec.id for rec in records]
        assert 'REF1' in record_ids
        assert 'REF2' in record_ids
    
    def test_no_rescued_transcripts(self, test_fasta_file, temp_output_dir):
        """Test with no rescued transcripts, only good ones."""
        prefix = str(Path(temp_output_dir) / "test_output")
        
        good_transcripts = ['REF1', 'REF2', 'REF3']
        inclusion_list = []
        
        write_rescue_fasta(
            str(test_fasta_file),
            str(test_fasta_file),
            good_transcripts,
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.fasta")
        assert output_file.exists()
        
        records = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(records) == 3
    
    def test_overlapping_transcripts(self, test_fasta_file, temp_output_dir):
        """Test when same transcript is in both good and rescued lists."""
        prefix = str(Path(temp_output_dir) / "test_output")
        
        good_transcripts = ['REF1', 'REF2']
        inclusion_list = ['REF2', 'REF3']  # REF2 is in both
        
        write_rescue_fasta(
            str(test_fasta_file),
            str(test_fasta_file),
            good_transcripts,
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.fasta")
        assert output_file.exists()
        
        records = list(SeqIO.parse(str(output_file), "fasta"))
        # Should have duplicates if both lists contain same ID
        record_ids = [rec.id for rec in records]
        assert 'REF1' in record_ids
        assert 'REF3' in record_ids


class TestWriteRescueGtf:
    """Test suite for write_rescue_gtf function."""
    
    def test_basic_gtf_writing(self, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test basic GTF writing with inclusion list."""
        prefix = str(Path(temp_output_dir) / "test_output")
        inclusion_list = ['REF1', 'REF2']
        
        write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.gtf")
        assert output_file.exists()
        
        # Read output and verify
        with open(output_file) as f:
            lines = f.readlines()
        
        # Should have input lines plus rescued lines
        # test_isoforms.gtf has 240+ lines
        assert len(lines) > 219
        
        # Check that input transcripts are preserved (some from test_isoforms.gtf)
        transcript_lines = [l for l in lines if 'PB.124830.1' in l]
        assert len(transcript_lines) >= 2
    
    def test_gtf_filters_by_transcript_id(self, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test that only transcripts in inclusion list are added."""
        prefix = str(Path(temp_output_dir) / "test_output")
        inclusion_list = ['REF1']  # Only one transcript
        
        write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.gtf")
        
        with open(output_file) as f:
            content = f.read()
        
        # Should include the requested transcript
        assert 'REF1' in content
        # Should not include other transcripts from reference
        assert 'REF2' not in content
    
    def test_gtf_includes_transcript_and_exon_features(self, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test that both transcript and exon features are included."""
        prefix = str(Path(temp_output_dir) / "test_output")
        inclusion_list = ['REF1']
        
        write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.gtf")
        
        with open(output_file) as f:
            lines = [l for l in f if l.startswith('chr')]
        
        # Count transcript and exon features for the rescued transcript
        transcript_lines = [l for l in lines if 'REF1' in l and '\ttranscript\t' in l]
        exon_lines = [l for l in lines if 'REF1' in l and '\texon\t' in l]
        print(exon_lines)
        assert len(transcript_lines) == 1
        assert len(exon_lines) == 2
    
    def test_gtf_preserves_input_lines(self, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test that input GTF lines are preserved."""
        prefix = str(Path(temp_output_dir) / "test_output")
        inclusion_list = ['REF1']
        
        write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.gtf")
        
        with open(output_file) as f:
            lines = f.readlines()
        
        # Input lines should be preserved - check for specific transcript from test_isoforms.gtf
        pb_lines = [l for l in lines if 'PB.124830.1' in l]
        assert len(pb_lines) == 3
    
    def test_gtf_empty_inclusion_list(self, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test with empty inclusion list."""
        prefix = str(Path(temp_output_dir) / "test_output")
        inclusion_list = []
        
        write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        output_file = Path(f"{prefix}_rescued.gtf")
        assert output_file.exists()
        
        with open(output_file) as f:
            lines = [l for l in f if not l.startswith('#')]
        
        # Should only have input lines from test_isoforms.gtf (219 lines)
        assert len(lines) == 219
    
    def test_gtf_returns_output_path(self, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test that function returns the output file path."""
        prefix = str(Path(temp_output_dir) / "test_output")
        inclusion_list = ['ENST00000456328']
        
        result = write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        expected_path = f"{prefix}_rescued.gtf"
        assert result == expected_path
        assert Path(result).exists()


class TestIntegrationRescueOutput:
    """Integration tests combining multiple output functions."""
    
    def test_complete_rescue_output_workflow(self, test_fasta_file, test_isoforms_gtf, test_gtf_file, temp_output_dir):
        """Test complete workflow of generating rescue outputs."""
        prefix = str(Path(temp_output_dir) / "test_rescue")
        
        # Define rescue scenario
        good_transcripts = ['REF1', 'REF2']
        inclusion_list = ['REF1']
        
        # Generate GTF output
        gtf_output = write_rescue_gtf(
            str(test_isoforms_gtf),
            str(test_gtf_file),
            inclusion_list,
            prefix
        )
        
        # Generate FASTA output
        write_rescue_fasta(
            str(test_fasta_file),
            str(test_fasta_file),
            good_transcripts,
            inclusion_list,
            prefix
        )
        
        # Verify both outputs exist
        assert Path(gtf_output).exists()
        assert Path(f"{prefix}_rescued.fasta").exists()
        
        # Verify GTF has correct content
        with open(gtf_output) as f:
            gtf_content = f.read()
        # Check for transcripts from test_isoforms.gtf
        assert 'PB.124830.1' in gtf_content
        assert 'REF1' in gtf_content
        
        # Verify FASTA has correct content
        fasta_records = list(SeqIO.parse(f"{prefix}_rescued.fasta", "fasta"))
        fasta_ids = [rec.id for rec in fasta_records]
        assert 'REF1' in fasta_ids
        assert 'REF2' in fasta_ids


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
