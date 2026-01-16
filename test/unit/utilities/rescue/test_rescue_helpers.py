"""
Unit tests for rescue_helpers.py functions.
Tests helper functions for target selection and GTF parsing.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.rescue_helpers import (
    get_rescue_gene_targets,
    parse_rescue_gtf,
    get_rescue_reference_targets
)


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(main_path) / "test" / "test_data" / "rescue"


@pytest.fixture
def sample_classification_df():
    """Sample classification DataFrame for testing."""
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.1.2', 'PB.2.1', 'PB.3.1'],
        'structural_category': ['full-splice_match', 'full-splice_match', 
                               'novel_in_catalog', 'incomplete-splice_match'],
        'exons': [2, 2, 3, 2],
        'filter_result': ['Isoform', 'Artifact', 'Artifact', 'Artifact'],
        'associated_transcript': ['REF1', 'REF1', 'REF2', 'REF3'],
        'associated_gene': ['GENE1', 'GENE1', 'GENE2', 'GENE3']
    })


@pytest.fixture
def mini_gtf_file(test_data_dir):
    """Return path to mini reference GTF file."""
    return test_data_dir / "reference_mini.gtf"

class TestGetRescueGeneTargets:
    """Test suite for get_rescue_gene_targets function."""
    
    def test_single_candidate_one_gene(self, sample_classification_df):
        """Single candidate should return its associated gene."""
        rescue_candidates = ['PB.2.1']
        
        result = get_rescue_gene_targets(sample_classification_df, rescue_candidates)
        
        assert len(result) == 1
        assert 'GENE2' in result
    
    def test_multiple_candidates_same_gene(self, sample_classification_df):
        """Multiple candidates from same gene should return gene only once."""
        rescue_candidates = ['PB.1.1', 'PB.1.2']
        
        result = get_rescue_gene_targets(sample_classification_df, rescue_candidates)
        
        assert len(result) == 1
        assert 'GENE1' in result
    
    def test_multiple_candidates_different_genes(self, sample_classification_df):
        """Multiple candidates from different genes should return all genes."""
        rescue_candidates = ['PB.2.1', 'PB.3.1']
        
        result = get_rescue_gene_targets(sample_classification_df, rescue_candidates)
        
        assert len(result) == 2
        assert 'GENE2' in result
        assert 'GENE3' in result
    
    def test_candidate_not_in_dataframe(self, sample_classification_df):
        """Candidate not in DataFrame should return empty."""
        rescue_candidates = ['PB.999.1']
        
        result = get_rescue_gene_targets(sample_classification_df, rescue_candidates)
        
        assert len(result) == 0
    
    def test_empty_candidate_list(self, sample_classification_df):
        """Empty candidate list should return empty result."""
        rescue_candidates = []
        
        result = get_rescue_gene_targets(sample_classification_df, rescue_candidates)
        
        assert len(result) == 0


class TestParseRescueGtf:
    """Test suite for parse_rescue_gtf function."""
    
    def test_parse_valid_gtf(self, mini_gtf_file):
        """Test parsing a valid GTF file."""
        result = parse_rescue_gtf(str(mini_gtf_file))
        
        assert isinstance(result, dict)
        assert len(result) > 0
        
        # Check that genes are keys
        assert 'GENE1' in result
        assert 'GENE2' in result
        assert 'GENE3' in result
    
    def test_gene_with_multiple_transcripts(self, mini_gtf_file):
        """Test that gene with multiple transcripts lists all of them."""
        result = parse_rescue_gtf(str(mini_gtf_file))
        
        # GENE1 should have 2 transcripts (REF1 and REF1b)
        assert 'GENE1' in result
        assert isinstance(result['GENE1'], list)
        assert len(result['GENE1']) == 2
        assert 'REF1' in result['GENE1']
        assert 'REF1b' in result['GENE1']
    
    def test_gene_with_single_transcript(self, mini_gtf_file):
        """Test that gene with single transcript returns list with one item."""
        result = parse_rescue_gtf(str(mini_gtf_file))
        
        # GENE2 should have 1 transcript (REF2)
        assert 'GENE2' in result
        assert isinstance(result['GENE2'], list)
        assert len(result['GENE2']) == 1
        assert 'REF2' in result['GENE2']
    
    def test_all_transcripts_captured(self, mini_gtf_file):
        """Test that all unique transcripts are captured."""
        result = parse_rescue_gtf(str(mini_gtf_file))
        
        # Flatten all transcripts
        all_transcripts = [t for transcripts in result.values() for t in transcripts]
        
        # Should have 4 unique transcripts: REF1, REF1b, REF2, REF3
        assert len(all_transcripts) == 4
        assert 'REF1' in all_transcripts
        assert 'REF1b' in all_transcripts
        assert 'REF2' in all_transcripts
        assert 'REF3' in all_transcripts


class TestGetRescueReferenceTargets:
    """Test suite for get_rescue_reference_targets function."""
    
    def test_target_genes_in_gtf(self, mini_gtf_file):
        """Test getting transcripts for genes that exist in GTF."""
        target_genes = ['GENE1', 'GENE2']
        
        result = get_rescue_reference_targets(str(mini_gtf_file), target_genes)
        
        assert isinstance(result, pd.Series)
        assert len(result) > 0
        
        # Should contain transcripts from both genes
        assert 'REF1' in result.values
        assert 'REF1b' in result.values
        assert 'REF2' in result.values
        
        # Should NOT contain transcripts from GENE3
        assert 'REF3' not in result.values
    
    def test_target_genes_not_in_gtf(self, mini_gtf_file):
        """Test getting transcripts for genes that don't exist in GTF."""
        target_genes = ['GENE999']
        
        result = get_rescue_reference_targets(str(mini_gtf_file), target_genes)
        
        assert isinstance(result, pd.Series)
        assert len(result) == 0
    
    def test_single_target_gene(self, mini_gtf_file):
        """Test getting transcripts for a single target gene."""
        target_genes = ['GENE1']
        
        result = get_rescue_reference_targets(str(mini_gtf_file), target_genes)
        
        assert len(result) == 2  # GENE1 has 2 transcripts
        assert 'REF1' in result.values
        assert 'REF1b' in result.values
    
    def test_empty_target_genes(self, mini_gtf_file):
        """Test with empty target genes list."""
        target_genes = []
        
        result = get_rescue_reference_targets(str(mini_gtf_file), target_genes)
        
        assert len(result) == 0
    
    def test_mixed_present_and_absent_genes(self, mini_gtf_file):
        """Test with mix of genes that do and don't exist in GTF."""
        target_genes = ['GENE1', 'GENE999', 'GENE3']
        
        result = get_rescue_reference_targets(str(mini_gtf_file), target_genes)
        
        # Should only get transcripts from GENE1 and GENE3
        assert len(result) == 3  # GENE1 (2 transcripts) + GENE3 (1 transcript)
        assert 'REF1' in result.values
        assert 'REF1b' in result.values
        assert 'REF3' in result.values
        assert 'REF2' not in result.values


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
