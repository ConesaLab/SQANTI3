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
    get_rescue_reference_targets,
    identify_rescue_candidates
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
def extended_classification_df(test_data_dir):
    """Extended classification DataFrame for comprehensive testing - loaded from file."""
    file_path = test_data_dir / "extended_classification.tsv"
    return pd.read_csv(file_path, sep="\t")


@pytest.fixture
def mini_gtf_file(test_data_dir):
    """Return path to mini reference GTF file."""
    return test_data_dir / "reference_mini.gtf"


class TestIdentifyRescueCandidates:
    """Test suite for identify_rescue_candidates function."""
    
    def test_nic_artifacts_included(self, extended_classification_df):
        """NIC artifacts should be included in rescue candidates."""
        result = identify_rescue_candidates(extended_classification_df)
        
        # PB.2.1 is NIC artifact
        assert 'PB.2.1' in result['isoform'].values
    
    def test_nnc_artifacts_included(self, extended_classification_df):
        """NNC artifacts should be included in rescue candidates."""
        result = identify_rescue_candidates(extended_classification_df)
        
        # PB.2.2 is NNC artifact
        assert 'PB.2.2' in result['isoform'].values
    
    def test_ism_with_lost_reference_included(self, extended_classification_df):
        """ISM with lost reference (no FSM) should be included."""
        result = identify_rescue_candidates(extended_classification_df)
        
        # PB.3.1 has reference REF5 which has no FSM
        assert 'PB.3.1' in result['isoform'].values
        # PB.8.1 has reference REF10 which has no FSM
        assert 'PB.8.1' in result['isoform'].values
    
    def test_ism_with_existing_fsm_included(self, extended_classification_df):
        """ISM with a reference that kept an FSM isoform should be included."""
        result = identify_rescue_candidates(extended_classification_df)
        
        # PB.3.2 has reference REF2 which has FSM isoform (PB.1.2)
        assert 'PB.3.2' in result['isoform'].values
    
    def test_fsm_isoforms_excluded(self, extended_classification_df):
        """FSM transcripts that passed the filter should never be included."""
        result = identify_rescue_candidates(extended_classification_df)
        
        assert 'PB.1.1' not in result['isoform'].values
        assert 'PB.1.2' not in result['isoform'].values
    
    def test_artifacts_of_represented_reference_included(self):
        """FSM and ISM artifacts of a reference with FSM isoforms should be included."""
        df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2', 'PB.1.3', 'PB.1.4'],
            'structural_category': ['full-splice_match', 'full-splice_match',
                                    'full-splice_match', 'incomplete-splice_match'],
            'exons': [3, 3, 3, 2],
            'filter_result': ['Isoform', 'Isoform', 'Artifact', 'Artifact'],
            'associated_transcript': ['REF1', 'REF1', 'REF1', 'REF1'],
            'associated_gene': ['GENE1', 'GENE1', 'GENE1', 'GENE1']
        })
        result = identify_rescue_candidates(df)
        
        assert set(result['isoform']) == {'PB.1.3', 'PB.1.4'}
    
    def test_artifacts_of_lost_reference_excluded(self):
        """FSM and ISM artifacts of a lost reference with FSM go to automatic rescue."""
        df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'structural_category': ['full-splice_match', 'incomplete-splice_match'],
            'exons': [3, 2],
            'filter_result': ['Artifact', 'Artifact'],
            'associated_transcript': ['REF1', 'REF1'],
            'associated_gene': ['GENE1', 'GENE1']
        })
        result = identify_rescue_candidates(df)
        
        assert result.empty
    
    def test_ism_of_reference_with_only_ism_isoforms_excluded(self):
        """ISM artifacts of a reference represented only by ISM isoforms should be excluded."""
        df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'structural_category': ['incomplete-splice_match', 'incomplete-splice_match'],
            'exons': [2, 2],
            'filter_result': ['Isoform', 'Artifact'],
            'associated_transcript': ['REF1', 'REF1'],
            'associated_gene': ['GENE1', 'GENE1']
        })
        result = identify_rescue_candidates(df)
        
        assert result.empty
    
    def test_monoexonic_fsm_of_represented_reference(self):
        """Mono-exonic FSM artifacts follow the same mono-exonic filter as other candidates."""
        df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'structural_category': ['full-splice_match', 'full-splice_match'],
            'exons': [1, 1],
            'filter_result': ['Isoform', 'Artifact'],
            'associated_transcript': ['REF1', 'REF1'],
            'associated_gene': ['GENE1', 'GENE1']
        })
        
        assert 'PB.1.2' in identify_rescue_candidates(df, monoexons='all')['isoform'].values
        assert identify_rescue_candidates(df, monoexons='fsm').empty
    
    def test_non_artifacts_excluded(self, extended_classification_df):
        """Non-artifact NIC/NNC should be excluded."""
        result = identify_rescue_candidates(extended_classification_df)
        
        # PB.5.1 is NIC but filter_result is 'Isoform'
        assert 'PB.5.1' not in result['isoform'].values
        # PB.6.1 is NNC but filter_result is 'Isoform'
        assert 'PB.6.1' not in result['isoform'].values
    
    def test_monoexon_filtering_all(self, extended_classification_df):
        """With monoexons='all', monoexon candidates should be included."""
        result = identify_rescue_candidates(extended_classification_df, monoexons='all')
        
        # PB.4.1 is ISM artifact with 1 exon
        assert 'PB.4.1' in result['isoform'].values
        # PB.7.1 is NIC artifact with 1 exon
        assert 'PB.7.1' in result['isoform'].values
    
    def test_monoexon_filtering_excluded(self, extended_classification_df):
        """With monoexons != 'all', monoexon candidates should be excluded."""
        result = identify_rescue_candidates(extended_classification_df, monoexons='fsm')
        
        # PB.4.1 is ISM artifact with 1 exon - should be excluded
        assert 'PB.4.1' not in result['isoform'].values
        # PB.7.1 is NIC artifact with 1 exon - should be excluded
        assert 'PB.7.1' not in result['isoform'].values
        
        # Multi-exon candidates should still be included
        assert 'PB.2.1' in result['isoform'].values  # NIC, 2 exons
        assert 'PB.3.1' in result['isoform'].values  # ISM, 2 exons
    
    def test_empty_dataframe(self):
        """Empty DataFrame should return empty result."""
        empty_df = pd.DataFrame({
            'isoform': [],
            'structural_category': [],
            'exons': [],
            'filter_result': [],
            'associated_transcript': [],
            'associated_gene': []
        })
        
        result = identify_rescue_candidates(empty_df)
        
        assert len(result) == 0
    
    def test_no_candidates(self):
        """DataFrame with no valid candidates should return empty result."""
        no_candidates_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'structural_category': ['full-splice_match', 'full-splice_match'],
            'exons': [2, 3],
            'filter_result': ['Isoform', 'Isoform'],
            'associated_transcript': ['REF1', 'REF2'],
            'associated_gene': ['GENE1', 'GENE2']
        })
        
        result = identify_rescue_candidates(no_candidates_df)
        
        assert len(result) == 0
    
    def test_only_nic_nnc_artifacts(self):
        """DataFrame with only NIC/NNC artifacts should return all of them."""
        nic_nnc_only = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'structural_category': ['novel_in_catalog', 'novel_not_in_catalog'],
            'exons': [2, 3],
            'filter_result': ['Artifact', 'Artifact'],
            'associated_transcript': ['REF1', 'REF2'],
            'associated_gene': ['GENE1', 'GENE2']
        })
        
        result = identify_rescue_candidates(nic_nnc_only)
        
        assert len(result) == 2
        assert 'PB.1.1' in result['isoform'].values
        assert 'PB.1.2' in result['isoform'].values
    
    def test_result_is_dataframe(self, extended_classification_df):
        """Result should be a pandas DataFrame."""
        result = identify_rescue_candidates(extended_classification_df)
        
        assert isinstance(result, pd.DataFrame)
    
    def test_result_columns_preserved(self, extended_classification_df):
        """Result should have same columns as input DataFrame."""
        result = identify_rescue_candidates(extended_classification_df)
        
        if len(result) > 0:
            assert list(result.columns) == list(extended_classification_df.columns)
    
    def test_candidate_count(self, extended_classification_df):
        """Verify expected number of candidates from extended test data."""
        result = identify_rescue_candidates(extended_classification_df, monoexons='all')
        
        # Expected candidates:
        # - PB.2.1 (NIC artifact, REF3)
        # - PB.2.2 (NNC artifact, REF4)
        # - PB.3.1 (ISM with lost reference REF5)
        # - PB.3.2 (ISM of REF2, which has FSM isoform PB.1.2)
        # - PB.4.1 (ISM with lost reference REF6, monoexon)
        # - PB.7.1 (NIC artifact REF9, monoexon)
        # - PB.8.1 (ISM with lost reference REF10)
        assert len(result) == 7
    
    def test_candidate_count_no_monoexons(self, extended_classification_df):
        """Verify expected number of candidates when monoexons are filtered."""
        result = identify_rescue_candidates(extended_classification_df, monoexons='exclude')
        
        # Expected candidates (excluding monoexons):
        # - PB.2.1 (NIC artifact, REF3, 2 exons)
        # - PB.2.2 (NNC artifact, REF4, 3 exons)
        # - PB.3.1 (ISM with lost reference REF5, 2 exons)
        # - PB.3.2 (ISM of REF2, which has FSM isoform PB.1.2, 2 exons)
        # - PB.8.1 (ISM with lost reference REF10, 4 exons)
        assert len(result) == 5


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
