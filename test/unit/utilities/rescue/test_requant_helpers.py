"""
Unit tests for requantification helper functions.
Tests individual functions in src.utilities.rescue.requant_helpers module.
"""
import pytest
import pandas as pd
import numpy as np
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.requant_helpers import (
    get_unrescued_artifacts,
    prepare_count_matrices,
    calculate_distribution_fractions,
    distribute_integer_counts,
    export_counts,
    calculate_tpm
)
from src.utilities.rescue.sq_requant import build_artifact_table


class TestGetUnrescuedArtifacts:
    """Test suite for get_unrescued_artifacts function."""
    
    def test_basic_unrescued_artifacts(self):
        """Test identifying artifacts not in rescue table."""
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1', 'PB.4.1'],
            'filter_result': ['Isoform', 'Artifact', 'Artifact', 'Artifact'],
            'associated_gene': ['GENE1', 'GENE1', 'GENE2', 'novel_gene_1']
        })
        
        rescue_df = pd.DataFrame({
            'artifact': ['PB.2.1'],
            'assigned_transcript': ['PB.1.1']
        })
        
        result = get_unrescued_artifacts(classif_df, rescue_df)
        
        # Should identify PB.3.1 and PB.4.1 as unrescued
        assert len(result) == 2
        assert set(result['artifact']) == {'PB.3.1', 'PB.4.1'}
        
        # Check TD assignments
        assert result[result['artifact'] == 'PB.3.1']['assigned_transcript'].values[0] == 'GENE2_TD'
        assert result[result['artifact'] == 'PB.4.1']['assigned_transcript'].values[0] == 'general_TD'
    
    def test_all_rescued(self):
        """Test when all artifacts are rescued."""
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'filter_result': ['Isoform', 'Artifact'],
            'associated_gene': ['GENE1', 'GENE1']
        })
        
        rescue_df = pd.DataFrame({
            'artifact': ['PB.2.1'],
            'assigned_transcript': ['PB.1.1']
        })
        
        result = get_unrescued_artifacts(classif_df, rescue_df)
        
        assert len(result) == 0
    
    def test_no_artifacts(self):
        """Test when there are no artifacts."""
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'filter_result': ['Isoform', 'Isoform'],
            'associated_gene': ['GENE1', 'GENE2']
        })
        
        rescue_df = pd.DataFrame({
            'artifact': [],
            'assigned_transcript': []
        })
        
        result = get_unrescued_artifacts(classif_df, rescue_df)
        
        assert len(result) == 0

# Perhaps this can be deleted or merged with the test in test_requantification.py
class TestBuildArtifactTable:
    """Test suite for build_artifact_table function."""
    
    def test_combine_rescued_and_unrescued(self):
        """Test combining rescued and unrescued artifacts."""
        rescue_df = pd.DataFrame({
            'artifact': ['PB.2.1', 'PB.3.1'],
            'assigned_transcript': ['PB.1.1', 'PB.1.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1', 'PB.4.1'],
            'filter_result': ['Isoform', 'Artifact', 'Artifact', 'Artifact'],
            'associated_gene': ['GENE1', 'GENE1', 'GENE1', 'GENE2']
        })
        
        result = build_artifact_table(rescue_df, classif_df)
        
        # Should have all 3 artifacts
        assert len(result) == 3
        assert 'isoform' in result.columns
        assert 'assigned_transcript' in result.columns
        
        # Check that PB.4.1 was assigned to GENE2_TD
        unrescued_row = result[result['isoform'] == 'PB.4.1']
        assert unrescued_row['assigned_transcript'].values[0] == 'GENE2_TD'


class TestPrepareCountMatrices:
    """Test suite for prepare_count_matrices function."""
    
    def test_basic_matrix_split(self):
        """Test splitting counts into base and source matrices."""
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1'],
            'sample1': [100, 50, 25],
            'sample2': [200, 75, 30]
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact']
        })
        
        base_df, source_df = prepare_count_matrices(old_counts, classif_df)
        
        # Check base (valid isoforms)
        assert len(base_df) == 2
        assert 'PB.1.1' in base_df.index
        assert 'PB.2.1' in base_df.index
        assert base_df.loc['PB.1.1', 'sample1'] == 100
        
        # Check source (artifacts)
        assert len(source_df) == 1
        assert 'PB.3.1' in source_df.index
        assert source_df.loc['PB.3.1', 'sample1'] == 25
    
    def test_missing_isoforms_filled_with_zero(self):
        """Test that missing isoforms are filled with 0."""
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1'],
            'sample1': [100]
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'filter_result': ['Isoform', 'Artifact']
        })
        
        _, source_df = prepare_count_matrices(old_counts, classif_df)
        
        # PB.2.1 should be in source with 0 counts
        assert 'PB.2.1' in source_df.index
        assert source_df.loc['PB.2.1', 'sample1'] == 0


class TestCalculateDistributionFractions:
    """Test suite for calculate_distribution_fractions function."""
    
    def test_proportional_distribution(self):
        """Test proportional distribution based on target abundance."""
        exploded_map = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1'],
            'assigned_transcript': ['target1', 'target2']
        })
        
        base_df = pd.DataFrame({
            'sample1': [100, 200]
        }, index=['target1', 'target2'])
        
        sample_cols = ['sample1']
        
        result = calculate_distribution_fractions(exploded_map, base_df, sample_cols)
        
        # Should split proportionally: 100/(100+200) = 1/3, 200/(100+200) = 2/3
        assert len(result) == 2
        assert pytest.approx(result.loc[0, 'sample1'], abs=0.01) == 1/3
        assert pytest.approx(result.loc[1, 'sample1'], abs=0.01) == 2/3
    
    def test_uniform_distribution_when_all_zero(self):
        """Test uniform distribution when all targets have zero counts."""
        exploded_map = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1'],
            'assigned_transcript': ['target1', 'target2']
        })
        
        base_df = pd.DataFrame({
            'sample1': [0, 0]
        }, index=['target1', 'target2'])
        
        sample_cols = ['sample1']
        
        result = calculate_distribution_fractions(exploded_map, base_df, sample_cols)
        
        # Should split uniformly: 1/2 each
        assert pytest.approx(result.loc[0, 'sample1'], abs=0.01) == 0.5
        assert pytest.approx(result.loc[1, 'sample1'], abs=0.01) == 0.5
    
    def test_multiple_samples(self):
        """Test distribution with multiple samples."""
        exploded_map = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1'],
            'assigned_transcript': ['target1', 'target2']
        })
        
        base_df = pd.DataFrame({
            'sample1': [100, 200],
            'sample2': [50, 150]
        }, index=['target1', 'target2'])
        
        sample_cols = ['sample1', 'sample2']
        
        result = calculate_distribution_fractions(exploded_map, base_df, sample_cols)
        
        # Sample1: 1/3 and 2/3
        assert pytest.approx(result.loc[0, 'sample1'], abs=0.01) == 1/3
        assert pytest.approx(result.loc[1, 'sample1'], abs=0.01) == 2/3
        
        # Sample2: 50/200=1/4 and 150/200=3/4
        assert pytest.approx(result.loc[0, 'sample2'], abs=0.01) == 0.25
        assert pytest.approx(result.loc[1, 'sample2'], abs=0.01) == 0.75


class TestDistributeIntegerCounts:
    """Test suite for distribute_integer_counts function."""
    
    def test_integer_conservation(self):
        """Test that integer counts are conserved after distribution."""
        source_df = pd.DataFrame({
            'sample1': [10]
        }, index=['artifact1'])
        
        exploded_map = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1', 'artifact1'],
            'assigned_transcript': ['target1', 'target2', 'target3']
        })
        
        fractions = pd.DataFrame({
            'sample1': [0.33, 0.33, 0.34]  # Adds up to 1.0
        })
        
        sample_cols = ['sample1']
        
        result = distribute_integer_counts(source_df, exploded_map, fractions, sample_cols)
        
        # Total should still be 10
        assert result['sample1'].sum() == 10
    
    def test_remainder_goes_to_first_target(self):
        """Test that remainder is added to first target."""
        source_df = pd.DataFrame({
            'sample1': [10]
        }, index=['artifact1'])
        
        exploded_map = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1'],
            'assigned_transcript': ['target1', 'target2']
        })
        
        # Equal split: 10 * 0.5 = 5.0 each (no remainder expected with equal split)
        fractions = pd.DataFrame({
            'sample1': [1/3, 2/3]  # 10 * 1/3 = 3.33 -> 3, 10 * 2/3 = 6.67 -> 6, remainder = 1
        })
        
        sample_cols = ['sample1']
        
        result = distribute_integer_counts(source_df, exploded_map, fractions, sample_cols)
        
        # Should have target1=4 (3+1 remainder), target2=6
        assert result.loc['target1', 'sample1'] == 4
        assert result.loc['target2', 'sample1'] == 6
        assert result['sample1'].sum() == 10


class TestExportCounts:
    """Test suite for export_counts function."""
    
    def test_basic_export(self, tmp_path):
        """Test basic count export functionality."""
        old_counts_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1'],
            'sample1': [100, 50, 25],
            'sample2': [200, 75, 30]
        })
        
        new_counts_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'GENE1_TD'],
            'sample1': [100, 75, 0],
            'sample2': [200, 105, 0]
        })
        
        prefix = str(tmp_path / "test_output")
        
        result = export_counts(old_counts_df, new_counts_df, prefix)
        
        # Check that files were created
        assert os.path.exists(f"{prefix}_reassigned_counts.tsv")
        assert os.path.exists(f"{prefix}_reassigned_counts_extended.tsv")
        
        # Check result dataframe
        assert len(result) == 3
        assert 'isoform' in result.columns
        assert set(result['isoform']) == {'PB.1.1', 'PB.2.1', 'GENE1_TD'}
    
    def test_extended_table_contains_all_isoforms(self, tmp_path):
        """Test that extended table contains all isoforms with old and new counts."""
        old_counts_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1'],
            'sample1': [100, 50, 25]
        })
        
        new_counts_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'GENE1_TD'],
            'sample1': [100, 75, 0]
        })
        
        prefix = str(tmp_path / "test_output")
        
        export_counts(old_counts_df, new_counts_df, prefix)
        
        # Read extended table
        extended_df = pd.read_csv(f"{prefix}_reassigned_counts_extended.tsv", sep='\t')
        
        # Should have all 4 unique isoforms
        assert len(extended_df) == 4
        assert set(extended_df['isoform']) == {'PB.1.1', 'PB.2.1', 'PB.3.1', 'GENE1_TD'}
        
        # Check column structure
        assert 'old_sample1' in extended_df.columns
        assert 'new_sample1' in extended_df.columns
        
        # Check specific values
        pb31_row = extended_df[extended_df['isoform'] == 'PB.3.1'].iloc[0]
        assert pb31_row['old_sample1'] == 25
        assert pb31_row['new_sample1'] == 0  # Artifact lost its counts
        
        gene1_td_row = extended_df[extended_df['isoform'] == 'GENE1_TD'].iloc[0]
        assert gene1_td_row['old_sample1'] == 0  # New isoform
        assert gene1_td_row['new_sample1'] == 0


class TestCalculateTPM:
    """Test suite for calculate_tpm function."""
    
    def test_basic_tpm_calculation(self):
        """Test basic TPM calculation."""
        counts = pd.Series([100, 200, 300])
        lengths = pd.Series([1000, 2000, 3000])
        
        result = calculate_tpm(counts, lengths)
        
        # TPM calculation:
        # RPK = counts / (length/1000)
        # RPK = [100/1, 200/2, 300/3] = [100, 100, 100]
        # Total RPK = 300
        # TPM = (RPK / Total RPK) * 1e6
        # TPM = [100/300, 100/300, 100/300] * 1e6 = [333333.33, 333333.33, 333333.33]
        
        expected = pd.Series([1e6/3, 1e6/3, 1e6/3])
        pd.testing.assert_series_equal(result, expected, check_exact=False, rtol=0.01)
    
    def test_tpm_sum_to_million(self):
        """Test that TPM values sum to 1 million."""
        counts = pd.Series([50, 100, 150, 200])
        lengths = pd.Series([1000, 2000, 1500, 2500])
        
        result = calculate_tpm(counts, lengths)
        
        # Sum should be approximately 1 million
        assert pytest.approx(result.sum(), abs=1) == 1e6
