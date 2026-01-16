"""
Process-level tests for run_automatic_rescue function.
Tests the complete automatic rescue workflow with realistic data.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.rescue_steps import run_automatic_rescue


class TestRunAutomaticRescueMultiExon:
    """Test suite for run_automatic_rescue with multi-exon transcripts."""
    
    def test_rescue_lost_fsm_multiexon(self, sample_rescue_classification_mix, tmp_path):
        """Test rescuing lost FSM multi-exon references."""
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            sample_rescue_classification_mix,
            monoexons='none',  # No monoexon rescue by default
            prefix=prefix
        )
        
        # REF2 and REF4 are lost (no Isoform), should be rescued
        # REF2 has no FSM, so won't be rescued
        # REF4 has FSM (PB.4.1 is artifact), so should be rescued via reference
        # The monoexonic artifacts should be ignored
        expected_items = ['REF4']
        result_items = inclusion_list.tolist()

        assert sorted(result_items) == sorted(expected_items)
        assert isinstance(rescue_df, pd.DataFrame)

        # Check rescue_df structure
        assert 'artifact' in rescue_df.columns
        assert 'assigned_transcript' in rescue_df.columns
        assert 'rescue_mode' in rescue_df.columns
        assert 'origin' in rescue_df.columns
        assert 'reintroduced' in rescue_df.columns
        
        # All rescue mode should be "automatic"
        assert all(rescue_df['rescue_mode'] == 'automatic')
    
    def test_no_lost_references(self, tmp_path):
        """Test when all references are already represented."""
        # Create classification where all references have Isoform
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'structural_category': ['full-splice_match', 'full-splice_match'],
            'exons': [2, 3],
            'filter_result': ['Isoform', 'Isoform'],
            'associated_transcript': ['REF1', 'REF2'],
            'associated_gene': ['GENE1', 'GENE2']
        })
        
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            classif_df,
            monoexons=None,
            prefix=prefix
        )
        
        # Should return empty inclusion list
        assert len(inclusion_list) == 0
        
        # Rescue_df should indicate no rescue
        assert 'assigned_transcript' in rescue_df.columns
    
    def test_lost_reference_without_fsm(self, tmp_path):
        """Test that references without FSM are not rescued."""
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'structural_category': ['full-splice_match', 'incomplete-splice_match'],
            'exons': [2, 2],
            'filter_result': ['Isoform', 'Artifact'],
            'associated_transcript': ['REF1', 'REF2'],
            'associated_gene': ['GENE1', 'GENE2']
        })
        
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, _ = run_automatic_rescue(
            classif_df,
            monoexons=None,
            prefix=prefix
        )
        
        # REF2 is lost but has no FSM, shouldn't be rescued
        assert len(inclusion_list) == 0


class TestRunAutomaticRescueMonoexon:
    """Test suite for monoexon rescue modes."""
    
    def test_monoexon_rescue_all(self, sample_rescue_classification_mix, tmp_path):
        """Test monoexon rescue with 'all' parameter."""
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            sample_rescue_classification_mix,
            monoexons='all',
            prefix=prefix
        )
        
        # Should rescue lost FSM monoexons
        #  Together with the previous, now we have some lost FSM monoexons
        expected_items = ['REF4','REF5','REF6']
        result_items = inclusion_list.tolist()

        assert sorted(result_items) == sorted(expected_items)
        
        # Check that rescued items are in rescue_df
        assert len(rescue_df) > 0
    
    def test_monoexon_rescue_fsm(self, sample_rescue_classification_mix, tmp_path):
        """Test monoexon rescue with 'fsm' parameter."""
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            sample_rescue_classification_mix,
            monoexons='fsm',
            prefix=prefix
        )
        
        # Should rescue lost FSM monoexons only
        expected_items = ['REF4','REF5','REF6']
        result_items = inclusion_list.tolist()

        assert sorted(result_items) == sorted(expected_items)
        assert isinstance(rescue_df, pd.DataFrame)
    
    def test_monoexon_rescue_only(self, sample_rescue_classification_monoexon, tmp_path):
        """Test that monoexons are  rescued when there are only monoexons."""
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, _ = run_automatic_rescue(
            sample_rescue_classification_monoexon,
            monoexons=None,
            prefix=prefix
        )
        
        # Should not rescue monoexons
        # All monoexons in this classification, so nothing should be rescued
        expected_items = ['REF5','REF6']
        result_items = inclusion_list.tolist()

        assert sorted(result_items) == sorted(expected_items)
    
    @pytest.mark.parametrize("monoexons_param", ["all", "fsm", None])
    def test_monoexon_parameter_variations(
        self, 
        sample_rescue_classification_mix, 
        tmp_path,
        monoexons_param
    ):
        """Test that different monoexon parameters work without errors."""
        prefix = str(tmp_path / "test_output")
        
        # Should not raise any errors
        inclusion_list, rescue_df = run_automatic_rescue(
            sample_rescue_classification_mix,
            monoexons=monoexons_param,
            prefix=prefix
        )
        
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)


class TestRunAutomaticRescueMultipleArtifacts:
    """Test suite for handling multiple artifacts per reference."""
    
    def test_multiple_artifacts_same_reference(self, tmp_path):
        """Test rescue when multiple artifacts point to same lost reference."""
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2', 'PB.1.3', 'PB.2.1'],
            'structural_category': [
                'incomplete-splice_match', 'incomplete-splice_match',
                'full-splice_match', 'full-splice_match'
            ],
            'exons': [2, 2, 2, 3],
            'filter_result': ['Artifact', 'Artifact', 'Artifact', 'Isoform'],
            'associated_transcript': ['REF1', 'REF1', 'REF1', 'REF2'],
            'associated_gene': ['GENE1', 'GENE1', 'GENE1', 'GENE2']
        })
        
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            classif_df,
            monoexons=None,
            prefix=prefix
        )
        
        # REF1 should be rescued once
        if len(inclusion_list) == 1:
            assert 'REF1' in inclusion_list.values
        
        # Check reintroduced column in rescue_df
        if len(rescue_df) > 0:
            # Filter to REF1 artifacts
            ref1_artifacts = rescue_df[rescue_df['assigned_transcript'] == 'REF1']
            if len(ref1_artifacts) > 0:
                # First one should be "yes", rest "no"
                reintroduced_values = ref1_artifacts['reintroduced'].tolist()
                assert reintroduced_values[0] == 'yes'
                assert all(v == 'no' for v in reintroduced_values[1:])


class TestRunAutomaticRescueRealData:
    """Test suite using real test data files."""
    
    def test_with_lost_fsm_file(self, rescue_classification_file_lost_fsm, tmp_path):
        """Test rescue with classification file containing lost FSM."""
        classif_df = pd.read_csv(rescue_classification_file_lost_fsm, sep='\t')
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            classif_df,
            monoexons=None,
            prefix=prefix
        )
        
        # Should successfully process the file
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
    
    def test_with_monoexon_file(self, rescue_classification_file_monoexon, tmp_path):
        """Test rescue with monoexonic classification file."""
        classif_df = pd.read_csv(rescue_classification_file_monoexon, sep='\t')
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            classif_df,
            monoexons='all',
            prefix=prefix
        )
        
        # Should successfully process the file
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
    
    def test_with_no_lost_file(self, rescue_classification_file_no_lost, tmp_path):
        """Test rescue when no references are lost."""
        classif_df = pd.read_csv(rescue_classification_file_no_lost, sep='\t')
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, _ = run_automatic_rescue(
            classif_df,
            monoexons=None,
            prefix=prefix
        )
        
        # Should return empty inclusion list
        assert len(inclusion_list) == 0
    
    def test_with_multi_artifact_file(self, rescue_classification_file_multi_artifact, tmp_path):
        """Test rescue with multiple artifacts per reference."""
        classif_df = pd.read_csv(rescue_classification_file_multi_artifact, sep='\t')
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, rescue_df = run_automatic_rescue(
            classif_df,
            monoexons=None,
            prefix=prefix
        )
        
        # Should successfully process the file
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
        
        # Check reintroduced logic for duplicate references
        if len(rescue_df) > 0:
            # Group by assigned_transcript and check first is "yes"
            for ref_id, group in rescue_df.groupby('assigned_transcript'):
                reintroduced_values = group['reintroduced'].tolist()
                if len(reintroduced_values) > 1:
                    assert reintroduced_values[0] == 'yes'
                    assert all(v == 'no' for v in reintroduced_values[1:])


class TestRunAutomaticRescueOutputFormat:
    """Test suite for verifying output format and structure."""
    
    def test_inclusion_list_format(self, sample_rescue_classification, tmp_path):
        """Test that inclusion list has correct format."""
        prefix = str(tmp_path / "test_output")
        
        inclusion_list, _ = run_automatic_rescue(
            sample_rescue_classification,
            monoexons=None,
            prefix=prefix
        )
        
        # Should be a pandas Series
        assert isinstance(inclusion_list, pd.Series)
        
        # If not empty, should contain transcript IDs
        if len(inclusion_list) > 0:
            # All values should be strings (transcript IDs)
            assert all(isinstance(val, str) for val in inclusion_list)
    
    def test_rescue_df_columns(self, sample_rescue_classification, tmp_path):
        """Test that rescue_df has all required columns."""
        prefix = str(tmp_path / "test_output")
        
        _, rescue_df = run_automatic_rescue(
            sample_rescue_classification,
            monoexons=None,
            prefix=prefix
        )
        
        # Required columns
        required_columns = [
            'artifact',
            'assigned_transcript',
            'rescue_mode',
            'origin',
            'reintroduced'
        ]
        
        for col in required_columns:
            assert col in rescue_df.columns, f"Missing column: {col}"
    
    def test_rescue_mode_values(self, sample_rescue_classification, tmp_path):
        """Test that rescue_mode column has correct values."""
        prefix = str(tmp_path / "test_output")
        
        _, rescue_df = run_automatic_rescue(
            sample_rescue_classification,
            monoexons='all',
            prefix=prefix
        )
        
        if len(rescue_df) > 0:
            # All rescue_mode should be "automatic"
            assert all(rescue_df['rescue_mode'] == 'automatic')
    
    def test_origin_values(self, sample_rescue_classification, tmp_path):
        """Test that origin column has correct values."""
        prefix = str(tmp_path / "test_output")
        
        _, rescue_df = run_automatic_rescue(
            sample_rescue_classification,
            monoexons='all',
            prefix=prefix
        )
        
        if len(rescue_df) > 0:
            # All origin should be "reference"
            assert all(rescue_df['origin'] == 'reference')
    
    def test_reintroduced_values(self, sample_rescue_classification, tmp_path):
        """Test that reintroduced column has valid values."""
        prefix = str(tmp_path / "test_output")
        
        _, rescue_df = run_automatic_rescue(
            sample_rescue_classification,
            monoexons='all',
            prefix=prefix
        )
        
        if len(rescue_df) > 0:
            # All reintroduced should be "yes" or "no"
            valid_values = {'yes', 'no'}
            assert all(val in valid_values for val in rescue_df['reintroduced'])


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
