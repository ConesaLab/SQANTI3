"""
Unit tests for automatic_rescue.py functions.
Tests individual functions with simple, focused test cases.
"""
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.automatic_rescue import (
    get_lost_reference_id,
    rescue_lost_reference,
    save_automatic_rescue
)


class TestGetLostReferenceId:
    """Test suite for get_lost_reference_id function."""
    
    def test_all_references_represented(self):
        """When all references have at least one Isoform, return empty array."""
        df = pd.DataFrame({
            'associated_transcript': ['REF1', 'REF1', 'REF2'],
            'filter_result': ['Isoform', 'Artifact', 'Isoform']
        })
        
        result = get_lost_reference_id(df)
        
        assert len(result) == 0
    
    def test_some_references_lost(self):
        """When some references have no Isoform, return those reference IDs."""
        df = pd.DataFrame({
            'associated_transcript': ['REF1', 'REF1', 'REF2', 'REF2'],
            'filter_result': ['Isoform', 'Artifact', 'Artifact', 'Artifact']
        })
        
        result = get_lost_reference_id(df)
        
        assert len(result) == 1
        assert 'REF2' in result
        assert 'REF1' not in result
    
    def test_all_references_lost(self):
        """When no references have Isoform, return all unique reference IDs."""
        df = pd.DataFrame({
            'associated_transcript': ['REF1', 'REF2', 'REF3'],
            'filter_result': ['Artifact', 'Artifact', 'Artifact']
        })
        
        result = get_lost_reference_id(df)
        
        assert len(result) == 3
        assert set(result) == {'REF1', 'REF2', 'REF3'}
    
    def test_multiple_artifacts_same_reference(self):
        """Multiple artifacts for same lost reference should return only once."""
        df = pd.DataFrame({
            'associated_transcript': ['REF1', 'REF1', 'REF1', 'REF2'],
            'filter_result': ['Artifact', 'Artifact', 'Artifact', 'Isoform']
        })
        
        result = get_lost_reference_id(df)
        
        assert len(result) == 1
        assert result[0] == 'REF1'


class TestRescueLostReference:
    """Test suite for rescue_lost_reference function."""
    
    def test_reference_with_fsm(self):
        """When reference has FSM, return the reference ID."""
        classif = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'associated_transcript': ['REF1', 'REF1'],
            'structural_category': ['full-splice_match', 'incomplete-splice_match']
        })
        
        result = rescue_lost_reference('REF1', classif)
        
        assert result is not None
        assert len(result) == 1
        assert result.iloc[0, 0] == 'REF1'
    
    def test_reference_without_fsm(self):
        """When reference has no FSM, return None."""
        classif = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2'],
            'associated_transcript': ['REF1', 'REF1'],
            'structural_category': ['incomplete-splice_match', 'novel_in_catalog']
        })
        
        result = rescue_lost_reference('REF1', classif)
        
        assert result is None
    
    def test_reference_not_in_classification(self):
        """When reference ID not in classification, return None."""
        classif = pd.DataFrame({
            'isoform': ['PB.1.1'],
            'associated_transcript': ['REF1'],
            'structural_category': ['full-splice_match']
        })
        
        result = rescue_lost_reference('REF2', classif)
        
        # Should return None or empty DataFrame
        assert result is None or len(result) == 0
 
class TestSaveAutomaticRescue:
    """Test suite for save_automatic_rescue function."""
    
    def test_normal_rescue(self, tmp_path):
        """Test normal rescue with valid inclusion list."""
        inclusion_df = pd.DataFrame({
            'isoform': ['REF1', 'REF2']
        })
        
        class_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2', 'PB.2.1'],
            'associated_transcript': ['REF1', 'REF1', 'REF2']
        })
        
        prefix = str(tmp_path / "test_output")
        
        result = save_automatic_rescue(inclusion_df, class_df, prefix)
        
        # Check result structure
        assert 'artifact' in result.columns
        assert 'assigned_transcript' in result.columns
        assert 'rescue_mode' in result.columns
        assert 'origin' in result.columns
        assert 'reintroduced' in result.columns
        
        # Check values
        assert all(result['rescue_mode'] == 'automatic')
        assert all(result['origin'] == 'reference')
        assert len(result) == 3  # Three artifacts matched
    
    def test_empty_inclusion_list(self, tmp_path):
        """Test with 'none' in inclusion list."""
        inclusion_df = pd.DataFrame({
            'associated_transcript': ['none']
        })
        
        class_df = pd.DataFrame({
            'isoform': ['PB.1.1'],
            'associated_transcript': ['REF1']
        })
        
        prefix = str(tmp_path / "test_output")
        
        result = save_automatic_rescue(inclusion_df, class_df, prefix)
        
        # Should handle 'none' gracefully
        assert 'isoform' in result.columns or 'artifact' in result.columns
    
    def test_reintroduced_column_logic(self, tmp_path):
        """Test that 'reintroduced' column correctly marks first occurrence."""
        inclusion_df = pd.DataFrame({
            'isoform': ['REF1']
        })
        
        class_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.1.2', 'PB.1.3'],
            'associated_transcript': ['REF1', 'REF1', 'REF1']
        })
        
        prefix = str(tmp_path / "test_output")
        
        result = save_automatic_rescue(inclusion_df, class_df, prefix)
        
        # First artifact should be marked "yes", rest "no"
        reintroduced_values = result['reintroduced'].tolist()
        assert reintroduced_values.count('yes') == 1
        assert reintroduced_values.count('no') == 2
        assert reintroduced_values[0] == 'yes'
    
    def test_multiple_different_references(self, tmp_path):
        """Test with multiple different rescued references."""
        inclusion_df = pd.DataFrame({
            'isoform': ['REF1', 'REF2', 'REF3']
        })
        
        class_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1', 'PB.3.2'],
            'associated_transcript': ['REF1', 'REF2', 'REF3', 'REF3']
        })
        
        prefix = str(tmp_path / "test_output")
        
        result = save_automatic_rescue(inclusion_df, class_df, prefix)
        
        # Check all references present
        assert len(result) == 4
        assert set(result['assigned_transcript'].unique()) == {'REF1', 'REF2', 'REF3'}
        
        # Each reference should have first occurrence marked "yes"
        yes_count = (result['reintroduced'] == 'yes').sum()
        assert yes_count == 3


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
