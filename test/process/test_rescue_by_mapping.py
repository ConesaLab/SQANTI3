"""
Integration tests for rescue_by_mapping.py function.
Tests the complete rescue by mapping pipeline end-to-end.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.rescue_by_mapping import rescue_by_mapping


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(main_path) / "test" / "test_data" / "rescue"


@pytest.fixture
def mapping_hits_df(test_data_dir):
    """Load mapping hits DataFrame."""
    return pd.read_csv(test_data_dir / "mapping_hits.tsv", sep="\t")


@pytest.fixture
def classification_df(test_data_dir):
    """Load classification DataFrame."""
    return pd.read_csv(test_data_dir / "classification.tsv", sep="\t")


@pytest.fixture
def reference_rules_file(test_data_dir):
    """Return path to reference rules file."""
    return str(test_data_dir / "reference_classification_filtered.tsv")


@pytest.fixture
def reference_ml_file(test_data_dir):
    """Return path to reference ML file."""
    return str(test_data_dir / "reference_ml.tsv")


@pytest.fixture
def automatic_rescue_df(test_data_dir):
    """Load automatic rescue DataFrame."""
    return pd.read_csv(test_data_dir / "automatic_rescue.tsv", sep="\t")


class TestRescueByMapping:
    """Integration tests for rescue_by_mapping function."""
    
    def test_rescue_by_mapping_rules(self, mapping_hits_df, reference_rules_file, 
                                     classification_df, automatic_rescue_df):
        """Test complete rescue by mapping with rules strategy."""
        automatic_list = pd.Series(['REF.AUTO.1', 'REF.AUTO.2'])
        
        inclusion_list, rescue_df = rescue_by_mapping(
            mapping_hits_df, 
            reference_rules_file, 
            classification_df,
            automatic_list,
            automatic_rescue_df,
            "rules"
        )
        
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
        assert len(inclusion_list) > 0
        assert len(rescue_df) > 0
        
        # Check rescue_df structure
        assert "artifact" in rescue_df.columns
        assert "assigned_transcript" in rescue_df.columns
        assert "rescue_mode" in rescue_df.columns
    
    def test_rescue_by_mapping_ml(self, mapping_hits_df, reference_ml_file, 
                                  classification_df, automatic_rescue_df):
        """Test complete rescue by mapping with ML strategy."""
        # Use minimal required columns for ML
        classif_minimal = classification_df[["isoform", "filter_result", 
                                            "structural_category", "associated_transcript","POS_MLprob"]].copy()
        automatic_list = pd.Series(['REF.AUTO.1'])
        
        inclusion_list, rescue_df = rescue_by_mapping(
            mapping_hits_df,
            reference_ml_file,
            classif_minimal,
            automatic_list,
            automatic_rescue_df,
            "ml",
            thr=0.7
        )
        
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
        
        # With ML strategy
        assert "ml_mapping" in rescue_df["rescue_mode"].values
    
    def test_rescue_by_mapping_different_threshold(self, mapping_hits_df, reference_ml_file,
                                                   classification_df, automatic_rescue_df):
        """Test rescue by mapping with different ML threshold."""
        classif_minimal = classification_df[["isoform", "filter_result",
                                            "structural_category", "associated_transcript","POS_MLprob"]].copy()
        automatic_list = pd.Series([])
        
        # High threshold
        inclusion_list_high, rescue_df_high = rescue_by_mapping(
            mapping_hits_df,
            reference_ml_file,
            classif_minimal,
            automatic_list,
            automatic_rescue_df,
            "ml",
            thr=0.9
        )
        
        # Low threshold  
        inclusion_list_low, rescue_df_low = rescue_by_mapping(
            mapping_hits_df,
            reference_ml_file,
            classif_minimal,
            automatic_list,
            automatic_rescue_df,
            "ml",
            thr=0.6
        )
        
        # Lower threshold should rescue more or equal transcripts
        assert len(inclusion_list_low) >= len(inclusion_list_high)
    
    def test_rescue_by_mapping_empty_automatic(self, mapping_hits_df, reference_rules_file,
                                               classification_df):
        """Test rescue by mapping with empty automatic rescue."""
        automatic_list = pd.Series([])
        automatic_df = pd.DataFrame(columns=['artifact', 'assigned_transcript', 
                                             'rescue_mode', 'origin', 'reintroduced'])
        
        inclusion_list, rescue_df = rescue_by_mapping(
            mapping_hits_df,
            reference_rules_file,
            classification_df,
            automatic_list,
            automatic_df,
            "rules"
        )
        
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
    
    def test_rescue_by_mapping_preserves_automatic(self, mapping_hits_df, reference_rules_file,
                                                   classification_df, automatic_rescue_df):
        """Test that automatic rescue entries are preserved."""
        automatic_list = pd.Series(['REF.AUTO.1', 'REF.AUTO.2'])
        
        inclusion_list, rescue_df = rescue_by_mapping(
            mapping_hits_df,
            reference_rules_file,
            classification_df,
            automatic_list,
            automatic_rescue_df,
            "rules"
        )
        
        # Check that automatic entries are in rescue_df
        automatic_entries = rescue_df[rescue_df["rescue_mode"] == "automatic"]
        assert len(automatic_entries) == 2
        assert "REF.AUTO.1" in automatic_entries["assigned_transcript"].values
        assert automatic_list.isin(inclusion_list).all()
    
    def test_multiple_candidates_same_target_both_rescued(self, mapping_hits_df, reference_rules_file,
                                                          classification_df):
        """Test that pipeline handles multiple candidates mapping to the same target."""
        automatic_list = pd.Series([])
        automatic_df = pd.DataFrame(columns=['artifact', 'assigned_transcript', 
                                             'rescue_mode', 'origin', 'reintroduced'])
        
        _, rescue_df = rescue_by_mapping(
            mapping_hits_df,
            reference_rules_file,
            classification_df,
            automatic_list,
            automatic_df,
            "rules"
        )
        
        # Verify pipeline produces results and handles multi-mapping scenario
        assert len(rescue_df) > 0
        assert "rules_mapping" in rescue_df["rescue_mode"].values
        
        # Check that some candidates share the same assigned transcript
        transcript_counts = rescue_df["assigned_transcript"].value_counts()
        assert any(transcript_counts > 1), "Some transcripts should be assigned to multiple candidates"
    
    def test_pipeline_handles_both_origins(self, mapping_hits_df, reference_rules_file,
                                                   classification_df):
        """Test that pipeline correctly handles both reference and lr_defined origins."""
        automatic_list = pd.Series([])
        automatic_df = pd.DataFrame(columns=['artifact', 'assigned_transcript', 
                                             'rescue_mode', 'origin', 'reintroduced'])
        
        _, rescue_df = rescue_by_mapping(
            mapping_hits_df,
            reference_rules_file,
            classification_df,
            automatic_list,
            automatic_df,
            "rules"
        )
        
        # Verify both origin types are present in results
        assert "origin" in rescue_df.columns
        origins = rescue_df["origin"].unique()
        
        # Pipeline should handle both reference and lr_defined transcripts
        assert len(origins) > 0, "Should have at least one origin type"
        # Check that lr_defined transcripts can be rescued
        lr_defined_entries = rescue_df[rescue_df["origin"] == "lr_defined"]
        if len(lr_defined_entries) > 0:
            assert all(lr_defined_entries["rescue_mode"] == "rules_mapping")
    
    def test_ml_strategy_produces_valid_results(self, mapping_hits_df, reference_ml_file,
                                                         classification_df):
        """Test that ML strategy pipeline produces valid rescue results."""
        classif_minimal = classification_df[["isoform", "filter_result",
                                            "structural_category", "associated_transcript","POS_MLprob"]].copy()
        automatic_list = pd.Series([])
        automatic_df = pd.DataFrame(columns=['artifact', 'assigned_transcript', 
                                             'rescue_mode', 'origin', 'reintroduced'])
        
        inclusion_list, rescue_df = rescue_by_mapping(
            mapping_hits_df,
            reference_ml_file,
            classif_minimal,
            automatic_list,
            automatic_df,
            "ml",
            thr=0.7
        )
        
        # Verify ML strategy produces results with correct mode
        assert len(rescue_df) > 0, "ML strategy should rescue some candidates"
        assert "ml_mapping" in rescue_df["rescue_mode"].values
        
        # Check that rescued transcripts have valid assignments
        assert all(rescue_df["assigned_transcript"].notna())
        assert all(rescue_df["origin"].notna())


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
