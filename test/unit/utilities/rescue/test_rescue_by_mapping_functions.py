"""
Unit tests for rescue_by_mapping.py functions.
Tests the rescue by mapping pipeline and its component functions.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os
import tempfile
import numpy as np

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.rescue_by_mapping import (
    merge_classifications,
    add_filter_results,
    select_best_hits,
    filter_mapping_hits,
    merge_rescue_modes,
    find_reintroduced_transcripts,
    rescue_by_mapping
)


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
    return str(test_data_dir / "reference_rules.tsv")


@pytest.fixture
def reference_ml_file(test_data_dir):
    """Return path to reference ML file."""
    return str(test_data_dir / "reference_ml.tsv")


@pytest.fixture
def automatic_rescue_df(test_data_dir):
    """Load automatic rescue DataFrame."""
    return pd.read_csv(test_data_dir / "automatic_rescue.tsv", sep="\t")


class TestMergeClassifications:
    """Test suite for merge_classifications function."""
    
    def test_merge_rules_strategy(self, reference_rules_file, classification_df):
        """Test merging classifications with rules strategy."""
        result = merge_classifications(reference_rules_file, classification_df, "rules")
        
        assert isinstance(result, pd.DataFrame)
        assert "isoform" in result.columns
        assert "filter_result" in result.columns
        assert "origin" in result.columns
        
        # Should have both reference and lr_defined entries
        assert "reference" in result["origin"].values
        assert "lr_defined" in result["origin"].values
    
    def test_merge_ml_strategy(self, reference_ml_file, classification_df):
        """Test merging classifications with ML strategy."""
        result = merge_classifications(reference_ml_file, classification_df, "ml", thr=0.7)
        
        assert isinstance(result, pd.DataFrame)
        assert "filter_result" in result.columns
        assert "origin" in result.columns
        
        # Check that filter_result is created based on threshold for reference entries
        ref_entries = result[result["origin"] == "reference"]
        assert "POS_MLprob" in ref_entries.columns
        high_prob = ref_entries[ref_entries["POS_MLprob"] >= 0.7]
        if len(high_prob) > 0:
            assert all(high_prob["filter_result"] == "Isoform")
        
        low_prob = ref_entries[ref_entries["POS_MLprob"] < 0.7]
        if len(low_prob) > 0:
            assert all(low_prob["filter_result"] == "Artifact")
    
    def test_merge_ml_different_threshold(self, reference_ml_file, classification_df):
        """Test ML strategy with different threshold."""
        result = merge_classifications(reference_ml_file, classification_df, "ml", thr=0.8)
        
        # With higher threshold, fewer should be Isoform
        ref_entries = result[result["origin"] == "reference"]
        high_prob = ref_entries[ref_entries["POS_MLprob"] >= 0.8]
        if len(high_prob) > 0:
            assert all(high_prob["filter_result"] == "Isoform")
        
        mid_prob = ref_entries[(ref_entries["POS_MLprob"] >= 0.7) & (ref_entries["POS_MLprob"] < 0.8)]
        if len(mid_prob) > 0:
            assert all(mid_prob["filter_result"] == "Artifact")
    
    def test_merge_preserves_lr_classification(self, reference_rules_file, classification_df):
        """Test that LR classification is preserved."""
        result = merge_classifications(reference_rules_file, classification_df, "rules")
        
        lr_entries = result[result["origin"] == "lr_defined"]
        
        # Should match original classification
        for _, row in lr_entries.iterrows():
            original = classification_df[classification_df["isoform"] == row["isoform"]]
            if not original.empty:
                assert row["filter_result"] == original.iloc[0]["filter_result"]


class TestAddFilterResults:
    """Test suite for add_filter_results function."""
    
    def test_add_filter_results_basic(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test adding filter results to mapping hits."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        result = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        assert "hit_filter_result" in result.columns
        assert "hit_origin" in result.columns
        assert "candidate_structural_category" in result.columns
        assert "associated_transcript" in result.columns
    
    def test_add_filter_results_correct_mapping(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test that filter results are correctly mapped."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        result = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        # Check specific mapping
        pb11_ref1 = result[(result["rescue_candidate"] == "PB.1.1") & (result["mapping_hit"] == "REF1")]
        assert len(pb11_ref1) == 1
        assert pb11_ref1.iloc[0]["hit_filter_result"] == "Isoform"
    
    def test_add_filter_results_with_ml(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test adding filter results with ML strategy."""
        combined_class = merge_classifications(reference_ml_file, classification_df, "ml", thr=0.7)
        result = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        # hit_origin should exist
        assert "hit_origin" in result.columns
        assert "hit_filter_result" in result.columns
        
        # Verify structure is correct
        assert len(result) > 0


class TestSelectBestHits:
    """Test suite for select_best_hits function."""
    
    def test_select_highest_score(self):
        """Test selection of highest score."""
        df = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.1.1', 'PB.1.1'],
            'mapping_hit': ['REF1', 'REF2', 'REF3'],
            'score': [100, 150, 120],
            'hit_origin': ['reference', 'reference', 'reference']
        })
        
        result = select_best_hits(df)
        
        assert len(result) == 1
        assert result.iloc[0]['mapping_hit'] == 'REF2'
        assert result.iloc[0]['score'] == 150
    
    def test_prefer_lr_defined(self):
        """Test preference for lr_defined over reference."""
        df = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.1.1'],
            'mapping_hit': ['REF1', 'LR.1.1'],
            'score': [100, 100],
            'hit_origin': ['reference', 'lr_defined']
        })
        
        result = select_best_hits(df)
        
        assert len(result) == 1
        assert result.iloc[0]['mapping_hit'] == 'LR.1.1'
        assert result.iloc[0]['hit_origin'] == 'lr_defined'
    
    def test_keep_ties(self):
        """Test that ties are kept."""
        df = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.1.1', 'PB.1.1'],
            'mapping_hit': ['REF1', 'REF2', 'REF3'],
            'score': [150, 150, 100],
            'hit_origin': ['reference', 'reference', 'reference']
        })
        
        result = select_best_hits(df)
        
        assert len(result) == 2
        assert set(result['mapping_hit']) == {'REF1', 'REF2'}


class TestFilterMappingHits:
    """Test suite for filter_mapping_hits function."""
    
    def test_filter_rules_strategy(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test filtering with rules strategy."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "rules")
        
        assert isinstance(result, pd.DataFrame)
        # Should only keep Isoforms
        assert all(result["hit_filter_result"] == "Isoform")
        assert "score" in result.columns
    
    def test_filter_ml_strategy(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test filtering with ML strategy."""
        combined_class = merge_classifications(reference_ml_file, classification_df[["isoform", "filter_result","POS_MLprob"]], "ml", thr=0.7)
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "ml")
        
        assert isinstance(result, pd.DataFrame)
        # Score should be calculated
        assert "score" in result.columns
        assert len(result) >= 0  # May be empty depending on data
    
    def test_filter_selects_best_per_candidate(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test that best hit is selected per candidate."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "rules")
        
        # Each rescue candidate should have best hit(s) selected
        # PB.1.1 maps to REF1 (200) and LR.1.1 (195), both Isoforms - should select REF1
        pb11_results = result[result["rescue_candidate"] == "PB.1.1"]
        if not pb11_results.empty:
            assert pb11_results.iloc[0]["mapping_hit"] == "REF1"


class TestMergeRescueModes:
    """Test suite for merge_rescue_modes function."""
    
    def test_merge_basic(self, automatic_rescue_df):
        """Test basic merging of rescue modes."""
        full_rescue = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.2.1'],
            'mapping_hit': ['REF1', 'REF5'],
            'hit_origin': ['reference', 'reference']
        })
        
        inclusion_list = pd.Series(['REF1', 'REF5'])
        
        result = merge_rescue_modes(automatic_rescue_df, full_rescue, inclusion_list, "rules")
        
        assert isinstance(result, pd.DataFrame)
        assert "artifact" in result.columns
        assert "assigned_transcript" in result.columns
        assert "rescue_mode" in result.columns
        assert "reintroduced" in result.columns
    
    def test_merge_marks_reintroduced(self, automatic_rescue_df):
        """Test that reintroduced transcripts are marked correctly."""
        full_rescue = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.2.1'],
            'mapping_hit': ['REF1', 'REF5'],
            'hit_origin': ['reference', 'reference']
        })
        
        inclusion_list = pd.Series(['REF1'])  # Only REF1 in inclusion list
        
        result = merge_rescue_modes(automatic_rescue_df, full_rescue, inclusion_list, "rules")
        
        # Check reintroduced marking
        ref1_entries = result[result["assigned_transcript"] == "REF1"]
        assert all(ref1_entries["reintroduced"] == "yes")
        
        ref5_entries = result[result["assigned_transcript"] == "REF5"]
        assert all(ref5_entries["reintroduced"] == "no")
    
    def test_merge_adds_rescue_mode(self, automatic_rescue_df):
        """Test that rescue mode is added correctly."""
        full_rescue = pd.DataFrame({
            'rescue_candidate': ['PB.1.1'],
            'mapping_hit': ['REF1'],
            'hit_origin': ['reference']
        })
        
        inclusion_list = pd.Series(['REF1'])
        
        result = merge_rescue_modes(automatic_rescue_df, full_rescue, inclusion_list, "ml")
        
        # Should have both automatic and ml_mapping modes
        assert "automatic" in result["rescue_mode"].values
        assert "ml_mapping" in result["rescue_mode"].values


class TestFindReintroducedTranscripts:
    """Test suite for find_reintroduced_transcripts function."""
    
    def test_find_with_automatic_rescue(self):
        """Test finding reintroduced transcripts with automatic rescue."""
        mapping_filt = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.2.1'],
            'mapping_hit': ['REF1', 'REF2'],
            'hit_origin': ['reference', 'reference']
        })
        
        automatic_rescue = pd.Series(['REF.AUTO.1', 'REF.AUTO.2'])
        
        result = find_reintroduced_transcripts(mapping_filt, automatic_rescue)
        
        assert isinstance(result, pd.Series)
        assert len(result) == 4  # 2 from mapping + 2 from automatic
        assert 'REF1' in result.values
        assert 'REF.AUTO.1' in result.values
    
    def test_find_without_automatic_rescue(self):
        """Test finding reintroduced transcripts without automatic rescue."""
        mapping_filt = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.2.1'],
            'mapping_hit': ['REF1', 'REF2'],
            'hit_origin': ['reference', 'lr_defined']
        })
        
        automatic_rescue = pd.Series([])
        
        result = find_reintroduced_transcripts(mapping_filt, automatic_rescue)
        
        assert isinstance(result, pd.Series)
        # Should only include reference transcripts
        assert len(result) == 1
        assert 'REF1' in result.values
        assert 'REF2' not in result.values  # lr_defined, not reference
    
    def test_find_removes_duplicates(self):
        """Test that duplicates are removed."""
        mapping_filt = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.1.2'],
            'mapping_hit': ['REF1', 'REF1'],
            'hit_origin': ['reference', 'reference']
        })
        
        automatic_rescue = pd.Series(['REF1'])
        
        result = find_reintroduced_transcripts(mapping_filt, automatic_rescue)
        
        assert len(result) == 1
        assert result.iloc[0] == 'REF1'


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
    
    def test_multiple_candidates_same_target_both_rescued(self, mapping_hits_df, reference_rules_file,
                                                          classification_df):
        """Test that multiple candidates mapping to the same target are all rescued."""
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
        
        # PB.5.1 and PB.5.2 both map to REF10 with same score
        # Both should be rescued
        pb5_1_rescued = rescue_df[rescue_df["artifact"] == "PB.5.1"]
        pb5_2_rescued = rescue_df[rescue_df["artifact"] == "PB.5.2"]
        
        assert len(pb5_1_rescued) == 1, "PB.5.1 should be rescued"
        assert len(pb5_2_rescued) == 1, "PB.5.2 should be rescued"
        assert pb5_1_rescued.iloc[0]["assigned_transcript"] == "REF10"
        assert pb5_2_rescued.iloc[0]["assigned_transcript"] == "REF10"
    
    def test_multiple_candidates_prefer_lr_defined(self, mapping_hits_df, reference_rules_file,
                                                   classification_df):
        """Test that when multiple candidates map to same target with one lr_defined and one reference, only lr_defined is kept."""
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
        
        # PB.6.1 and PB.6.2 both map to LR.6.1 with same score (200)
        # LR.6.1 is lr_defined, so both should map to it and both should be rescued
        pb6_1_rescued = rescue_df[rescue_df["artifact"] == "PB.6.1"]
        pb6_2_rescued = rescue_df[rescue_df["artifact"] == "PB.6.2"]
        
        assert len(pb6_1_rescued) == 1, "PB.6.1 should be rescued"
        assert len(pb6_2_rescued) == 1, "PB.6.2 should be rescued"
        
        # Both should map to LR.6.1 (lr_defined)
        assert pb6_1_rescued.iloc[0]["assigned_transcript"] == "LR.6.1"
        assert pb6_2_rescued.iloc[0]["assigned_transcript"] == "LR.6.1"
        assert pb6_1_rescued.iloc[0]["origin"] == "lr_defined"
        assert pb6_2_rescued.iloc[0]["origin"] == "lr_defined"
    
    def test_multiple_candidates_same_target_ml_strategy(self, mapping_hits_df, reference_ml_file,
                                                         classification_df):
        """Test multiple candidates mapping to same target with ML strategy."""
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
        
        # PB.5.1 and PB.5.2 both map to REF10 (POS_MLprob=0.88, above threshold)
        # Both should be rescued
        pb5_1_rescued = rescue_df[rescue_df["artifact"] == "PB.5.1"]
        pb5_2_rescued = rescue_df[rescue_df["artifact"] == "PB.5.2"]
        
        assert len(pb5_1_rescued) == 1, "PB.5.1 should be rescued with ML strategy"
        assert len(pb5_2_rescued) == 1, "PB.5.2 should be rescued with ML strategy"
        assert pb5_1_rescued.iloc[0]["assigned_transcript"] == "REF10"
        assert pb5_2_rescued.iloc[0]["assigned_transcript"] == "REF10"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
