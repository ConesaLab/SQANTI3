"""
Unit tests for rescue_by_mapping.py functions.
Tests individual functions in the rescue by mapping module.
Integration tests are in test/process/test_rescue_by_mapping.py.
"""
import pytest
import pandas as pd
from pathlib import Path
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.rescue_by_mapping import (
    merge_classifications,
    add_filter_results,
    select_best_hits,
    filter_mapping_hits,
    merge_rescue_modes,
    find_reintroduced_transcripts
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
    """Return path to reference filtered classification file."""
    return str(test_data_dir / "reference_classification_filtered.tsv")


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
        assert "hit_POS_MLprob" not in pb11_ref1.columns  # rules strategy should not have this
    
    def test_add_filter_results_with_ml(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test adding filter results with ML strategy."""
        combined_class = merge_classifications(reference_ml_file, classification_df, "ml", thr=0.7)
        result = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        # hit_origin should exist
        assert "hit_origin" in result.columns
        assert "hit_filter_result" in result.columns
        assert "hit_POS_MLprob" in result.columns
        
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
    
    def test_filter_rules_strategy_basic(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test filtering with rules strategy - basic checks."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "rules")
        
        assert isinstance(result, pd.DataFrame)
        # Should only keep Isoforms
        assert all(result["hit_filter_result"] == "Isoform")
        assert "score" in result.columns
        # Score should equal alignment_score for rules strategy
        assert all(result["score"] == result["alignment_score"])
    
    def test_filter_rules_keeps_only_isoforms(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test that only mapping hits with Isoform filter_result are kept (rules strategy)."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "rules")
        
        # Based on reference_rules.tsv, these should be kept (Isoform):
        # REF1, REF2, REF4, REF6, REF7, REF8, LR.1.1, LR.2.1, REF10, LR.6.1
        # These should be excluded (Artifact): REF3, REF5, REF9, LR.3.1
        
        expected_isoforms = {"REF1", "REF2", "REF4", "REF6", "REF7", "REF8", 
                            "LR.1.1", "LR.2.1", "REF10", "LR.6.1"}
        excluded_artifacts = {"REF3", "REF5", "REF9", "LR.3.1"}
        
        result_hits = set(result["mapping_hit"].unique())
        
        # All results should be from expected isoforms
        assert result_hits.issubset(expected_isoforms)
        # No artifacts should be in results
        assert result_hits.isdisjoint(excluded_artifacts)
    
    def test_filter_rules_specific_candidates(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test specific rescue candidates with rules strategy."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "rules")
        
        # PB.1.1 maps to REF1 (200, Isoform) and LR.1.1 (195, Isoform)
        # Should select REF1 (higher score)
        pb11 = result[result["rescue_candidate"] == "PB.1.1"]
        assert len(pb11) == 1
        assert pb11.iloc[0]["mapping_hit"] == "REF1"
        assert pb11.iloc[0]["score"] == 200
        
        # PB.1.2 maps to REF2 (180, Isoform) and REF3 (180, Artifact)
        # Only REF2 should be kept (Artifact REF3 filtered out)
        pb12 = result[result["rescue_candidate"] == "PB.1.2"]
        assert len(pb12) == 1
        assert pb12.iloc[0]["mapping_hit"] == "REF2"
        assert pb12.iloc[0]["score"] == 180
        
        # PB.1.3 maps to REF4 (150, Isoform)
        # Should be kept
        pb13 = result[result["rescue_candidate"] == "PB.1.3"]
        assert len(pb13) == 1
        assert pb13.iloc[0]["mapping_hit"] == "REF4"
        
        # PB.2.1 maps to REF5 (190, Artifact) and LR.2.1 (190, Isoform)
        # Only LR.2.1 should be kept (REF5 is Artifact)
        pb21 = result[result["rescue_candidate"] == "PB.2.1"]
        assert len(pb21) == 1
        assert pb21.iloc[0]["mapping_hit"] == "LR.2.1"
        
        # PB.2.2 maps to REF6 (175, Isoform)
        # Should be kept
        pb22 = result[result["rescue_candidate"] == "PB.2.2"]
        assert len(pb22) == 1
        assert pb22.iloc[0]["mapping_hit"] == "REF6"
        
        # PB.3.1 maps to REF7 (160, Isoform)
        # Should be kept
        pb31 = result[result["rescue_candidate"] == "PB.3.1"]
        assert len(pb31) == 1
        assert pb31.iloc[0]["mapping_hit"] == "REF7"
        
        # PB.3.2 maps to LR.3.1 (185, Artifact) and REF8 (170, Isoform)
        # Only REF8 should be kept (LR.3.1 is Artifact)
        pb32 = result[result["rescue_candidate"] == "PB.3.2"]
        assert len(pb32) == 1
        assert pb32.iloc[0]["mapping_hit"] == "REF8"
        
        # PB.4.1 maps to REF9 (140, Artifact)
        # Should be filtered out completely
        pb41 = result[result["rescue_candidate"] == "PB.4.1"]
        assert len(pb41) == 0
        
        # PB.5.1 and PB.5.2 both map to REF10 (195, Isoform)
        # Both should be kept
        pb51 = result[result["rescue_candidate"] == "PB.5.1"]
        assert len(pb51) == 1
        assert pb51.iloc[0]["mapping_hit"] == "REF10"
        
        pb52 = result[result["rescue_candidate"] == "PB.5.2"]
        assert len(pb52) == 1
        assert pb52.iloc[0]["mapping_hit"] == "REF10"
        
        # PB.6.1 and PB.6.2 both map to LR.6.1 (200, Isoform)
        # Both should be kept
        pb61 = result[result["rescue_candidate"] == "PB.6.1"]
        assert len(pb61) == 1
        assert pb61.iloc[0]["mapping_hit"] == "LR.6.1"
        
        pb62 = result[result["rescue_candidate"] == "PB.6.2"]
        assert len(pb62) == 1
        assert pb62.iloc[0]["mapping_hit"] == "LR.6.1"
    
    def test_filter_rules_prefers_lr_defined(self, mapping_hits_df, reference_rules_file, classification_df):
        """Test that lr_defined is preferred over reference when scores are equal (rules)."""
        combined_class = merge_classifications(reference_rules_file, classification_df, "rules")
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "rules")
        
        # PB.2.1 maps to both REF5 (190, but Artifact) and LR.2.1 (190, Isoform)
        # After filtering artifacts, only LR.2.1 remains with lr_defined origin
        pb21 = result[result["rescue_candidate"] == "PB.2.1"]
        assert len(pb21) == 1
        assert pb21.iloc[0]["mapping_hit"] == "LR.2.1"
        assert pb21.iloc[0]["hit_origin"] == "lr_defined"
    
    def test_filter_ml_strategy_basic(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test filtering with ML strategy - basic checks."""
        combined_class = merge_classifications(reference_ml_file, classification_df[["isoform", "filter_result","POS_MLprob"]], "ml", thr=0.7)
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "ml")
        
        assert isinstance(result, pd.DataFrame)
        # Score should be calculated as alignment_score * hit_POS_MLprob
        assert "score" in result.columns
        assert "hit_POS_MLprob" in result.columns
        # Should only keep those classified as Isoform based on threshold
        assert all(result["hit_filter_result"] == "Isoform")
    
    def test_filter_ml_threshold_filtering(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test that ML strategy respects the probability threshold."""
        # With threshold 0.7, these should be kept (>= 0.7):
        # REF1 (0.95), REF2 (0.85), REF4 (0.75), REF6 (0.80), REF7 (0.72),
        # REF8 (0.88), LR.1.1 (0.90), LR.2.1 (0.92), REF10 (0.88), LR.6.1 (0.85)
        # These should be filtered out (< 0.7):
        # REF3 (0.45), REF5 (0.55), REF9 (0.60), LR.3.1 (0.50)
        
        combined_class = merge_classifications(reference_ml_file, classification_df[["isoform", "filter_result","POS_MLprob"]], "ml", thr=0.7)
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "ml")
        
        expected_kept = {"REF1", "REF2", "REF4", "REF6", "REF7", "REF8", 
                        "LR.1.1", "LR.2.1", "REF10", "LR.6.1"}
        expected_filtered = {"REF3", "REF5", "REF9", "LR.3.1"}
        
        result_hits = set(result["mapping_hit"].unique())
        
        # All results should be from expected kept transcripts
        assert result_hits.issubset(expected_kept)
        # No low-probability transcripts should be in results
        assert result_hits.isdisjoint(expected_filtered)
    
    def test_filter_ml_specific_candidates(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test specific rescue candidates with ML strategy."""
        combined_class = merge_classifications(reference_ml_file, classification_df[["isoform", "filter_result","POS_MLprob"]], "ml", thr=0.7)
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "ml")
        
        # PB.1.1 maps to REF1 (200, 0.95) and LR.1.1 (195, 0.90)
        # Scores: REF1 = 200*0.95 = 190, LR.1.1 = 195*0.90 = 175.5
        # Should select REF1 (higher weighted score)
        pb11 = result[result["rescue_candidate"] == "PB.1.1"]
        assert len(pb11) == 1
        assert pb11.iloc[0]["mapping_hit"] == "REF1"
        assert pb11.iloc[0]["score"] == pytest.approx(200 * 0.95)
        
        # PB.1.2 maps to REF2 (180, 0.85) and REF3 (180, 0.45)
        # REF3 filtered out (< 0.7), only REF2 remains
        pb12 = result[result["rescue_candidate"] == "PB.1.2"]
        assert len(pb12) == 1
        assert pb12.iloc[0]["mapping_hit"] == "REF2"
        assert pb12.iloc[0]["score"] == pytest.approx(180 * 0.85)
        
        # PB.2.1 maps to REF5 (190, 0.55) and LR.2.1 (190, 0.92)
        # REF5 filtered out (< 0.7), only LR.2.1 remains
        pb21 = result[result["rescue_candidate"] == "PB.2.1"]
        assert len(pb21) == 1
        assert pb21.iloc[0]["mapping_hit"] == "LR.2.1"
        assert pb21.iloc[0]["score"] == pytest.approx(190 * 0.92)
        
        # PB.3.2 maps to LR.3.1 (185, 0.50) and REF8 (170, 0.88)
        # LR.3.1 filtered out (< 0.7), only REF8 remains
        pb32 = result[result["rescue_candidate"] == "PB.3.2"]
        assert len(pb32) == 1
        assert pb32.iloc[0]["mapping_hit"] == "REF8"
        assert pb32.iloc[0]["score"] == pytest.approx(170 * 0.88)
        
        # PB.4.1 maps to REF9 (140, 0.60)
        # Filtered out (< 0.7)
        pb41 = result[result["rescue_candidate"] == "PB.4.1"]
        assert len(pb41) == 0
    
    def test_filter_ml_score_calculation(self, mapping_hits_df, reference_ml_file, classification_df):
        """Test that ML scores are correctly calculated as alignment_score * ML_prob."""
        combined_class = merge_classifications(reference_ml_file, classification_df[["isoform", "filter_result","POS_MLprob"]], "ml", thr=0.7)
        mapping_with_results = add_filter_results(mapping_hits_df, combined_class, classification_df)
        
        result = filter_mapping_hits(mapping_with_results, "ml")
        
        # Verify score calculation for all results
        for _, row in result.iterrows():
            expected_score = row["alignment_score"] * row["hit_POS_MLprob"]
            assert row["score"] == pytest.approx(expected_score)


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
        
    def test_merge_reintroduces_propperly(self, automatic_rescue_df):
        """ Test that a reference transcript is not reintroduced if it 
        was already reintroduced by automatic rescue."""
        full_rescue = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.2.1'],
            'mapping_hit': ['REF.AUTO.1', 'REF2'],
            'hit_origin': ['reference', 'reference']
        })
        
        inclusion_list = pd.Series(['REF.AUTO.1', 'REF2'])
        
        result = merge_rescue_modes(automatic_rescue_df, full_rescue, inclusion_list, "rules")
        
        # For the isoform PB.1.1 the reintroduced column should be 'no' 
        ref1_entries = result[result["artifact"] == "PB.1.1"]
        assert ref1_entries["reintroduced"].iloc[0] == "no"
        assert ref1_entries["rescue_mode"].iloc[0] == "rules_mapping"

    def test_merge_multiple_reintroduced(self, automatic_rescue_df):
        """ Test that multiple different reference transcripts are marked as reintroduced correctly."""
        full_rescue = pd.DataFrame({
            'rescue_candidate': ['PB.1.1', 'PB.2.1'],
            'mapping_hit': ['REF1', 'REF1'],
            'hit_origin': ['reference', 'reference']
        })
        
        inclusion_list = pd.Series(['REF1'])

        result = merge_rescue_modes(automatic_rescue_df, full_rescue, inclusion_list, "rules")
        # REF2 should be marked as reintroduced by mapping rescue
        ref_entries = result[result["assigned_transcript"] == "REF1"]
        assert ref_entries["reintroduced"].iloc[0] == "yes"
        assert ref_entries["reintroduced"].iloc[1] == "no"


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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
