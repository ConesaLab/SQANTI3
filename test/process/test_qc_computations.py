"""
Process-level tests for qc_computations functions.
Tests the various QC computation functions that modify isoforms_info dictionaries.
"""
import pytest
from src.qc_computations import (
    classify_fsm,
    full_length_quantification,
    isoforms_junctions,
    ratio_TSS_dict_reading,
    process_rts
)


def test_ratio_tss_basic_assignment(sample_isoforms_info):
    """Test basic ratio TSS assignment."""
    ratio_TSS_dict = {
        "PB.124830.1": {"return_ratio": 0.8},
        "PB.103714.1": {"return_ratio": 0.95},
        "PB.103724.1": {"return_ratio": 0.5}
    }
    
    result = ratio_TSS_dict_reading(sample_isoforms_info, ratio_TSS_dict)
    
    assert result["PB.124830.1"].ratio_TSS == 0.8
    assert result["PB.103714.1"].ratio_TSS == 0.95
    assert result["PB.103724.1"].ratio_TSS == 0.5


def test_ratio_tss_missing_isoforms(sample_isoforms_info):
    """Test that missing isoforms get ratio_TSS = 1."""
    ratio_TSS_dict = {
        "PB.124830.1": {"return_ratio": 0.8}
        # Other isoforms not in dict
    }
    
    result = ratio_TSS_dict_reading(sample_isoforms_info, ratio_TSS_dict)
    
    # Present isoform gets actual value
    assert result["PB.124830.1"].ratio_TSS == 0.8
    
    # Missing isoforms get default value of None
    assert result["PB.103714.1"].ratio_TSS is None
    assert result["PB.103724.1"].ratio_TSS is None
    assert result["PB.103781.1"].ratio_TSS is None


def test_ratio_tss_nan_handling(sample_isoforms_info):
    """Test that NaN values are converted to None."""
    ratio_TSS_dict = {
        "PB.124830.1": {"return_ratio": "nan"},
        "PB.103714.1": {"return_ratio": 0.95}
    }
    
    result = ratio_TSS_dict_reading(sample_isoforms_info, ratio_TSS_dict)
    
    # NaN should become None
    assert result["PB.124830.1"].ratio_TSS is None
    # Regular value should work
    assert result["PB.103714.1"].ratio_TSS == 0.95


def test_ratio_tss_extra_isoforms_in_file(sample_isoforms_info):
    """Test handling of isoforms in ratio file but not in isoforms_info."""
    ratio_TSS_dict = {
        "PB.124830.1": {"return_ratio": 0.8},
        "PB.EXTRA.1": {"return_ratio": 0.7}  # Not in isoforms_info
    }
    
    # Should not raise error, just log warning
    result = ratio_TSS_dict_reading(sample_isoforms_info, ratio_TSS_dict)
    
    assert result["PB.124830.1"].ratio_TSS == 0.8

## RTS is checked on its own test file

## classify_fsm tests

def test_classify_fsm_single_isoform_per_gene(sample_isoforms_info):
    """Test FSM classification when each gene has only one isoform."""
    # Modify sample data to have single isoforms per gene
    test_info = {
        "PB.124830.1": sample_isoforms_info["PB.124830.1"],
        "PB.103714.1": sample_isoforms_info["PB.103714.1"],
    }
    
    result = classify_fsm(test_info)
    
    # When gene has only one isoform, FSM_class should be "A"
    assert result["PB.124830.1"].FSM_class == "A"
    assert result["PB.103714.1"].FSM_class == "A"


def test_classify_fsm_gene_with_fsm(sample_isoforms_info):
    """Test FSM classification when gene has a full-splice_match isoform."""
    # Create multiple isoforms for same gene, one with FSM
    test_info = {
        "PB.1": sample_isoforms_info["PB.103714.1"],  # FSM
        "PB.2": sample_isoforms_info["PB.103724.1"],  # Not FSM
    }
    # Set them to same gene
    test_info["PB.1"].genes = ["GENE_SHARED"]
    test_info["PB.2"].genes = ["GENE_SHARED"]
    test_info["PB.1"].structural_category = "full-splice_match"
    test_info["PB.2"].structural_category = "novel_in_catalog"
    
    result = classify_fsm(test_info)
    
    # When gene has FSM, all isoforms get "C"
    assert result["PB.1"].FSM_class == "C"
    assert result["PB.2"].FSM_class == "C"


def test_classify_fsm_gene_without_fsm(sample_isoforms_info):
    """Test FSM classification when gene has no full-splice_match."""
    # Create multiple isoforms for same gene, none with FSM
    test_info = {
        "PB.1": sample_isoforms_info["PB.103724.1"],
        "PB.2": sample_isoforms_info["PB.103781.1"],
    }
    # Set them to same gene
    test_info["PB.1"].genes = ["GENE_SHARED"]
    test_info["PB.2"].genes = ["GENE_SHARED"]
    test_info["PB.1"].structural_category = "novel_in_catalog"
    test_info["PB.2"].structural_category = "novel_not_in_catalog"
    
    result = classify_fsm(test_info)
    
    # When gene has no FSM but multiple isoforms, get "B"
    assert result["PB.1"].FSM_class == "B"
    assert result["PB.2"].FSM_class == "B"


def test_fl_count_single_sample(sample_isoforms_info, fl_count_file_single, fields_class_cur):
    """Test FL count parsing for single-sample format."""
    result_info, result_fields = full_length_quantification(
        fl_count_file_single,
        sample_isoforms_info,
        fields_class_cur.copy()
    )
    
    # Check that FL counts were assigned
    assert result_info["PB.124830.1"].FL == 150
    assert result_info["PB.103714.1"].FL == 200
    assert result_info["PB.103724.1"].FL == 50
    assert result_info["PB.103781.1"].FL == 300
    assert result_info["PB.23068.2"].FL == 25
    
    # Fields should not change for single sample
    assert len(result_fields) == len(fields_class_cur)


def test_fl_count_multi_sample(sample_isoforms_info, fl_count_file_multi, fields_class_cur):
    """Test FL count parsing for multi-sample format."""
    result_info, result_fields = full_length_quantification(
        fl_count_file_multi,
        sample_isoforms_info,
        fields_class_cur.copy()
    )
    
    # Check that FL_dict was populated for multi-sample
    assert result_info["PB.124830.1"].FL_dict["sample1"] == 150
    assert result_info["PB.124830.1"].FL_dict["sample2"] == 200
    assert result_info["PB.124830.1"].FL_dict["sample3"] == 180
    
    assert result_info["PB.103714.1"].FL_dict["sample1"] == 200
    assert result_info["PB.103714.1"].FL_dict["sample2"] == 220
    assert result_info["PB.103714.1"].FL_dict["sample3"] == 190
    
    # Fields should be extended with sample names
    assert "FL.sample1" in result_fields
    assert "FL.sample2" in result_fields
    assert "FL.sample3" in result_fields


def test_fl_count_missing_isoforms(sample_isoforms_info, fl_count_file_single, fields_class_cur):
    """Test that missing isoforms get FL count of 0."""
    # Add an isoform not in the FL count file
    sample_isoforms_info["PB.MISSING.1"] = sample_isoforms_info["PB.124830.1"]
    sample_isoforms_info["PB.MISSING.1"].isoform = "PB.MISSING.1"
    
    result_info, _ = full_length_quantification(
        fl_count_file_single,
        sample_isoforms_info,
        fields_class_cur.copy()
    )
    
    # Missing isoform should get FL = 0
    assert result_info["PB.MISSING.1"].FL == 0


def test_fl_count_preserves_other_attributes(sample_isoforms_info, fl_count_file_single, fields_class_cur):
    """Test that FL quantification doesn't modify other isoform attributes."""
    original_category = sample_isoforms_info["PB.124830.1"].structural_category
    original_length = sample_isoforms_info["PB.124830.1"].length
    
    result_info, _ = full_length_quantification(
        fl_count_file_single,
        sample_isoforms_info,
        fields_class_cur.copy()
    )
    
    # Other attributes should remain unchanged
    assert result_info["PB.124830.1"].structural_category == original_category
    assert result_info["PB.124830.1"].length == original_length


@pytest.fixture
def mock_junction_reader():
    """Create a mock junction reader with test data."""
    # Simulate DictReader output
    junction_data = [
        {
            'isoform': 'PB.124830.1',
            'canonical': 'canonical',
            'bite_junction': 'FALSE',
            'indel_near_junct': 'FALSE',
            'sample_with_cov': '3',
            'total_coverage_unique': '100',
            'junction_number': '1'
        },
        {
            'isoform': 'PB.124830.1',
            'canonical': 'canonical',
            'bite_junction': 'FALSE',
            'indel_near_junct': 'FALSE',
            'sample_with_cov': '2',
            'total_coverage_unique': '80',
            'junction_number': '2'
        },
        {
            'isoform': 'PB.103724.1',
            'canonical': 'non_canonical',
            'bite_junction': 'TRUE',
            'indel_near_junct': 'TRUE',
            'sample_with_cov': '1',
            'total_coverage_unique': '50',
            'junction_number': '1'
        }
    ]
    return junction_data


def test_canonical_junction_detection(sample_isoforms_info, mock_junction_reader):
    """Test that canonical junctions are properly detected."""
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # PB.124830.1 has all canonical junctions
    assert result["PB.124830.1"].all_canonical == "canonical"
    
    # PB.103724.1 has non-canonical junction
    assert result["PB.103724.1"].all_canonical == "non_canonical"


def test_bite_junction_detection(sample_isoforms_info, mock_junction_reader):
    """Test that bite junctions are properly detected."""
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # PB.124830.1 has no bite junctions
    assert result["PB.124830.1"].bite == "FALSE"
    
    # PB.103724.1 has bite junction
    assert result["PB.103724.1"].bite == "TRUE"


def test_indel_near_junction_counting(sample_isoforms_info, mock_junction_reader):
    """Test that indels near junctions are counted."""
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # PB.124830.1 has no indels near junctions
    assert result["PB.124830.1"].n_indels_junc is None
    
    # PB.103724.1 has 1 indel near junction
    assert result["PB.103724.1"].n_indels_junc == 1


def test_min_coverage_calculation(sample_isoforms_info, mock_junction_reader):
    """Test minimum coverage calculations."""
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # PB.124830.1 should have min_cov = 80 (minimum of 100 and 80)
    assert result["PB.124830.1"].min_cov == 80
    assert result["PB.124830.1"].min_cov_pos == '2'  # Junction 2 has the minimum
    
    # PB.103724.1 should have min_cov = 50
    assert result["PB.103724.1"].min_cov == 50
    assert result["PB.103724.1"].min_cov_pos == '1'


def test_min_sample_coverage(sample_isoforms_info, mock_junction_reader):
    """Test minimum sample coverage across junctions."""
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # PB.124830.1 should have min_sample_cov = 2 (minimum of 3 and 2)
    assert result["PB.124830.1"].min_sample_cov == 2
    
    # PB.103724.1 should have min_sample_cov = 1
    assert result["PB.103724.1"].min_sample_cov == 1


def test_standard_deviation_calculation(sample_isoforms_info, mock_junction_reader):
    """Test standard deviation of junction coverage."""
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # PB.124830.1 has coverages [100, 80], should have sd calculated
    assert result["PB.124830.1"].sd is not None
    assert result["PB.124830.1"].sd > 0  # Should be ~10
    
    # PB.103724.1 has only one junction, sd should be 0
    assert result["PB.103724.1"].sd == 0


def test_isoforms_without_junctions(sample_isoforms_info, mock_junction_reader):
    """Test isoforms that have no junctions remain unmodified."""
    # PB.103714.1 and others are not in mock_junction_reader
    result = isoforms_junctions(sample_isoforms_info, mock_junction_reader)
    
    # These should remain None since they weren't in the junction file
    assert result["PB.103714.1"].all_canonical is None
    assert result["PB.103714.1"].bite is None
    assert result["PB.103714.1"].min_cov is None
