import sys,os,pytest
import json
import pandas as pd
from unittest.mock import mock_open, patch
from typing import Any

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)
from src.utilities.filter.sqanti3_rules_filter import read_json_rules


## Reading JSON rules
@pytest.fixture
def json_default():    return os.path.join(main_path, "test","test_data","filter_rules.json")


def test_read_json_rules(json_default: str):
    res_list = read_json_rules(json_default)
    assert len(res_list) == 2
    assert list(res_list.keys()) == ['full-splice_match','rest']
    assert len(res_list['full-splice_match']) == 1
    assert len(res_list['rest']) == 3


@pytest.fixture
def sample_json_data_good():
    return {
        "full-splice_match": [
            {"col1": [1, 5], "col2": "value"},
            {"col3": 10}
        ],
        "genic": [
            {"col4": ["A", "B", "C"]}
        ]
    }

@pytest.fixture
def sample_json_data_bad():
    return {
        "category1": [
            {"col1": [1, 5], "col2": "value"},
            {"col3": 10}
        ],
        "category2": [
            {"col4": ["A", "B", "C"]}
        ]
    }

def test_read_json_rules_structure(sample_json_data_good):
    with patch("builtins.open", mock_open(read_data=json.dumps(sample_json_data_good))):
        result = read_json_rules("dummy_path.json")
    
    assert isinstance(result, dict)
    assert set(result.keys()) == {"full-splice_match", "genic","rest"}
    assert all(isinstance(v, list) for v in result.values())
    assert all(isinstance(df, pd.DataFrame) for category in result.values() for df in category)

def test_read_json_rules_invalid_structure(sample_json_data_bad):
    with patch("builtins.open", mock_open(read_data=json.dumps(sample_json_data_bad))):
        with patch("sys.exit") as mock_exit:
            read_json_rules("dummy_path.json")
            mock_exit.assert_called_once_with(1)
    

def test_read_json_rules_dataframe_columns(sample_json_data_good):
    with patch("builtins.open", mock_open(read_data=json.dumps(sample_json_data_good))):
        result = read_json_rules("dummy_path.json")
    
    expected_columns = ['column', 'type', 'rule']
    for category in result.values():
        for df in category:
            assert list(df.columns) == expected_columns

def test_read_json_rules_numeric_range():
    json_data = {"genic": [{"col": [1, 5]}]}
    with patch("builtins.open", mock_open(read_data=json.dumps(json_data))):
        result = read_json_rules("dummy_path.json")
    
    df = result["genic"][0]
    assert df.shape == (2, 3)
    assert df.iloc[0].tolist() == ["col", "Min_Threshold", 1]
    assert df.iloc[1].tolist() == ["col", "Max_Threshold", 5]

def test_read_json_rules_category_list():
    json_data = {"genic": [{"col": ["A", "B", "C"]}]}
    with patch("builtins.open", mock_open(read_data=json.dumps(json_data))):
        result = read_json_rules("dummy_path.json")
    
    df = result["genic"][0]
    assert df.shape == (1, 3)
    assert all(df["type"] == "Category")
    assert df.iloc[0]["rule"] == ["a", "b", "c"]

def test_read_json_rules_single_numeric():
    json_data = {"genic": [{"col": 10}]}
    with patch("builtins.open", mock_open(read_data=json.dumps(json_data))):
        result = read_json_rules("dummy_path.json")
    
    df = result["genic"][0]
    assert df.shape == (1, 3)
    assert df.iloc[0].tolist() == ["col", "Min_Threshold", 10]

def test_read_json_rules_single_string():
    json_data = {"genic": [{"col": "value"}]}
    with patch("builtins.open", mock_open(read_data=json.dumps(json_data))):
        result = read_json_rules("dummy_path.json")
    
    df = result["genic"][0]
    assert df.shape == (1, 3)
    assert df.iloc[0].tolist() == ["col", "Category", "value"]

## Applying rules
from src.utilities.filter.sqanti3_rules_filter import apply_rules

@pytest.fixture
def classification_df():
    file = os.path.join(main_path, "test", "test_data", "test_isoforms_classification.tsv")
    return pd.read_csv(file, sep="\t")

@pytest.fixture
def rules_dict(json_default: str):
    return read_json_rules(json_default)

def test_apply_rules_isoform_rest(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[1]
    result = apply_rules(row, False, rules_dict)
    assert result == "Isoform"
def test_apply_rules_isoform_monoexon(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[0]
    result = apply_rules(row, False, rules_dict)
    assert result == "Isoform"

def test_apply_rules_isoform_fsm(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[16]
    result = apply_rules(row, False, rules_dict)
    assert result == "Isoform"
    
def test_apply_rules_isoform_fsm_monoexon(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[15]
    result = apply_rules(row, False, rules_dict)
    assert result == "Isoform"

def test_apply_rules_artifact_rest(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[21]
    print(row['isoform'])
    result = apply_rules(row, False, rules_dict)
    assert result == "Artifact"

def test_apply_rules_artifact_fsm(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[39]
    result = apply_rules(row, False, rules_dict)
    assert result == "Artifact"

def test_apply_rules_force_monoexon(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    #REST
    row = classification_df.iloc[0]
    result = apply_rules(row, False, rules_dict)
    assert result == "Isoform"
    #FSM
    row = classification_df.iloc[15]
    result = apply_rules(row, True, rules_dict)
    assert result == "Artifact"
    
# Reasons for filtering
from src.utilities.filter.sqanti3_rules_filter import get_reasons

def test_get_reasons_isoform_passes(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[1]  # Should be an Isoform
    result = get_reasons(row, False, rules_dict)
    
    assert result["isoform"] == row["isoform"]
    assert result["structural_category"] == row["structural_category"]
    assert result["filter_reason"] == ""

def test_get_reasons_artifact_due_to_monoexon(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[0]  # 1-exon FSM or rest
    result = get_reasons(row, True, rules_dict)
    
    assert "Mono-exonic" in result["filter_reason"]
    assert result["isoform"] == row["isoform"]

def test_get_reasons_threshold_violation(classification_df: pd.DataFrame, rules_dict: dict[str, Any]):
    row = classification_df.iloc[39]  # Should be an Artifact by threshold
    result = get_reasons(row, False, rules_dict)

    assert result["isoform"] == row["isoform"]
    assert result["filter_reason"] != ""  # At least one reason should be present
    assert any(x in result["filter_reason"] for x in ["<", ">"])  # Threshold message present

def test_get_reasons_categorical_mismatch():
    row = pd.Series({
        "isoform": "fake_iso",
        "structural_category": "category",
        "col": "Mismatch"
    })
    rules_dict = {
        "category": [pd.DataFrame([["col", "Category", "expected"]], columns=["column", "type", "rule"])],
        "rest": []
    }
    result = get_reasons(row, False, rules_dict)
    
    assert "col: Mismatch" in result["filter_reason"]

def test_get_reasons_multiple_failures():
    row = pd.Series({
        "isoform": "fake_iso",
        "structural_category": "category",
        "col1": 3,      # should fail min
        "col2": 20,     # should fail max
        "col3": "bad"   # should fail category
    })
    rules_df = pd.DataFrame([
        ["col1", "Min_Threshold", 5],
        ["col2", "Max_Threshold", 10],
        ["col3", "Category", "good"]
    ], columns=["column", "type", "rule"])
    rules_dict = {
        "category": [rules_df],
        "rest": []
    }
    result = get_reasons(row, False, rules_dict)

    assert "col1: 3 < 5" in result["filter_reason"]
    assert "col2: 20 > 10" in result["filter_reason"]
    assert "col3: bad" in result["filter_reason"]
    assert result["isoform"] == "fake_iso"