"""Rules that rescue removes before filtering the reference transcriptome.

Requisites on long-read counts (FL, FL.<sample>, prevalence, prevalence_<group>)
and on subcategory cannot be evaluated on the reference, so rescue writes a copy
of the JSON file without them. The rules filter itself is not changed.
"""
import sys, os, json
import pytest
import pandas as pd
from unittest.mock import patch

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.utilities.filter.sqanti3_rules_filter import read_json_rules, apply_rules
from src.rescue_steps import is_ignored_rule, write_reference_rules


def isoform():
    return pd.Series({"isoform": "t1", "structural_category": "genic", "exons": 3})



@pytest.mark.parametrize("column", ["FL", "FL.K1", "prevalence", "prevalence_K", "prevalence_B_2", "subcategory"])
def test_ignored_columns(column):
    assert is_ignored_rule(column)


@pytest.mark.parametrize("column", ["FLAIR_score", "prevalent", "min_cov", "all_canonical", "perc_A_downstream_TTS"])
def test_kept_columns(column):
    assert not is_ignored_rule(column)


def test_write_reference_rules(tmp_path):
    """Only the ignored requisites are removed; categories and rules keep their structure."""
    rules = {
        "full-splice_match": [{"subcategory": "reference_match", "perc_A_downstream_TTS": [0, 59]}],
        "rest": [
            {"all_canonical": "canonical", "prevalence_K": 2},
            {"all_canonical": "canonical", "prevalence_B": 2},
            {"FL": 5},
        ],
    }
    src_json, out_json = tmp_path / "rules.json", tmp_path / "reference_rules.json"
    src_json.write_text(json.dumps(rules))

    removed = write_reference_rules(str(src_json), str(out_json))

    assert removed == ["FL", "prevalence_B", "prevalence_K", "subcategory"]
    assert json.loads(out_json.read_text()) == {
        "full-splice_match": [{"perc_A_downstream_TTS": [0, 59]}],
        "rest": [{"all_canonical": "canonical"}, {"all_canonical": "canonical"}, {}],
    }


def test_emptied_rule_accepts_reference(tmp_path):
    """A rule made only of ignored requisites accepts every reference transcript."""
    src_json, out_json = tmp_path / "rules.json", tmp_path / "reference_rules.json"
    src_json.write_text(json.dumps({"genic": [{"prevalence": 2}]}))
    write_reference_rules(str(src_json), str(out_json))
    rules = read_json_rules(str(out_json))
    assert apply_rules(isoform(), False, rules) == "Isoform"


def test_rescue_filters_reference_with_cleaned_rules(tmp_path):
    from src import rescue_steps
    src_json = tmp_path / "rules.json"
    src_json.write_text(json.dumps({"rest": [{"prevalence": 2, "min_cov": 3}]}))
    with patch.object(rescue_steps, "run_command") as mock_cmd, \
         patch.object(rescue_steps, "rescue_by_mapping", return_value=(None, None)):
        rescue_steps.run_rules_rescue("filt.txt", "ref.txt", None, None, None,
                                      str(tmp_path), str(src_json))
    ref_json = tmp_path / "reference_rules_filter" / "reference_rules.json"
    assert f"-j {ref_json} " in mock_cmd.call_args.args[0]
    assert json.loads(ref_json.read_text()) == {"rest": [{"min_cov": 3}]}
