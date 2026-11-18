"""Prevalence columns in the rules filter.

The rules filter has no code specific to prevalence: prevalence and
prevalence_<group> are numeric columns of the classification, evaluated as
Min_Threshold like any other. These tests pin that contract.
"""
import sys, os, json
import pytest
import pandas as pd
from unittest.mock import mock_open, patch

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.utilities.filter.sqanti3_rules_filter import read_json_rules, apply_rules, get_reasons


# ---------------------------------------------------------------- helpers

def rules_from(json_data):
    """Parse a rules dict through read_json_rules() without touching disk."""
    with patch("builtins.open", mock_open(read_data=json.dumps(json_data))):
        return read_json_rules("dummy_path.json")


def isoform(**values):
    """Multi-exon genic isoform with the given prevalence values."""
    row = {"isoform": "t1", "structural_category": "genic", "exons": 3}
    row.update(values)
    return pd.Series(row)


def result(row, json_data):
    return apply_rules(row, False, rules_from(json_data))


# ---------------------------------------------------------------- filter

def test_prevalence_is_a_min_threshold():
    """A prevalence requisite is parsed like any other numeric requisite."""
    rules = rules_from({"genic": [{"prevalence": 2}]})["genic"][0]
    assert rules.to_dict("records") == [{"column": "prevalence", "type": "Min_Threshold", "rule": 2}]


@pytest.mark.parametrize("prevalence, expected", [(3, "Isoform"), (2, "Isoform"), (1, "Artifact"), (0, "Artifact")])
def test_global_prevalence(prevalence, expected):
    assert result(isoform(prevalence=prevalence), {"genic": [{"prevalence": 2}]}) == expected


def test_prevalence_one_still_filters():
    """prevalence 0 is possible with fractional EM counts, so a threshold of 1 is not a no-op."""
    assert result(isoform(prevalence=0), {"genic": [{"prevalence": 1}]}) == "Artifact"


def test_na_prevalence_fails():
    """An isoform missing from --fl_count has NA prevalence: it was detected in no sample."""
    row = isoform(prevalence=float("nan"))
    assert result(row, {"genic": [{"prevalence": 1}]}) == "Artifact"
    reasons = get_reasons(row, False, rules_from({"genic": [{"prevalence": 1}]}))
    assert reasons["filter_reason"] == "NA value in prevalence"


AND_RULES = {"genic": [{"prevalence_K": 2, "prevalence_B": 2}]}
OR_RULES = {"genic": [{"prevalence_K": 2}, {"prevalence_B": 2}]}


@pytest.mark.parametrize("k, b, and_expected, or_expected", [
    (2, 2, "Isoform", "Isoform"),
    (2, 0, "Artifact", "Isoform"),
    (0, 3, "Artifact", "Isoform"),
    (1, 1, "Artifact", "Artifact"),
])
def test_groups_and_or(k, b, and_expected, or_expected):
    """Requisites of one rule are AND; rules of one structural category are OR."""
    row = isoform(prevalence_K=k, prevalence_B=b)
    assert result(row, AND_RULES) == and_expected
    assert result(row, OR_RULES) == or_expected


def test_global_or_groups():
    """The global prevalence can be combined with the per-group ones."""
    rules = {"genic": [{"prevalence": 4}, {"prevalence_K": 2, "prevalence_B": 2}]}
    assert result(isoform(prevalence=4, prevalence_K=1, prevalence_B=1), rules) == "Isoform"
    assert result(isoform(prevalence=3, prevalence_K=2, prevalence_B=1), rules) == "Artifact"
    assert result(isoform(prevalence=3, prevalence_K=2, prevalence_B=2), rules) == "Isoform"


def test_group_failure_reason():
    row = isoform(prevalence_K=1, prevalence_B=0)
    reasons = get_reasons(row, False, rules_from(AND_RULES))["filter_reason"].split("; ")
    assert sorted(reasons) == ["prevalence_B: 0 < 2", "prevalence_K: 1 < 2"]


def test_missing_group_column_exits():
    """A group column absent from the classification is reported like any misspelled column."""
    with pytest.raises(SystemExit):
        result(isoform(prevalence_K=2), {"genic": [{"prevalence_X": 2}]})
