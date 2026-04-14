import pandas as pd
from src.utilities.filter.sqanti3_rules_filter import apply_rules, get_reasons, junction_related_columns

class MockRulesList:
    def __init__(self, elements):
        self.elements = elements
        self.keys_list = list(elements.keys())
    def __iter__(self):
        return iter(self.keys_list)
    def __getitem__(self, idx):
        return self.elements[self.keys_list[idx]]

def test_monoexon_junction_failsafe_ignored():
    isoform_info = {
        "isoform": "PB.1.1",
        "exons": 1,
        "structural_category": "rest",
        "min_cov": 0,
        "bite": "False"
    }

    rules_df = pd.DataFrame({
        "column": ["min_cov", "bite"],
        "type": ["Min_Threshold", "Category"],
        "rule": [10, "True"]
    })

    json_df = {"rest": True}
    rules_list = MockRulesList({"rest": rules_df})

    # Since it's a monoexon and both rules are junction related, it should ignore them and return Isoform
    assert apply_rules(isoform_info, False, json_df, rules_list) == "Isoform"
    
    reasons_df = get_reasons(isoform_info, False, json_df, rules_list)
    assert len(reasons_df) == 0

def test_multiexon_junction_failsafe_enforced():
    isoform_info = {
        "isoform": "PB.1.2",
        "exons": 2,
        "structural_category": "rest",
        "min_cov": 0,
        "bite": "False"
    }

    rules_df = pd.DataFrame({
        "column": ["min_cov", "bite"],
        "type": ["Min_Threshold", "Category"],
        "rule": [10, "true"]
    })

    json_df = {"rest": True}
    rules_list = MockRulesList({"rest": rules_df})

    # Since it's multiexon, rules are enforced -> fails min_cov and bite -> Artifact
    assert apply_rules(isoform_info, False, json_df, rules_list) == "Artifact"
    
    reasons_df = get_reasons(isoform_info, False, json_df, rules_list)
    assert len(reasons_df) == 2
    assert set(reasons_df["reasons"].tolist()) == {"Low min_cov", "Out bite"}

