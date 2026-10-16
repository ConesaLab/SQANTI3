"""Tests for the YAML-to-CLI translation in the SQANTI3 wrapper."""
from src.wrapper_utils import format_options


class TestFormatOptionsBooleans:

    def test_true_emits_bare_flag(self):
        assert format_options({"requant": True}) == "--requant"

    def test_false_negatable_option_emits_no_flag(self):
        """requant: false must actively disable a now-default option."""
        assert format_options({"requant": False}) == "--no-requant"

    def test_yaml_string_booleans_are_normalised(self):
        assert format_options({"requant": "true"}) == "--requant"
        assert format_options({"requant": "false"}) == "--no-requant"

    def test_false_store_true_option_is_dropped(self):
        """Non-negatable flags have no --no- form and must be omitted."""
        assert format_options({"skipORF": False}) == ""

    def test_empty_and_none_values_are_dropped(self):
        assert format_options({"refClassif": "", "counts": None}) == ""

    def test_value_options_are_passed_through(self):
        assert format_options({"counts": "counts.tsv"}) == "--counts counts.tsv"

    def test_zero_is_kept_as_a_value(self):
        assert format_options({"threshold": 0}) == "--threshold 0"
