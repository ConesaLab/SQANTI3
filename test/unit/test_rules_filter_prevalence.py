import sys, os, json, logging
import pytest
import pandas as pd
from unittest.mock import mock_open, patch

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.module_logging import filter_logger
from src.filter_argparse import filter_argparse
from src.utilities.filter.sqanti3_rules_filter import (
    read_json_rules, apply_rules, get_reasons, rules_filter,
    get_highest_min_prevalence, drop_prevalence_rules, check_prevalence_column
)


# ---------------------------------------------------------------- fixtures

@pytest.fixture
def filter_caplog(caplog, monkeypatch):
    """caplog does not see filter_logger because it has propagate=False."""
    monkeypatch.setattr(filter_logger, "propagate", True)
    caplog.set_level(logging.INFO)
    return caplog


def rules_from(json_data):
    """Parse a rules dict through read_json_rules() without touching disk."""
    with patch("builtins.open", mock_open(read_data=json.dumps(json_data))):
        return read_json_rules("dummy_path.json")


def make_classif(prevalence="multisample"):
    """Toy classification with three genic isoforms.

    t1 passes everything, t2 fails only prevalence, t3 fails prevalence and min_cov.
    prevalence:
        "multisample" -> prevalence column plus FL.s1 / FL.s2
        "na"          -> prevalence column all NA, no FL.* columns (reference or single-sample QC)
        "missing"     -> no prevalence column (QC from a version without it)
    """
    df = pd.DataFrame({
        "isoform": ["t1", "t2", "t3"],
        "structural_category": ["genic"] * 3,
        "exons": [2, 2, 2],
        "min_cov": [5, 5, 1],
    })
    if prevalence == "multisample":
        df["FL.s1"] = [3, 4, 0]
        df["FL.s2"] = [2, 0, 0]
        df["prevalence"] = [2, 1, 0]
    elif prevalence == "na":
        df["prevalence"] = pd.NA
    return df


def write_inputs(tmp_path, classif, json_data):
    class_file = tmp_path / "classification.txt"
    json_file = tmp_path / "rules.json"
    classif.to_csv(class_file, sep="\t", index=False)
    json_file.write_text(json.dumps(json_data))
    return str(class_file), str(json_file)


def read_pass(prefix):
    with open(f"{prefix}_pass_isoforms.txt") as f:
        return [line.strip() for line in f]


PREV_RULES = {"genic": [{"min_prevalence": 2, "min_cov": 3}], "rest": []}


# ---------------------------------------------------------------- read_json_rules

@pytest.mark.parametrize("key", ["min_prevalence", "prevalence", "Min_Prevalence"])
def test_read_json_rules_prevalence_keys(key):
    rules = rules_from({"genic": [{key: 2}]})
    assert rules["genic"][0].iloc[0].tolist() == ["prevalence", "min_prevalence", 2]


@pytest.mark.parametrize("value", [0, -1, 1.5, "2", True, [1, 2]])
def test_read_json_rules_invalid_prevalence(value):
    with pytest.raises(SystemExit):
        rules_from({"genic": [{"min_prevalence": value}]})


# ---------------------------------------------------------------- get_highest_min_prevalence

def test_highest_min_prevalence_none_without_rules():
    assert get_highest_min_prevalence(rules_from({"genic": [{"min_cov": 3}]})) is None


def test_highest_min_prevalence_across_branches():
    rules = rules_from({
        "genic": [{"min_prevalence": 2}, {"min_prevalence": 4, "min_cov": 3}],
        "full-splice_match": [{"min_prevalence": 3}],
    })
    assert get_highest_min_prevalence(rules) == 4


# ---------------------------------------------------------------- drop_prevalence_rules

def test_drop_prevalence_rules_keeps_other_rules():
    rules = rules_from({"genic": [{"min_prevalence": 2, "min_cov": 3}]})
    dropped = drop_prevalence_rules(rules)
    assert dropped["genic"][0]["type"].tolist() == ["Min_Threshold"]
    assert get_highest_min_prevalence(dropped) is None


def test_drop_prevalence_rules_does_not_modify_input():
    rules = rules_from({"genic": [{"min_prevalence": 2, "min_cov": 3}]})
    drop_prevalence_rules(rules)
    assert "min_prevalence" in rules["genic"][0]["type"].tolist()


def test_drop_prevalence_rules_empty_branch_accepts_all():
    rules = drop_prevalence_rules(rules_from({"genic": [{"min_prevalence": 2}]}))
    assert rules["genic"][0].empty
    row = make_classif("missing").iloc[2]
    assert apply_rules(row, False, rules) == "Isoform"


# ---------------------------------------------------------------- apply_rules / get_reasons

@pytest.mark.parametrize("idx, expected", [(0, "Isoform"), (1, "Artifact"), (2, "Artifact")])
def test_apply_rules_min_prevalence(idx, expected):
    rules = rules_from(PREV_RULES)
    assert apply_rules(make_classif().iloc[idx], False, rules) == expected


def test_apply_rules_min_prevalence_na_is_artifact():
    """An isoform missing from the --fl_count file has NA prevalence:
    it was detected in no sample, so it fails the requisite."""
    rules = rules_from({"genic": [{"min_prevalence": 1}]})
    row = make_classif().iloc[0].copy()
    row["prevalence"] = pd.NA
    assert apply_rules(row, False, rules) == "Artifact"


def test_get_reasons_min_prevalence_na():
    rules = rules_from({"genic": [{"min_prevalence": 1}]})
    row = make_classif().iloc[0].copy()
    row["prevalence"] = pd.NA
    assert get_reasons(row, False, rules)["filter_reason"] == "NA value in prevalence"


def test_get_reasons_min_prevalence():
    rules = rules_from(PREV_RULES)
    classif = make_classif()
    assert get_reasons(classif.iloc[1], False, rules)["filter_reason"] == \
        "prevalence: 1 < 2 Multisample-artifact"
    reasons = set(get_reasons(classif.iloc[2], False, rules)["filter_reason"].split("; "))
    assert reasons == {"prevalence: 0 < 2 Multisample-artifact", "min_cov: 1 < 3"}


# ---------------------------------------------------------------- check_prevalence_column

@pytest.mark.parametrize("prevalence", ["missing", "na"])
def test_check_prevalence_column_no_values(prevalence, filter_caplog):
    with pytest.raises(SystemExit):
        check_prevalence_column(make_classif(prevalence), 2)
    errors = [r for r in filter_caplog.records if r.levelno == logging.ERROR]
    # One message pair only: the sample-count error must not appear as well
    assert len(errors) == 2
    assert "has no prevalence values" in errors[0].getMessage()


def test_check_prevalence_column_threshold_above_samples(filter_caplog):
    with pytest.raises(SystemExit):
        check_prevalence_column(make_classif(), 3)
    assert "min_prevalence 3 is larger than the number of samples" in filter_caplog.text
    assert "(2)" in filter_caplog.text


@pytest.mark.parametrize("threshold", [1, 2])
def test_check_prevalence_column_valid(threshold):
    check_prevalence_column(make_classif(), threshold)  # must not exit


# ---------------------------------------------------------------- rules_filter scenarios

@pytest.mark.parametrize("prevalence, ignore_prevalence, rules, expected", [
    # 1. Old QC (no column), user filter         -> error
    ("missing", False, PREV_RULES, SystemExit),
    # 2. Old QC used as reference by rescue       -> prevalence ignored
    ("missing", True,  PREV_RULES, ["t1", "t2"]),
    # 3. Current reference (all NA) by rescue      -> prevalence ignored
    ("na",      True,  PREV_RULES, ["t1", "t2"]),
    # 4. Current QC without multisample, user     -> error
    ("na",      False, PREV_RULES, SystemExit),
    # 5. Multisample, threshold 3 with 2 samples  -> error
    ("multisample", False, {"genic": [{"min_prevalence": 3}], "rest": []}, SystemExit),
    # Normal multisample case
    ("multisample", False, PREV_RULES, ["t1"]),
    # Multisample with the flag on behaves as if the rule were not there
    ("multisample", True,  PREV_RULES, ["t1", "t2"]),
])
def test_rules_filter_scenarios(tmp_path, prevalence, ignore_prevalence, rules, expected):
    class_file, json_file = write_inputs(tmp_path, make_classif(prevalence), rules)
    prefix = str(tmp_path / "out")
    run = lambda: rules_filter(class_file, json_file, False, prefix, filter_logger,
                               ignore_prevalence=ignore_prevalence)
    if expected is SystemExit:
        with pytest.raises(SystemExit):
            run()
        assert not os.path.exists(f"{prefix}_pass_isoforms.txt")
    else:
        run()
        assert read_pass(prefix) == expected


def test_rules_filter_isoform_missing_from_counts(tmp_path):
    """Only t1 lacks prevalence (absent from --fl_count); the rest have values."""
    classif = make_classif()
    classif["prevalence"] = [pd.NA, 1, 0]
    rules = {"genic": [{"min_prevalence": 1, "min_cov": 3}], "rest": []}
    class_file, json_file = write_inputs(tmp_path, classif, rules)
    prefix = str(tmp_path / "out")
    rules_filter(class_file, json_file, False, prefix, filter_logger)
    assert read_pass(prefix) == ["t2"]
    

def test_rules_filter_ignore_is_logged(tmp_path, filter_caplog):
    class_file, json_file = write_inputs(tmp_path, make_classif("na"), PREV_RULES)
    rules_filter(class_file, json_file, False, str(tmp_path / "out"), filter_logger,
                 ignore_prevalence=True)
    assert "Ignoring min_prevalence rules (--ignore_prevalence)" in filter_caplog.text


def test_rules_filter_reasons_file(tmp_path):
    class_file, json_file = write_inputs(tmp_path, make_classif(), PREV_RULES)
    prefix = str(tmp_path / "out")
    rules_filter(class_file, json_file, False, prefix, filter_logger)
    reasons = pd.read_csv(f"{prefix}_filtering_reasons.txt", sep="\t").set_index("isoform")
    assert reasons.loc["t2", "filter_reason"] == "prevalence: 1 < 2 Multisample-artifact"
    assert set(reasons.loc["t3", "filter_reason"].split("; ")) == \
        {"prevalence: 0 < 2 Multisample-artifact", "min_cov: 1 < 3"}


def test_rules_filter_without_prevalence_rules_unchanged(tmp_path):
    """6. Rules with no min_prevalence on a classification without the column:
    the flag must make no difference and nothing must fail."""
    class_file = os.path.join(main_path, "test", "test_data", "isoforms", "test_isoforms_classification.tsv")
    json_file = os.path.join(main_path, "test", "test_data", "other", "filter_rules.json")
    assert "prevalence" not in pd.read_csv(class_file, sep="\t", nrows=0).columns

    outputs = {}
    for flag in (True, False):
        prefix = str(tmp_path / f"out_{flag}")
        rules_filter(class_file, json_file, False, prefix, filter_logger, ignore_prevalence=flag)
        outputs[flag] = (read_pass(prefix),
                         pd.read_csv(f"{prefix}_filtering_reasons.txt", sep="\t"))
    assert outputs[True][0] == outputs[False][0]
    pd.testing.assert_frame_equal(outputs[True][1], outputs[False][1])


# ---------------------------------------------------------------- argparse and callers

@pytest.mark.parametrize("extra, expected", [([], False), (["--ignore_prevalence"], True)])
def test_argparse_ignore_prevalence(extra, expected):
    args = filter_argparse().parse_args(["rules", "--sqanti_class", "x"] + extra)
    assert args.ignore_prevalence is expected


@pytest.mark.parametrize("flag", [True, False])
def test_run_rules_passes_flag(tmp_path, flag):
    from src import filter_steps
    args = filter_argparse().parse_args(
        ["rules", "--sqanti_class", "x", "-d", str(tmp_path), "-o", "out", "--skip_report"]
        + (["--ignore_prevalence"] if flag else []))
    (tmp_path / "out_pass_isoforms.txt").write_text("")
    with patch.object(filter_steps, "rules_filter") as mock_rf:
        filter_steps.run_rules(args)
    assert mock_rf.call_args.kwargs["ignore_prevalence"] is flag


@pytest.mark.parametrize("value, expected", [(False, ""), (True, "--ignore_prevalence")])
def test_wrapper_serializes_ignore_prevalence(value, expected):
    """store_true: false in the YAML must be omitted, true must emit the bare flag."""
    from src.wrapper_utils import format_options
    assert format_options({"ignore_prevalence": value}) == expected


def test_rescue_disables_prevalence_on_reference(tmp_path):
    from src import rescue_steps
    with patch.object(rescue_steps, "run_command") as mock_cmd, \
         patch.object(rescue_steps, "rescue_by_mapping", return_value=(None, None)):
        rescue_steps.run_rules_rescue("filt.txt", "ref.txt", None, None, None,
                                      str(tmp_path), "rules.json")
    assert mock_cmd.call_args.args[0].endswith("--ignore_prevalence")