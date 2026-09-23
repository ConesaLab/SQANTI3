import sys, os, json, logging
import pytest
import pandas as pd
from unittest.mock import mock_open, patch

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.module_logging import filter_logger
from src.filter_argparse import filter_argparse
from src.utilities.filter import sqanti3_rules_filter
from src.utilities.filter.sqanti3_rules_filter import (
    read_json_rules, apply_rules, get_reasons, rules_filter,
    read_counts_design, check_counts_design, resolve_prevalence_rules,
    add_group_prevalence, passes_prevalence, has_prevalence_rules
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


def design_from(text):
    """Parse a design file content through read_counts_design() without touching disk."""
    with patch("builtins.open", mock_open(read_data=text)):
        return read_counts_design("dummy_design.json")


DESIGN = {"K": ["K1", "K2", "K3"], "B": ["B1", "B2"]}


def make_classif(fl=True):
    """Toy classification with five genic isoforms, two groups and one extra sample.

    Per-group prevalence (detection = at least one full read):

        isoform  K  B  all samples  min_cov
        t1       2  0  2            5
        t2       0  2  2            5
        t3       1  0  2 (MIX)      5    0.4 and 0.9 in K are not detections
        t4       3  0  4 (MIX)      1    fails min_cov
        t5       1  1  2            5    detected in both groups, enough in neither

    MIX is not in DESIGN: it stands for a sample that must not take part in
    the filter. With fl=False there are no FL.* columns and prevalence is NA.
    """
    df = pd.DataFrame({
        "isoform": ["t1", "t2", "t3", "t4", "t5"],
        "structural_category": ["genic"] * 5,
        "exons": [2] * 5,
        "min_cov": [5, 5, 5, 1, 5],
    })
    if fl:
        df["FL.K1"] = [3, 0, 0.4, 5, 1]
        df["FL.K2"] = [2, 0, 0.9, 5, 0]
        df["FL.K3"] = [0, 0, 1.0, 5, 0]
        df["FL.B1"] = [0, 4, 0, 0, 1]
        df["FL.B2"] = [0, 3, 0, 0, 0]
        df["FL.MIX"] = [0, 0, 9, 9, 0]
        df["prevalence"] = (df.filter(like="FL.") >= 1).sum(axis=1)
    else:
        df["prevalence"] = pd.NA
    return df


def write_inputs(tmp_path, classif, rules, design=None):
    class_file = tmp_path / "classification.txt"
    json_file = tmp_path / "rules.json"
    classif.to_csv(class_file, sep="\t", index=False)
    json_file.write_text(json.dumps(rules))
    design_file = None
    if design is not None:
        design_file = tmp_path / "design.json"
        design_file.write_text(json.dumps(design))
        design_file = str(design_file)
    return str(class_file), str(json_file), design_file


def read_pass(prefix):
    with open(f"{prefix}_pass_isoforms.txt") as f:
        return [line.strip() for line in f]


INT_RULES = {"genic": [{"min_prevalence": 2, "min_cov": 3}], "rest": []}
GROUP_RULES = {"genic": [{"min_prevalence": {"K": 3, "B": 1}, "min_cov": 3}], "rest": []}


# ---------------------------------------------------------------- check_min_prevalence (via read_json_rules)

def test_read_json_rules_keeps_group_thresholds():
    rules = rules_from({"genic": [{"min_prevalence": {"K": 2, "B": 3}}]})
    assert rules["genic"][0].iloc[0].tolist() == ["prevalence", "min_prevalence", {"K": 2, "B": 3}]


@pytest.mark.parametrize("value", [
    {}, {"K": 0}, {"K": -1}, {"K": 1.5}, {"K": "2"}, {"K": True}, {"K": [2]}, {"K": {"x": 2}},
    {"K": 2, "B": 0},
])
def test_read_json_rules_invalid_group_thresholds(value):
    with pytest.raises(SystemExit):
        rules_from({"genic": [{"min_prevalence": value}]})


def test_invalid_group_threshold_names_the_group(filter_caplog):
    with pytest.raises(SystemExit):
        rules_from({"genic": [{"min_prevalence": {"K": 2, "B": 0}}]})
    assert "genic (group 'B')" in filter_caplog.text


# ---------------------------------------------------------------- read_counts_design

def test_read_counts_design_valid():
    design = design_from(json.dumps(DESIGN))
    assert design == DESIGN
    assert list(design) == ["K", "B"]  # file order is kept


@pytest.mark.parametrize("text", [
    "not json",
    json.dumps(["K1", "K2"]),                       # not an object
    json.dumps({}),                                 # no groups
    json.dumps({"": ["K1"]}),                       # empty group name
    json.dumps({" ": ["K1"]}),                      # blank group name
    json.dumps({"K": []}),                          # empty group
    json.dumps({"K": "K1"}),                        # samples not a list
    json.dumps({"K": ["K1", 2]}),                   # non-string sample
    json.dumps({"K": ["K1", ""]}),                  # empty sample name
    json.dumps({"K": ["K1", "K1"]}),                # sample twice in a group
    json.dumps({"K": ["K1"], "B": ["K1"]}),         # sample in two groups
    '{"K": ["K1"], "K": ["K2"]}',                   # duplicated group
])
def test_read_counts_design_invalid(text):
    with pytest.raises(SystemExit):
        design_from(text)


def test_read_counts_design_duplicated_group_is_reported(filter_caplog):
    """json.load would silently keep the last "K"; the design must be rejected instead."""
    with pytest.raises(SystemExit):
        design_from('{"K": ["K1"], "B": ["B1"], "K": ["K2"]}')
    assert "Duplicated group names" in filter_caplog.text
    assert "['K']" in filter_caplog.text


def test_read_counts_design_sample_in_two_groups_is_reported(filter_caplog):
    with pytest.raises(SystemExit):
        design_from(json.dumps({"K": ["K1", "S"], "B": ["S"]}))
    assert "'S'" in filter_caplog.text
    assert "'K' and 'B'" in filter_caplog.text


# ---------------------------------------------------------------- check_counts_design

def test_check_counts_design_no_fl_columns(filter_caplog):
    with pytest.raises(SystemExit):
        check_counts_design(make_classif(fl=False), DESIGN)
    assert "no per-sample counts" in filter_caplog.text


def test_check_counts_design_missing_sample(filter_caplog):
    with pytest.raises(SystemExit):
        check_counts_design(make_classif(), {"K": ["K1", "K9"]})
    assert "['K9']" in filter_caplog.text
    assert "Available samples" in filter_caplog.text


def test_check_counts_design_expects_names_without_prefix():
    """Samples are named as in the --fl_count header, not as classification columns."""
    with pytest.raises(SystemExit):
        check_counts_design(make_classif(), {"K": ["FL.K1", "FL.K2"]})


def test_check_counts_design_unused_samples_warn(filter_caplog):
    check_counts_design(make_classif(), DESIGN)  # must not exit
    warnings = [r.getMessage() for r in filter_caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "['MIX']" in warnings[0]


def test_check_counts_design_all_samples_used_no_warning(filter_caplog):
    check_counts_design(make_classif(), {**DESIGN, "M": ["MIX"]})
    assert "not assigned to any group" not in filter_caplog.text


def test_check_counts_design_single_sample_group_warns(filter_caplog):
    check_counts_design(make_classif(), {"K": ["K1", "K2", "K3"], "B": ["B1"]})
    assert "Group 'B' has a single sample" in filter_caplog.text


# ---------------------------------------------------------------- resolve_prevalence_rules

def prevalence_values(rules_dict, sc="genic"):
    return [rules.loc[rules["type"] == "min_prevalence", "rule"].tolist()
            for rules in rules_dict[sc]]


def test_resolve_without_design_keeps_integers():
    rules = resolve_prevalence_rules(rules_from(INT_RULES), None)
    assert prevalence_values(rules) == [[2]]


def test_resolve_without_design_rejects_group_thresholds(filter_caplog):
    with pytest.raises(SystemExit):
        resolve_prevalence_rules(rules_from(GROUP_RULES), None)
    assert "requires a counts design file" in filter_caplog.text


def test_resolve_integer_applies_to_every_group():
    rules = resolve_prevalence_rules(rules_from(INT_RULES), DESIGN)
    assert prevalence_values(rules) == [[{"prevalence.K": 2, "prevalence.B": 2}]]


def test_resolve_group_thresholds():
    rules = resolve_prevalence_rules(rules_from(GROUP_RULES), DESIGN)
    assert prevalence_values(rules) == [[{"prevalence.K": 3, "prevalence.B": 1}]]


def test_resolve_keeps_other_rules():
    rules = resolve_prevalence_rules(rules_from(INT_RULES), DESIGN)
    other = rules["genic"][0][rules["genic"][0]["type"] != "min_prevalence"]
    assert other.values.tolist() == [["min_cov", "Min_Threshold", 3]]


def test_resolve_branch_with_only_integer_rules():
    """A branch whose rules are all numeric has an int64 rule column that cannot hold dicts."""
    raw = rules_from({"genic": [{"min_prevalence": 2}]})
    assert pd.api.types.is_integer_dtype(raw["genic"][0]["rule"])
    rules = resolve_prevalence_rules(raw, DESIGN)
    assert prevalence_values(rules) == [[{"prevalence.K": 2, "prevalence.B": 2}]]


def test_resolve_across_categories_and_branches():
    rules = resolve_prevalence_rules(rules_from({
        "genic": [{"min_prevalence": 1}, {"min_prevalence": {"K": 2, "B": 2}, "min_cov": 3}],
        "full-splice_match": [{"min_cov": 3}],
    }), DESIGN)
    assert prevalence_values(rules) == [[{"prevalence.K": 1, "prevalence.B": 1}],
                                        [{"prevalence.K": 2, "prevalence.B": 2}]]
    assert prevalence_values(rules, "full-splice_match") == [[]]


def test_resolve_does_not_modify_input():
    raw = rules_from(INT_RULES)
    resolve_prevalence_rules(raw, DESIGN)
    assert prevalence_values(raw) == [[2]]


@pytest.mark.parametrize("value, message", [
    ({"K": 2}, "no threshold for groups: ['B']"),
    ({"K": 2, "B": 1, "X": 1}, "groups not in the counts design file: ['X']"),
])
def test_resolve_groups_must_match_design(value, message, filter_caplog):
    with pytest.raises(SystemExit):
        resolve_prevalence_rules(rules_from({"genic": [{"min_prevalence": value}]}), DESIGN)
    assert message in filter_caplog.text
    assert "Groups in the counts design file: ['K', 'B']" in filter_caplog.text


@pytest.mark.parametrize("value, group, size", [
    (3, "B", 2),                    # integer larger than the smallest group
    ({"K": 4, "B": 1}, "K", 3),
    ({"K": 1, "B": 3}, "B", 2),
])
def test_resolve_threshold_above_group_size(value, group, size, filter_caplog):
    with pytest.raises(SystemExit):
        resolve_prevalence_rules(rules_from({"genic": [{"min_prevalence": value}]}), DESIGN)
    assert f"in group '{group}' ({size})" in filter_caplog.text


@pytest.mark.parametrize("value", [2, {"K": 3, "B": 2}])
def test_resolve_threshold_equal_to_group_size(value):
    resolve_prevalence_rules(rules_from({"genic": [{"min_prevalence": value}]}), DESIGN)  # must not exit


# ---------------------------------------------------------------- add_group_prevalence

def test_add_group_prevalence():
    classif = add_group_prevalence(make_classif(), DESIGN)
    assert classif["prevalence.K"].tolist() == [2, 0, 1, 3, 1]
    assert classif["prevalence.B"].tolist() == [0, 2, 0, 0, 1]


def test_add_group_prevalence_ignores_samples_outside_design():
    classif = add_group_prevalence(make_classif(), DESIGN)
    assert "prevalence.MIX" not in classif.columns
    assert classif["prevalence"].tolist() == [2, 2, 2, 4, 2]  # QC column untouched


def test_add_group_prevalence_missing_counts_are_not_detections():
    classif = make_classif()
    classif.loc[0, "FL.K1"] = float("nan")
    assert add_group_prevalence(classif, DESIGN)["prevalence.K"].iloc[0] == 1


def test_add_group_prevalence_uses_shared_detection_threshold(monkeypatch):
    """The filter must count detections with the same threshold as QC."""
    monkeypatch.setattr(sqanti3_rules_filter, "MIN_DETECTION_COUNT", 0.5)
    classif = add_group_prevalence(make_classif(), DESIGN)
    assert classif["prevalence.K"].iloc[2] == 2  # 0.9 and 1.0 now count, 0.4 does not


# ---------------------------------------------------------------- passes_prevalence / has_prevalence_rules

def row(**values):
    return pd.Series(values)


@pytest.mark.parametrize("prevalence, expected", [(3, True), (2, True), (1, False)])
def test_passes_prevalence_integer_uses_qc_column(prevalence, expected):
    assert passes_prevalence(row(prevalence=prevalence), 2) is expected


@pytest.mark.parametrize("k, b, expected", [
    (2, 0, True),    # only K
    (0, 3, True),    # only B
    (2, 3, True),    # both
    (1, 2, False),   # detected in both, enough in neither: OR is not a sum
    (0, 0, False),
])
def test_passes_prevalence_or_between_groups(k, b, expected):
    thresholds = {"prevalence.K": 2, "prevalence.B": 3}
    assert passes_prevalence(row(**{"prevalence.K": k, "prevalence.B": b}), thresholds) is expected


def test_passes_prevalence_all_na_fails():
    """NA prevalence means detected in no sample (isoform missing from --fl_count)."""
    assert passes_prevalence(row(prevalence=pd.NA), 2) is False
    assert passes_prevalence(row(**{"prevalence.K": pd.NA, "prevalence.B": pd.NA}),
                             {"prevalence.K": 2, "prevalence.B": 2}) is False


@pytest.mark.parametrize("b, expected", [(2, True), (1, False)])
def test_passes_prevalence_partial_na_uses_remaining_groups(b, expected):
    r = row(**{"prevalence.K": pd.NA, "prevalence.B": b})
    assert passes_prevalence(r, {"prevalence.K": 2, "prevalence.B": 2}) is expected


@pytest.mark.parametrize("rules, expected", [
    ({"genic": [{"min_cov": 3}]}, False),
    ({"genic": [{"min_cov": 3}], "rest": []}, False),
    ({"genic": [{"min_cov": 3}, {"min_prevalence": 2}]}, True),
    ({"genic": [{"min_cov": 3}], "antisense": [{"min_prevalence": {"K": 2}}]}, True),
])
def test_has_prevalence_rules(rules, expected):
    assert has_prevalence_rules(rules_from(rules)) is expected


# ---------------------------------------------------------------- apply_rules / get_reasons with groups

def resolved(rules):
    return resolve_prevalence_rules(rules_from(rules), DESIGN)


@pytest.mark.parametrize("idx, expected", [
    (0, "Isoform"),   # t1: K=2
    (1, "Isoform"),   # t2: B=2
    (2, "Artifact"),  # t3: K=1, MIX does not count
    (3, "Artifact"),  # t4: K=3 but fails min_cov
    (4, "Artifact"),  # t5: K=1, B=1
])
def test_apply_rules_group_prevalence(idx, expected):
    classif = add_group_prevalence(make_classif(), DESIGN)
    assert apply_rules(classif.iloc[idx], False, resolved(INT_RULES)) == expected


def test_get_reasons_group_prevalence():
    classif = add_group_prevalence(make_classif(), DESIGN)
    rules = resolved(INT_RULES)
    assert get_reasons(classif.iloc[4], False, rules)["filter_reason"] == \
        "prevalence.K: 1 < 2, prevalence.B: 1 < 2 Multisample-artifact"
    assert get_reasons(classif.iloc[3], False, rules)["filter_reason"] == "min_cov: 1 < 3"


def test_get_reasons_group_thresholds():
    classif = add_group_prevalence(make_classif(), DESIGN)
    assert get_reasons(classif.iloc[0], False, resolved(GROUP_RULES))["filter_reason"] == \
        "prevalence.K: 2 < 3, prevalence.B: 0 < 1 Multisample-artifact"


# ---------------------------------------------------------------- rules_filter scenarios

@pytest.mark.parametrize("rules, design, expected", [
    # Integer threshold, no design: as in the previous PR, over all samples (MIX counts)
    (INT_RULES,   None,   ["t1", "t2", "t3", "t5"]),
    # Same threshold, with design: per group, MIX left out, OR between groups
    (INT_RULES,   DESIGN, ["t1", "t2"]),
    # One threshold per group
    (GROUP_RULES, DESIGN, ["t2", "t5"]),
    # Per-group thresholds without design
    (GROUP_RULES, None,   SystemExit),
    # Design groups do not match the rules
    ({"genic": [{"min_prevalence": {"K": 2}}], "rest": []}, DESIGN, SystemExit),
    # Threshold valid for all samples (6) but not for group K (3)
    ({"genic": [{"min_prevalence": 4}], "rest": []}, None,   ["t4"]),
    ({"genic": [{"min_prevalence": 4}], "rest": []}, DESIGN, SystemExit),
    # Design sample not in the classification
    (INT_RULES, {"K": ["K1", "K9"]}, SystemExit),
])
def test_rules_filter_scenarios(tmp_path, rules, design, expected):
    class_file, json_file, design_file = write_inputs(tmp_path, make_classif(), rules, design)
    prefix = str(tmp_path / "out")
    run = lambda: rules_filter(class_file, json_file, False, prefix, filter_logger,
                               counts_design=design_file)
    if expected is SystemExit:
        with pytest.raises(SystemExit):
            run()
        assert not os.path.exists(f"{prefix}_pass_isoforms.txt")
    else:
        run()
        assert read_pass(prefix) == expected


def test_rules_filter_writes_group_prevalence_columns(tmp_path):
    class_file, json_file, design_file = write_inputs(tmp_path, make_classif(), INT_RULES, DESIGN)
    prefix = str(tmp_path / "out")
    rules_filter(class_file, json_file, False, prefix, filter_logger, counts_design=design_file)
    out = pd.read_csv(f"{prefix}_RulesFilter_classification.txt", sep="\t")
    assert out["prevalence.K"].tolist() == [2, 0, 1, 3, 1]
    assert out["prevalence.B"].tolist() == [0, 2, 0, 0, 1]
    assert out["prevalence"].tolist() == [2, 2, 2, 4, 2]


def test_rules_filter_without_design_adds_no_columns(tmp_path):
    class_file, json_file, _ = write_inputs(tmp_path, make_classif(), INT_RULES)
    prefix = str(tmp_path / "out")
    rules_filter(class_file, json_file, False, prefix, filter_logger)
    out = pd.read_csv(f"{prefix}_RulesFilter_classification.txt", sep="\t")
    assert not [c for c in out.columns if c.startswith("prevalence.")]


def test_rules_filter_reasons_file_with_design(tmp_path):
    class_file, json_file, design_file = write_inputs(tmp_path, make_classif(), INT_RULES, DESIGN)
    prefix = str(tmp_path / "out")
    rules_filter(class_file, json_file, False, prefix, filter_logger, counts_design=design_file)
    reasons = pd.read_csv(f"{prefix}_filtering_reasons.txt", sep="\t").set_index("isoform")
    assert reasons.loc["t3", "filter_reason"] == \
        "prevalence.K: 1 < 2, prevalence.B: 0 < 2 Multisample-artifact"


def test_rules_filter_logs_design(tmp_path, filter_caplog):
    class_file, json_file, design_file = write_inputs(tmp_path, make_classif(), INT_RULES, DESIGN)
    rules_filter(class_file, json_file, False, str(tmp_path / "out"), filter_logger,
                 counts_design=design_file)
    assert "Counts design: K (3 samples), B (2 samples)" in filter_caplog.text


def test_rules_filter_design_without_prevalence_rules(tmp_path, filter_caplog):
    """A design with no min_prevalence rules warns and changes nothing."""
    rules = {"genic": [{"min_cov": 3}], "rest": []}
    outputs = {}
    for with_design in (True, False):
        sub = tmp_path / str(with_design)
        sub.mkdir()
        class_file, json_file, design_file = write_inputs(
            sub, make_classif(), rules, DESIGN if with_design else None)
        prefix = str(sub / "out")
        rules_filter(class_file, json_file, False, prefix, filter_logger, counts_design=design_file)
        outputs[with_design] = pd.read_csv(f"{prefix}_RulesFilter_classification.txt", sep="\t")
    pd.testing.assert_frame_equal(outputs[True], outputs[False])
    assert "--counts_design is ignored" in filter_caplog.text


def test_rules_filter_ignore_prevalence_never_reads_design(tmp_path):
    """Rescue path: with --ignore_prevalence the design is not even opened."""
    class_file, json_file, _ = write_inputs(tmp_path, make_classif(fl=False), GROUP_RULES)
    prefix = str(tmp_path / "out")
    rules_filter(class_file, json_file, False, prefix, filter_logger, ignore_prevalence=True,
                 counts_design=str(tmp_path / "does_not_exist.json"))
    assert read_pass(prefix) == ["t1", "t2", "t3", "t5"]


# ---------------------------------------------------------------- argparse and callers

def test_argparse_counts_design_default():
    args = filter_argparse().parse_args(["rules", "--sqanti_class", "x"])
    assert args.counts_design is None


def test_argparse_counts_design_value():
    args = filter_argparse().parse_args(["rules", "--sqanti_class", "x", "--counts_design", "d.json"])
    assert args.counts_design == "d.json"


@pytest.mark.parametrize("design", [None, "design.json"])
def test_run_rules_passes_counts_design(tmp_path, design):
    from src import filter_steps
    args = filter_argparse().parse_args(
        ["rules", "--sqanti_class", "x", "-d", str(tmp_path), "-o", "out", "--skip_report"]
        + (["--counts_design", design] if design else []))
    (tmp_path / "out_pass_isoforms.txt").write_text("")
    with patch.object(filter_steps, "rules_filter") as mock_rf:
        filter_steps.run_rules(args)
    assert mock_rf.call_args.kwargs["counts_design"] == design


def validation_args(tmp_path, *extra):
    class_file = tmp_path / "class.txt"
    class_file.write_text("")
    return filter_argparse().parse_args(
        ["rules", "--sqanti_class", str(class_file), "-d", str(tmp_path / "out")] + list(extra))


def test_filter_args_validation_missing_design(tmp_path):
    from src.argparse_utils import filter_args_validation
    args = validation_args(tmp_path, "--counts_design", str(tmp_path / "missing.json"))
    with pytest.raises(SystemExit):
        filter_args_validation(args)


def test_filter_args_validation_existing_design(tmp_path):
    from src.argparse_utils import filter_args_validation
    design = tmp_path / "design.json"
    design.write_text(json.dumps(DESIGN))
    filter_args_validation(validation_args(tmp_path, "--counts_design", str(design)))  # must not exit


@pytest.mark.parametrize("design", [None, "design.json"])
def test_write_filter_parameters_counts_design(tmp_path, design):
    from src.write_parameters import write_filter_parameters
    args = filter_argparse().parse_args(
        ["rules", "--sqanti_class", "x", "-d", str(tmp_path), "-o", "out"]
        + (["--counts_design", design] if design else []))
    write_filter_parameters(args)
    params = dict(line.rstrip("\n").split("\t", 1)
                  for line in open(tmp_path / "out_params.txt"))
    assert params["CountsDesign"] == (os.path.abspath(design) if design else "NA")


@pytest.mark.parametrize("value, expected", [("", ""), ("d.json", "--counts_design d.json")])
def test_wrapper_serializes_counts_design(value, expected):
    """The empty string in the YAML must be omitted, a path must be passed through."""
    from src.wrapper_utils import format_options
    assert format_options({"counts_design": value}) == expected