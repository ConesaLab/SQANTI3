"""
Unit tests for evidence_check.py: the requisites on long-read counts of the
rules applied to the counts that rescued reference transcripts receive from
their artifacts, distributed as requantification does.
"""
import json
import os
import sys
import tempfile
from types import SimpleNamespace

import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.evidence_check import (
    is_count_rule, target_count_rules, count_samples, target_counts,
    add_count_columns, failed_targets, best_lr_fallback, write_evidence_rules,
    check_design_for_rules, apply_evidence_check, rows_for_requant
)
from src.rescue_steps import evidence_check_enabled, run_evidence_check, run_fallback_mapping

K = ["K1", "K2", "K3"]
B = ["B1", "B2", "B3"]
MIX = ["M1", "M2"]
SAMPLES = K + B + MIX
FL = [f"FL.{s}" for s in SAMPLES]
DESIGN = {"K": K, "B": B}

# Detected in at least 2 samples of any group (OR between rules)
ANY_GROUP = [{"prevalence_K": 2}, {"prevalence_B": 2}]
# Detected in at least 2 samples of every group (AND within a rule)
EVERY_GROUP = [{"prevalence_K": 2, "prevalence_B": 2}]
GLOBAL = [{"prevalence": 2}]


def parsed(rule_sets):
    """Count requisites parsed by the rules filter, as run_evidence_check() does."""
    with tempfile.TemporaryDirectory() as d:
        return write_evidence_rules(rule_sets, os.path.join(d, "rules.json"))


def classif(rows):
    """Filter classification: {isoform: (filter_result, {sample: count})}."""
    data = []
    for iso, (result, counts) in rows.items():
        rec = {"isoform": iso, "filter_result": result, "associated_gene": "G1"}
        rec.update({f"FL.{s}": counts.get(s, 0) for s in SAMPLES})
        data.append(rec)
    return pd.DataFrame(data, columns=["isoform", "filter_result", "associated_gene"] + FL)


def counts_of(c):
    """Count matrix of requantification (--counts) matching a classification."""
    return c[["isoform"] + FL].rename(columns=lambda col: col[len("FL."):] if col.startswith("FL.") else col)


def counts_file(tmp_path, c):
    path = tmp_path / "counts.tsv"
    counts_of(c).to_csv(path, sep="\t", index=False)
    return str(path)


def rescue_table(rows):
    """rows: (artifact, assigned_transcript, rescue_mode, origin, reintroduced)."""
    return pd.DataFrame(rows, columns=["artifact", "assigned_transcript", "rescue_mode",
                                       "origin", "reintroduced"])


def hits(rows):
    return pd.DataFrame(rows, columns=["rescue_candidate", "mapping_hit", "alignment_type",
                                       "alignment_score"])


class TestTargetCountRules:
    @pytest.mark.parametrize("column", ["FL", "FL.K1", "prevalence", "prevalence_K"])
    def test_count_columns(self, column):
        assert is_count_rule(column)

    @pytest.mark.parametrize("column", ["FLAIR", "min_cov", "subcategory", "prevalent"])
    def test_other_columns(self, column):
        assert not is_count_rule(column)

    def test_keeps_only_count_requisites(self):
        rules = {"full-splice_match": [{"perc_A_downstream_TTS": [0, 59], "prevalence_K": 2},
                                       {"prevalence_B": 2, "FL": 3}]}
        assert target_count_rules(rules) == [{"prevalence_K": 2}, {"prevalence_B": 2, "FL": 3}]

    def test_falls_back_to_rest(self):
        assert target_count_rules({"rest": [{"prevalence": 3}]}) == [{"prevalence": 3}]

    def test_rule_without_count_requisites_makes_them_unnecessary(self):
        # OR between rules: a rule without counts can be passed without them
        rules = {"full-splice_match": [{"prevalence": 2}, {"perc_A_downstream_TTS": [0, 59]}]}
        assert target_count_rules(rules) == []

    def test_no_rules(self):
        assert target_count_rules({"rest": []}) == []


class TestTargetCounts:
    def test_sums_artifacts_of_the_same_target(self):
        c = classif({"a1": ("Artifact", {"K1": 1}), "a2": ("Artifact", {"K2": 2})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a2", "REF1", "rules_mapping", "reference", "no")])
        agg = target_counts(r, c, counts_of(c))
        assert agg.loc["REF1", "FL.K1"] == 1
        assert agg.loc["REF1", "FL.K2"] == 2
        assert agg.loc["REF1"].sum() == 3

    def test_integer_counts_stay_integer(self):
        # As in requantification: floor of each share, remainder to the first target
        c = classif({"a1": ("Artifact", {"K1": 1, "K2": 3})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a1", "REF2", "rules_mapping", "reference", "yes")])
        agg = target_counts(r, c, counts_of(c))
        assert agg.loc["REF1", ["FL.K1", "FL.K2"]].tolist() == [1, 2]
        assert agg.loc["REF2", ["FL.K1", "FL.K2"]].tolist() == [0, 1]

    def test_fractional_counts_split_proportionally(self):
        c = classif({"a1": ("Artifact", {"K1": 0.6})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a1", "REF2", "rules_mapping", "reference", "yes")])
        agg = target_counts(r, c, counts_of(c))
        assert agg.loc["REF1", "FL.K1"] == pytest.approx(0.3)
        assert agg.loc["REF2", "FL.K1"] == pytest.approx(0.3)

    def test_tie_with_a_long_read_isoform_with_counts(self):
        # Shares follow the counts of the targets: the reference has none
        c = classif({"a1": ("Artifact", {"K1": 4, "K2": 4}), "LR1": ("Isoform", {"K1": 10})})
        r = rescue_table([("a1", "LR1", "rules_mapping", "lr_defined", "no"),
                          ("a1", "REF1", "rules_mapping", "reference", "yes")])
        agg = target_counts(r, c, counts_of(c))
        assert agg.loc["REF1", "FL.K1"] == 0      # LR1 has counts in K1: all to LR1
        assert agg.loc["REF1", "FL.K2"] == 2      # no target has counts in K2: even split

    def test_only_reference_targets_are_returned(self):
        c = classif({"a1": ("Artifact", {"K1": 5}), "LR1": ("Isoform", {"K1": 1})})
        r = rescue_table([("a1", "LR1", "rules_mapping", "lr_defined", "no")])
        assert target_counts(r, c, counts_of(c)).empty

    def test_artifact_without_counts_is_zero(self):
        c = classif({"a1": ("Artifact", {})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        assert target_counts(r, c, counts_of(c)).loc["REF1"].sum() == 0

    def test_single_sample_column_is_fl(self):
        c = pd.DataFrame({"isoform": ["a1"], "filter_result": ["Artifact"], "associated_gene": ["G1"]})
        counts = pd.DataFrame({"isoform": ["a1"], "count": [4]})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        agg = target_counts(r, c, counts)
        assert list(agg.columns) == ["FL"] and agg.loc["REF1", "FL"] == 4
        assert count_samples(counts) == ["count"]


class TestCountColumns:
    def test_columns_as_qc_computes_them(self):
        agg = pd.DataFrame({"FL.K1": [1.0], "FL.K2": [0.5], "FL.B1": [2.0], "FL.M1": [3.0]},
                           index=["REF1"])
        out = add_count_columns(agg, {"K": ["K1", "K2"], "B": ["B1"]}).loc["REF1"]
        assert out["FL"] == 6.5
        assert out["prevalence"] == 4          # any count above 0 is expression
        assert out["prevalence_K"] == 2        # half a read included
        assert out["prevalence_B"] == 1

    def test_min_expression_as_in_qc(self):
        agg = pd.DataFrame({"FL.K1": [1.0], "FL.K2": [0.5], "FL.B1": [2.0], "FL.M1": [0.0]},
                           index=["REF1"])
        out = add_count_columns(agg, {"K": ["K1", "K2"], "B": ["B1"]}, min_expression=1).loc["REF1"]
        assert out["prevalence"] == 2          # K1 and B1 reach 1
        assert (out["prevalence_K"], out["prevalence_B"]) == (1, 1)

    def test_single_sample(self):
        agg = add_count_columns(pd.DataFrame({"FL": [4.0]}, index=["REF1"]))
        assert agg.loc["REF1", "FL"] == 4 and agg.loc["REF1", "prevalence"] == 1


class TestFailedTargets:
    def test_detections_in_different_samples_add_up(self):
        # Neither ISM reaches prevalence_K = 2 alone, together they do
        c = classif({"ism1": ("Artifact", {"K1": 1}), "ism2": ("Artifact", {"K2": 2})})
        r = rescue_table([("ism1", "REF1", "rules_mapping", "reference", "yes"),
                          ("ism2", "REF1", "rules_mapping", "reference", "no")])
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN) == set()

    def test_counts_in_the_same_sample_do_not_create_reproducibility(self):
        c = classif({"ism1": ("Artifact", {"K1": 1}), "ism2": ("Artifact", {"K1": 2})})
        r = rescue_table([("ism1", "REF1", "rules_mapping", "reference", "yes"),
                          ("ism2", "REF1", "rules_mapping", "reference", "no")])
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN) == {"REF1"}

    def test_mixture_only_signal_fails_per_group_rules(self):
        # Detected in the mixtures, zero in both pure groups
        c = classif({"a1": ("Artifact", {"M1": 3, "M2": 4})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN) == {"REF1"}

    def test_mixture_only_signal_passes_global_prevalence(self):
        # The global prevalence counts every sample, mixtures included
        c = classif({"a1": ("Artifact", {"M1": 3, "M2": 4})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        assert failed_targets(r, c, counts_of(c), parsed(GLOBAL), None) == set()

    def test_or_between_groups(self):
        c = classif({"a1": ("Artifact", {"B1": 1, "B2": 1})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN) == set()

    def test_and_between_groups(self):
        c = classif({"a1": ("Artifact", {"B1": 1, "B2": 1}),
                     "a2": ("Artifact", {"B1": 1, "B2": 1, "K1": 1, "K3": 2})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes"),
                          ("a2", "REF2", "automatic", "reference", "yes")])
        assert failed_targets(r, c, counts_of(c), parsed(EVERY_GROUP), DESIGN) == {"REF1"}

    def test_fl_requisite_on_the_sum(self):
        c = classif({"a1": ("Artifact", {"K1": 2}), "a2": ("Artifact", {"B1": 2})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a2", "REF1", "rules_mapping", "reference", "no"),
                          ("a1", "REF2", "rules_mapping", "reference", "yes")])
        # a1 is split 1 + 1: REF1 receives 1 + 2 = 3 reads, REF2 1
        assert failed_targets(r, c, counts_of(c), parsed([{"FL": 3}]), None) == {"REF2"}

    def test_integer_read_is_not_split_into_detections(self):
        # One read in K1 and K2 tied between two references: requantification gives
        # it to the first one, so only REF1 is detected in two K samples
        c = classif({"a1": ("Artifact", {"K1": 1, "K2": 1})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a1", "REF2", "rules_mapping", "reference", "yes")])
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN) == {"REF2"}

    def test_fractional_shares_are_expression_by_default(self):
        # Fractional counts (e.g. bambu) are split proportionally: 0.3 each
        c = classif({"a1": ("Artifact", {"K1": 0.6, "K2": 0.6})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a1", "REF2", "rules_mapping", "reference", "yes")])
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN) == set()
        assert failed_targets(r, c, counts_of(c), parsed(ANY_GROUP), DESIGN,
                              min_expression=0.5) == {"REF1", "REF2"}


class TestDesignForRules:
    def test_group_rules_need_a_design(self):
        with pytest.raises(SystemExit):
            check_design_for_rules(ANY_GROUP, None)

    def test_group_missing_from_design(self):
        with pytest.raises(SystemExit):
            check_design_for_rules([{"prevalence_X": 2}], DESIGN)

    def test_global_rules_need_no_design(self):
        check_design_for_rules(GLOBAL, None)
        check_design_for_rules(EVERY_GROUP, DESIGN)


class TestBestLrFallback:
    def test_keeps_best_supported_hit_and_ties(self):
        h = hits([("a1", "LR1", 0, 90), ("a1", "LR2", 256, 90), ("a1", "LR3", 256, 80),
                  ("a1", "REF1", 0, 100), ("a1", "LR_ART", 256, 95)])
        fb = best_lr_fallback(["a1"], h, ["LR1", "LR2", "LR3"])
        assert set(fb["assigned_transcript"]) == {"LR1", "LR2"}

    def test_no_supported_hit(self):
        h = hits([("a1", "REF1", 0, 100)])
        assert best_lr_fallback(["a1"], h, ["LR1"]).empty


class TestApplyEvidenceCheck:
    @pytest.fixture
    def setup(self):
        c = classif({
            "fsm1": ("Artifact", {"M1": 2}),            # REF1: mixture only -> fails
            "ism1": ("Artifact", {"K1": 1}),            # REF2: K1 + K2 -> passes
            "ism2": ("Artifact", {"K2": 1}),
            "nic1": ("Artifact", {"M2": 1}),            # REF3: fails, no supported hit
            "LR1": ("Isoform", {"K1": 10, "K2": 8}),
            "LR_ART": ("Artifact", {}),
            "ism9": ("Artifact", {"K3": 1}),
        })
        r = rescue_table([
            ("fsm1", "REF1", "automatic", "reference", "yes"),
            ("ism1", "REF2", "rules_mapping", "reference", "yes"),
            ("ism2", "REF2", "rules_mapping", "reference", "no"),
            ("nic1", "REF3", "rules_mapping", "reference", "yes"),
            ("ism9", "LR1", "rules_mapping", "lr_defined", "no"),
        ])
        h = hits([("fsm1", "REF1", 0, 100), ("fsm1", "LR1", 256, 85),
                  ("nic1", "REF3", 0, 100), ("nic1", "LR_ART", 256, 90)])
        inclusion = pd.Series(["REF1", "REF2", "REF3"])
        return apply_evidence_check(r, inclusion, c, counts_of(c), h, parsed(ANY_GROUP), DESIGN)

    def test_failed_targets_leave_the_inclusion_list(self, setup):
        inclusion, _ = setup
        assert list(inclusion) == ["REF2"]

    def test_failed_rows_are_kept_and_not_reintroduced(self, setup):
        _, df = setup
        failed = df[df["evidence_check"] == "failed"]
        assert set(failed["assigned_transcript"]) == {"REF1", "REF3"}
        assert (failed["reintroduced"] == "no").all()

    def test_artifact_reassigned_to_supported_lr_hit(self, setup):
        _, df = setup
        re = df[df["evidence_check"] == "reassigned"]
        assert re[["artifact", "assigned_transcript", "rescue_mode", "origin"]].values.tolist() == \
            [["fsm1", "LR1", "automatic", "lr_defined"]]

    def test_artifact_without_supported_hit_goes_to_residual(self, setup):
        _, df = setup
        req = rows_for_requant(df)
        assert "nic1" not in set(req["artifact"])

    def test_passing_and_lr_rows_untouched(self, setup):
        _, df = setup
        assert set(df.loc[df["evidence_check"] == "pass", "assigned_transcript"]) == {"REF2"}
        assert df.loc[df["origin"].eq("lr_defined") & df["evidence_check"].eq("not_required"),
                      "artifact"].tolist() == ["ism9"]

    def test_automatic_fsm_without_fallback_mapping_goes_to_residual(self):
        # Without --map_automatic_fsm the FSM artifacts of automatic rescue have
        # no hits, even if a long-read isoform of the gene passed the filter
        c = classif({"fsm1": ("Artifact", {"M1": 2}), "LR1": ("Isoform", {"K1": 10, "K2": 8})})
        r = rescue_table([("fsm1", "REF1", "automatic", "reference", "yes")])
        inc, df = apply_evidence_check(r, pd.Series(["REF1"]), c, counts_of(c), hits([]),
                                       parsed(ANY_GROUP), DESIGN)
        assert inc.empty
        assert "reassigned" not in set(df["evidence_check"])
        assert "fsm1" not in set(rows_for_requant(df)["artifact"])

    def test_artifact_with_a_passing_reference_is_not_reassigned(self):
        # a1 is tied between REF1 and REF2; a3 adds K3 to REF1 only.
        # With prevalence_K >= 3, REF1 passes (K1, K2, K3) and REF2 fails (K1, K2).
        c = classif({"a1": ("Artifact", {"K1": 2, "K2": 2}), "a2": ("Artifact", {}),
                     "a3": ("Artifact", {"K3": 1}), "LR1": ("Isoform", {"K1": 1})})
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes"),
                          ("a1", "REF2", "rules_mapping", "reference", "yes"),
                          ("a2", "REF2", "rules_mapping", "reference", "no"),
                          ("a3", "REF1", "rules_mapping", "reference", "no")])
        h = hits([("a1", "LR1", 256, 90), ("a2", "LR1", 256, 90)])
        inc, df = apply_evidence_check(r, pd.Series(["REF1", "REF2"]), c, counts_of(c), h,
                                       parsed([{"prevalence_K": 3}]), DESIGN)
        assert list(inc) == ["REF1"]
        assert set(df.loc[df["evidence_check"] == "reassigned", "artifact"]) == {"a2"}

    def test_empty_automatic_inclusion_dataframe(self):
        # run_automatic_rescue returns an empty DataFrame when nothing is lost
        r = rescue_table([])
        c = classif({})
        out, df = apply_evidence_check(r, pd.DataFrame(columns=["isoform"]), c, counts_of(c), hits([]),
                                       parsed(ANY_GROUP), DESIGN)
        assert out.empty and "evidence_check" in df.columns


class TestRowsForRequant:
    def test_without_column_passes_through(self):
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        assert rows_for_requant(r) is r


class TestOrchestration:
    @pytest.fixture
    def json_file(self, tmp_path):
        p = tmp_path / "rules.json"
        p.write_text(json.dumps({"full-splice_match": [
            {"perc_A_downstream_TTS": [0, 59], "prevalence_K": 2},
            {"perc_A_downstream_TTS": [0, 59], "prevalence_B": 2}], "rest": []}))
        return str(p)

    @pytest.fixture
    def design_file(self, tmp_path):
        p = tmp_path / "design.json"
        p.write_text(json.dumps(DESIGN))
        return str(p)

    def args(self, json_file, **kw):
        base = dict(skip_evidence_check=False, strategy="rules", json_filter=json_file, counts_design=None,
                    counts="counts.tsv")
        base.update(kw)
        return SimpleNamespace(**base)

    def test_enabled(self, json_file):
        assert evidence_check_enabled(self.args(json_file))

    def test_disabled_by_flag(self, json_file):
        assert not evidence_check_enabled(self.args(json_file, skip_evidence_check=True))

    def test_disabled_for_ml(self, json_file):
        assert not evidence_check_enabled(self.args(json_file, strategy="ml"))

    def test_disabled_without_count_rules(self, tmp_path):
        p = tmp_path / "r.json"
        p.write_text(json.dumps({"rest": [{"min_cov": 3}]}))
        assert not evidence_check_enabled(self.args(str(p)))

    def test_count_rules_without_counts_file_exit(self, json_file):
        with pytest.raises(SystemExit):
            evidence_check_enabled(self.args(json_file, counts=None))

    def test_run_evidence_check_with_design(self, json_file, design_file, tmp_path):
        c = classif({"a1": ("Artifact", {"M1": 5, "M2": 5})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        inc, df = run_evidence_check(c, r, pd.Series(["REF1"]), hits([]), json_file, design_file,
                                     counts_file(tmp_path, c), str(tmp_path))
        assert inc.empty
        assert df["evidence_check"].tolist() == ["failed"]
        # Only the count requisites are evaluated, written for traceability
        written = json.loads((tmp_path / "evidence_check_rules.json").read_text())
        assert written == {"full-splice_match": [{"prevalence_K": 2}, {"prevalence_B": 2}]}

    def test_run_evidence_check_without_design_exits(self, json_file, tmp_path):
        c = classif({"a1": ("Artifact", {"K1": 1})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        with pytest.raises(SystemExit):
            run_evidence_check(c, r, pd.Series(["REF1"]), hits([]), json_file, None,
                               counts_file(tmp_path, c), str(tmp_path))

    def test_run_evidence_check_matches_main_unpacking(self, json_file, design_file, tmp_path):
        # main() unpacks (inclusion_list, rescue_df), as for the other rescue steps,
        # and save_rescue_results() calls to_frame() on the inclusion list
        c = classif({"a1": ("Artifact", {"K1": 5, "K2": 5})})
        r = rescue_table([("a1", "REF1", "automatic", "reference", "yes")])
        inclusion_list, rescue_df = run_evidence_check(c, r, pd.Series(["REF1"]), hits([]),
                                                       json_file, design_file, counts_file(tmp_path, c),
                                                       str(tmp_path))
        assert isinstance(inclusion_list, pd.Series)
        assert isinstance(rescue_df, pd.DataFrame)
        assert inclusion_list.to_frame(name="isoform")["isoform"].tolist() == ["REF1"]

    def test_fallback_mapping_without_automatic_rows(self, tmp_path):
        r = rescue_table([("a1", "REF1", "rules_mapping", "reference", "yes")])
        out = run_fallback_mapping(classif({}), r, None, None, None, str(tmp_path), "x")
        assert out.empty

    def test_fallback_mapping_reuses_existing_hits(self, tmp_path):
        hits([("fsm1", "LR1", 0, 50)]).to_csv(tmp_path / "x_fallback_rescue_mapping_hits.tsv",
                                              sep="\t", index=False)
        r = rescue_table([("fsm1", "REF1", "automatic", "reference", "yes")])
        out = run_fallback_mapping(classif({}), r, None, None, None, str(tmp_path), "x")
        assert out["mapping_hit"].tolist() == ["LR1"]
