"""
Tests for the rescue argument parser and its validation, focused on the
requantification flag defaults.
"""
import argparse
import json
import logging

import pytest

from src.rescue_argparse import rescue_argparse
import src.argparse_utils as au

MIN_ARGS = ["--filter_class", "class.txt", "-rg", "ref.gtf", "-rf", "ref.fasta"]


@pytest.fixture
def parser():
    return rescue_argparse()


class TestRequantFlagParsing:
    """Requantification is on by default and --skip_requant disables it."""

    def test_enabled_by_default(self, parser):
        args = parser.parse_args(MIN_ARGS)
        assert args.skip_requant is False

    def test_skip_requant_disables_it(self, parser):
        args = parser.parse_args(MIN_ARGS + ["--skip_requant"])
        assert args.skip_requant is True

    def test_legacy_requant_flag_is_rejected(self, parser):
        with pytest.raises(SystemExit):
            parser.parse_args(MIN_ARGS + ["--requant"])


def _rescue_namespace(**overrides):
    """Minimal namespace covering every attribute the rescue validation reads."""
    ns = argparse.Namespace(
        filter_class="class.txt",
        dir="rescue_out",
        output="rescue",
        refGTF="ref.gtf",
        refFasta="ref.fasta",
        corrected_isoforms_fasta=None,
        filtered_isoforms_gtf=None,
        refClassif=None,
        counts=None,
        mode="automatic",
        strategy="rules",
        json_filter="filter.json",
        counts_design=None,
        min_expression=0.0,
        map_automatic_fsm=False,
        random_forest="rf.RData",
        threshold=0.7,
        skip_requant=False,
        skip_evidence_check=True,
    )
    for key, value in overrides.items():
        setattr(ns, key, value)
    return ns


class TestRequantValidation:
    """Missing --counts stops the run at validation, before any step, unless --skip_requant."""

    @pytest.fixture(autouse=True)
    def _no_disk_access(self, monkeypatch, caplog):
        monkeypatch.setattr(au, "valid_file", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_gtf", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_fasta", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_dir", lambda *a, **k: True)
        # The rescue logger does not propagate to root, so caplog cannot see
        # its records unless propagation is temporarily restored.
        monkeypatch.setattr(au.rescue_logger, "propagate", True)
        caplog.set_level(logging.WARNING, logger=au.rescue_logger.name)

    def test_missing_counts_exits(self, caplog):
        args = _rescue_namespace(skip_requant=False, counts=None)
        with pytest.raises(SystemExit):
            au.rescue_args_validation(args)
        assert "--skip_requant" in caplog.text

    def test_counts_provided_keeps_requant_on(self, monkeypatch):
        checked = []
        monkeypatch.setattr(au, "valid_file", lambda filename, logger: checked.append(filename))
        args = _rescue_namespace(skip_requant=False, counts="counts.tsv")
        au.rescue_args_validation(args)
        assert args.skip_requant is False
        assert "counts.tsv" in checked

    def test_skip_requant_without_counts_passes(self, monkeypatch):
        checked = []
        monkeypatch.setattr(au, "valid_file", lambda filename, logger: checked.append(filename))
        args = _rescue_namespace(skip_requant=True, counts=None)
        au.rescue_args_validation(args)
        assert None not in checked

    def test_skip_requant_still_checks_given_counts(self, monkeypatch):
        """The evidence check also reads --counts, so a given file is validated even without requantification."""
        checked = []
        monkeypatch.setattr(au, "valid_file", lambda filename, logger: checked.append(filename))
        au.rescue_args_validation(_rescue_namespace(skip_requant=True, counts="counts.tsv"))
        assert "counts.tsv" in checked


class TestMapAutomaticFsm:
    """--map_automatic_fsm is off by default and needs the corrected FASTA."""

    @pytest.fixture(autouse=True)
    def _no_disk_access(self, monkeypatch):
        monkeypatch.setattr(au, "valid_file", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_gtf", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_fasta", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_dir", lambda *a, **k: True)

    def test_disabled_by_default(self, parser):
        args = parser.parse_args(MIN_ARGS)
        assert args.map_automatic_fsm is False

    def test_hidden_from_help_and_config(self, parser):
        """Beta feature: it works on the command line but is not shown to users."""
        from src.wrapper_utils import get_parser_specific_args_simple
        assert "--map_automatic_fsm" not in parser.format_help()
        assert "map_automatic_fsm" not in get_parser_specific_args_simple(parser, {})["options"]

    def test_flag_enables_it(self, parser):
        args = parser.parse_args(MIN_ARGS + ["--map_automatic_fsm"])
        assert args.map_automatic_fsm is True

    def test_without_corrected_fasta_exits(self):
        args = _rescue_namespace(map_automatic_fsm=True, counts="counts.tsv")
        with pytest.raises(SystemExit):
            au.rescue_args_validation(args)

    def test_with_corrected_fasta_passes(self):
        args = _rescue_namespace(map_automatic_fsm=True, counts="counts.tsv",
                                 corrected_isoforms_fasta="corrected.fasta")
        au.rescue_args_validation(args)
        assert args.map_automatic_fsm is True


class TestEvidenceCheckValidation:
    """The evidence check is decided at validation; it reads --counts even with --skip_requant."""

    @pytest.fixture(autouse=True)
    def _no_disk_access(self, monkeypatch):
        monkeypatch.setattr(au, "valid_gtf", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_fasta", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_dir", lambda *a, **k: True)
        monkeypatch.setattr(au, "valid_file", lambda *a, **k: True)

    @pytest.fixture
    def json_file(self, tmp_path):
        p = tmp_path / "rules.json"
        p.write_text(json.dumps({"full-splice_match": [{"prevalence_K": 2}, {"prevalence_B": 2}]}))
        return str(p)

    def args(self, json_file, **kw):
        base = dict(skip_evidence_check=False, json_filter=json_file, counts="counts.tsv")
        base.update(kw)
        return _rescue_namespace(**base)

    def test_enabled(self, json_file):
        assert au.rescue_args_validation(self.args(json_file)) is True

    def test_disabled_by_flag(self, json_file):
        assert au.rescue_args_validation(self.args(json_file, skip_evidence_check=True)) is False

    def test_disabled_for_ml(self, json_file):
        assert au.rescue_args_validation(self.args(json_file, strategy="ml")) is False

    def test_disabled_without_count_rules(self, tmp_path):
        p = tmp_path / "r.json"
        p.write_text(json.dumps({"rest": [{"min_cov": 3}]}))
        assert au.rescue_args_validation(self.args(str(p))) is False

    def test_count_rules_without_counts_exit(self, json_file):
        """With --skip_requant, --counts is still needed by the evidence check."""
        with pytest.raises(SystemExit):
            au.rescue_args_validation(self.args(json_file, counts=None, skip_requant=True))

    def test_skip_both_without_counts_passes(self, json_file):
        args = self.args(json_file, counts=None, skip_requant=True, skip_evidence_check=True)
        assert au.rescue_args_validation(args) is False
