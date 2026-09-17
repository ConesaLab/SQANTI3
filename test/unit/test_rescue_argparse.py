"""
Tests for the rescue argument parser and its validation, focused on the
requantification flag defaults.
"""
import argparse
import logging

import pytest

from src.rescue_argparse import rescue_argparse
import src.argparse_utils as au

MIN_ARGS = ["--filter_class", "class.txt", "-rg", "ref.gtf", "-rf", "ref.fasta"]


@pytest.fixture
def parser():
    return rescue_argparse()


class TestRequantFlagParsing:
    """--requant must be on by default and still accept its legacy forms."""

    def test_enabled_by_default(self, parser):
        args = parser.parse_args(MIN_ARGS)
        assert args.requant is True

    def test_long_flag_still_accepted(self, parser):
        args = parser.parse_args(MIN_ARGS + ["--requant"])
        assert args.requant is True

    def test_short_flag_still_accepted(self, parser):
        args = parser.parse_args(MIN_ARGS + ["-q"])
        assert args.requant is True

    def test_no_requant_disables_it(self, parser):
        args = parser.parse_args(MIN_ARGS + ["--no-requant"])
        assert args.requant is False


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
        random_forest="rf.RData",
        threshold=0.7,
        requant=True,
    )
    for key, value in overrides.items():
        setattr(ns, key, value)
    return ns


class TestRequantValidation:
    """Missing --counts must warn and downgrade, never abort."""

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

    def test_missing_counts_does_not_exit(self, caplog):
        args = _rescue_namespace(requant=True, counts=None)
        with caplog.at_level(logging.WARNING):
            au.rescue_args_validation(args)
        assert args.requant is False

    def test_missing_counts_warns_explicitly(self, caplog):
        args = _rescue_namespace(requant=True, counts=None)
        with caplog.at_level(logging.WARNING):
            au.rescue_args_validation(args)
        assert any("NOT requantified" in record.message for record in caplog.records)

    def test_counts_provided_keeps_requant_on(self, caplog):
        args = _rescue_namespace(requant=True, counts="counts.tsv")
        with caplog.at_level(logging.WARNING):
            au.rescue_args_validation(args)
        assert args.requant is True
        assert not any("NOT requantified" in r.message for r in caplog.records)

    def test_no_requant_skips_counts_check(self, monkeypatch):
        checked = []
        monkeypatch.setattr(au, "valid_file", lambda filename, logger: checked.append(filename))
        args = _rescue_namespace(requant=False, counts=None)
        au.rescue_args_validation(args)
        assert args.requant is False
        assert None not in checked
