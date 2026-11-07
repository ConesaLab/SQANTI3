"""
test_auto_filter.py

Unit tests of the automatic filter of sqanti3_filter.py: default
subcommand, graded pubLR, calibrated random forest cut and the
fallbacks to graded pubLR and to the rules.

Version: 1.0
Author: Adrian Berenguer-Agustí (I2SysBio-CSIC)
Date: 2026-11-07
"""

import json
import logging
import os
import sys

import numpy as np
import pandas as pd

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.filter_argparse import default_subcommand  # noqa: E402
from src.utilities.filter.sqanti3_auto_filter import (  # noqa: E402
    ROUTE_ML, ROUTE_PUBLR, ROUTE_RULES, auto_decisions, auto_filter,
    calibrated_cuts
)

NIC = "novel_in_catalog"
NNC = "novel_not_in_catalog"
LOGGER = logging.getLogger("test_auto_filter")


def classification(rows: list[tuple]) -> pd.DataFrame:
    """Classification of (category, exons, perc_A, pubLR_prob, label)."""
    table = pd.DataFrame(rows, columns=["structural_category", "exons",
                                        "perc_A_downstream_TTS",
                                        "pubLR_prob", "pantranscriptome"])
    table.insert(0, "isoform", [f"t{i}" for i in range(len(table))])
    result = table

    return result


def test_default_subcommand() -> None:
    assert default_subcommand(["--sqanti_class", "c.txt"])[0] == "auto"
    assert default_subcommand([]) == ["auto"]
    for argv in (["rules", "--sqanti_class", "c.txt"], ["ml"], ["-h"]):
        assert default_subcommand(argv) == argv


def test_graded_publr() -> None:
    classif = classification([
        (NIC, 3, 10, 0.80, "seen"),        # sure: kept without rules
        (NIC, 3, 10, 0.60, "seen"),        # unsure, passes the rules
        (NNC, 3, 10, 0.60, "unseen"),      # unsure, fails the rules
        (NNC, 3, 10, 0.40, "recurrent"),   # below the threshold
        (NIC, 3, 80, 0.90, "recurrent"),   # intra-primed
        (NIC, 1, 10, np.nan, None),        # mono-exon: rules
        ("full-splice_match", 3, 10, np.nan, None),
    ])
    keep_rules = pd.Series([False, True, False, True, True, True, False])
    decisions = auto_decisions(classif, keep_rules, None, 60, LOGGER)
    assert decisions["keep"].tolist() == [True, True, False, False, False,
                                          True, False]
    assert decisions["filter_route"].tolist() == [ROUTE_PUBLR] * 5 + [
        ROUTE_RULES] * 2
    assert decisions.loc[3, "reason"].startswith("pubLR_prob: 0.4 <")
    assert decisions.loc[4, "reason"] == "Intra-priming"
    assert decisions.loc[2, "reason"] == ""


def test_calibrated_ml() -> None:
    # 20 clean recurrent NIC with POS_MLprob 0.05-1.00: the cut (0.10)
    # keeps 19 of them, at least 90 %. NNC has too few recurrent isoforms and stays graded.
    rows = [(NIC, 3, 10, 0.9, "recurrent")] * 20
    rows += [(NIC, 3, 10, 0.9, "unseen"), (NNC, 3, 10, 0.8, "recurrent")]
    classif = classification(rows)
    ml_prob = pd.Series(list(np.arange(1, 21) / 20) + [0.05, 0.01])
    clean = pd.Series(True, index=classif.index)
    cuts = calibrated_cuts(classif, ml_prob, clean)
    assert set(cuts) == {NIC} and cuts[NIC] == 0.10
    keep_rules = pd.Series(False, index=classif.index)
    decisions = auto_decisions(classif, keep_rules, ml_prob, 60, LOGGER)
    assert (decisions["filter_route"][:21] == ROUTE_ML).all()
    assert decisions["keep"][:21].sum() == 19
    assert decisions.loc[20, "reason"] == "POS_MLprob: 0.05 < 0.1"
    assert decisions.loc[21, "filter_route"] == ROUTE_PUBLR
    assert decisions.loc[21, "keep"]


def test_untrained_forest_is_graded() -> None:
    classif = classification([(NIC, 3, 10, 0.3, "recurrent")] * 25)
    ml_prob = pd.Series(np.nan, index=classif.index)
    keep_rules = pd.Series(True, index=classif.index)
    decisions = auto_decisions(classif, keep_rules, ml_prob, 60, LOGGER)
    assert (decisions["filter_route"] == ROUTE_PUBLR).all()
    assert not decisions["keep"].any()


def test_auto_filter_files(tmp_path) -> None:
    classif = classification([(NIC, 3, 10, 0.3, "seen"),
                              (NIC, 3, 10, 0.9, "seen"),
                              ("full-splice_match", 3, 80, np.nan, None)])
    classif["all_canonical"] = "canonical"
    class_path = tmp_path / "x_classification.txt"
    classif.to_csv(class_path, sep="\t", index=False)
    json_path = tmp_path / "rules.json"
    json_path.write_text(json.dumps({
        "full-splice_match": [{"perc_A_downstream_TTS": [0, 59]}],
        "rest": [{"all_canonical": "canonical"}]}))
    prefix = str(tmp_path / "x")
    auto_filter(str(class_path), str(json_path), False, 60, None, prefix,
                LOGGER)
    out = pd.read_csv(f"{prefix}_AutoFilter_classification.txt", sep="\t")
    assert out["filter_result"].tolist() == ["Artifact", "Isoform",
                                             "Artifact"]
    assert out["filter_route"].tolist() == [ROUTE_PUBLR, ROUTE_PUBLR,
                                            ROUTE_RULES]
    kept = open(f"{prefix}_pass_isoforms.txt").read().split()
    assert kept == ["t1"]
    reasons = pd.read_csv(f"{prefix}_filtering_reasons.txt", sep="\t")
    assert reasons["filter_reason"].iloc[0].startswith("pubLR_prob")
    assert reasons["filter_reason"].iloc[1].startswith(
        "perc_A_downstream_TTS")


def test_without_publr_uses_rules(tmp_path, caplog) -> None:
    classif = classification([(NIC, 3, 10, np.nan, None)])
    classif["all_canonical"] = "non_canonical"
    class_path = tmp_path / "x_classification.txt"
    classif.to_csv(class_path, sep="\t", index=False)
    json_path = tmp_path / "rules.json"
    json_path.write_text(json.dumps(
        {"rest": [{"all_canonical": "canonical"}]}))
    prefix = str(tmp_path / "x")
    with caplog.at_level(logging.WARNING):
        auto_filter(str(class_path), str(json_path), False, 60, None,
                    prefix, LOGGER)
    assert "no pubLR probabilities" in caplog.text
    out = pd.read_csv(f"{prefix}_AutoFilter_classification.txt", sep="\t")
    assert out["filter_route"].tolist() == [ROUTE_RULES]
    assert out["filter_result"].tolist() == ["Artifact"]
