"""
sqanti3_auto_filter.py

Automatic SQANTI3 filter, the default of sqanti3_filter.py. It needs a
classification made with QC --pantranscriptome and picks the filter of
the multi-exon NIC/NNC from the data of the experiment:

  - with short-read junction coverage (min_cov), the random forest of
    the ML filter, cut within NIC and within NNC where it keeps
    ML_CALIBRATION_TARGET of the clean isoforms whose intron chain is
    recurrent in public data (calibrated ML);
  - otherwise, or when the random forest is not trained, or a category
    has fewer than ML_CALIBRATION_MIN_SUPPORTED clean recurrent
    isoforms, graded pubLR: an isoform with pubLR_prob of at least
    PUBLR_GRADED_CUT is kept, one between PUBLR_THRESHOLD and
    PUBLR_GRADED_CUT is kept only if it also passes the rules, and one
    below PUBLR_THRESHOLD is an artifact.

Intra-primed NIC/NNC are artifacts in both. Every other isoform, the
mono-exon NIC/NNC included, gets the rules filter. Without pubLR
probabilities the rules filter is applied to everything.

Version: 1.0
Author: Adrian Berenguer-Agustí (I2SysBio-CSIC)
Date: 2026-11-07
"""

import logging
import os

import numpy as np
import pandas as pd

from src.config import (
    ML_CALIBRATION_MIN_SUPPORTED, ML_CALIBRATION_TARGET,
    PANTRANSCRIPTOME_CATEGORIES, PUBLR_GRADED_CUT, PUBLR_THRESHOLD
)
from src.utilities.filter.sqanti3_rules_filter import (
    apply_rules_table, get_reasons_table, read_json_rules
)

# Routes written in the filter_route column.
ROUTE_RULES = "rules"
ROUTE_PUBLR = "pubLR_graded"
ROUTE_ML = "ML_calibrated"


def has_publr(classif: pd.DataFrame) -> bool:
    """Whether the classification carries pubLR probabilities.

    Args:
        classif: SQANTI3 classification.

    Returns:
        True if some isoform has a pubLR_prob (QC was run with
        --pantranscriptome and there are multi-exon NIC/NNC).
    """
    result = bool("pubLR_prob" in classif.columns
                  and pd.to_numeric(classif["pubLR_prob"],
                                    errors="coerce").notna().any())

    return result


def has_short_reads(classif: pd.DataFrame) -> bool:
    """Whether the junctions have short-read coverage.

    Args:
        classif: SQANTI3 classification.

    Returns:
        True if some isoform has a min_cov value (QC was run with
        --coverage or --short_reads).
    """
    result = bool("min_cov" in classif.columns
                  and pd.to_numeric(classif["min_cov"],
                                    errors="coerce").notna().any())

    return result


def calibrated_cuts(classif: pd.DataFrame, ml_prob: pd.Series,
                    clean: pd.Series) -> dict[str, float]:
    """Cut of the random forest probability within NIC and within NNC.

    Args:
        classif: SQANTI3 classification with structural_category and
            pantranscriptome.
        ml_prob: POS_MLprob of each isoform (NaN if not scored).
        clean: Whether each isoform is free of intra-priming.

    Returns:
        Category -> cut that keeps ML_CALIBRATION_TARGET of its clean
        recurrent isoforms; categories with fewer than
        ML_CALIBRATION_MIN_SUPPORTED of them are left out.
    """
    recurrent = classif["pantranscriptome"] == "recurrent"
    result = {}
    for category in PANTRANSCRIPTOME_CATEGORIES:
        labeled = ((classif["structural_category"] == category) & clean
                   & recurrent & ml_prob.notna())
        if labeled.sum() >= ML_CALIBRATION_MIN_SUPPORTED:
            result[category] = float(np.quantile(
                ml_prob[labeled], 1.0 - ML_CALIBRATION_TARGET,
                method="lower"))

    return result


def auto_decisions(classif: pd.DataFrame, keep_rules: pd.Series,
                   ml_prob: pd.Series | None, intrapriming: float,
                   logger: logging.Logger) -> pd.DataFrame:
    """Decide every isoform with the route that applies to it.

    Args:
        classif: SQANTI3 classification.
        keep_rules: Whether each isoform passes the rules filter.
        ml_prob: POS_MLprob of each isoform, or None when the random
            forest was not run.
        intrapriming: Percentage of A downstream of the TTS from which
            an isoform is intra-primed.
        logger: Logger for the routes and cuts.

    Returns:
        Columns filter_route, keep (bool) and reason (why a NIC/NNC is
        an artifact outside the rules; empty otherwise), same index.
    """
    probability = pd.Series(np.nan, index=classif.index)
    if "pubLR_prob" in classif.columns:
        probability = pd.to_numeric(classif["pubLR_prob"], errors="coerce")
    clean = pd.to_numeric(classif["perc_A_downstream_TTS"],
                          errors="coerce").fillna(0) < intrapriming
    result = pd.DataFrame({"filter_route": ROUTE_RULES, "keep": keep_rules,
                           "reason": ""}, index=classif.index)
    if probability.notna().any():
        scored = probability.notna()
        graded = (probability >= PUBLR_THRESHOLD) & clean & (
            (probability >= PUBLR_GRADED_CUT) | keep_rules)
        result.loc[scored, "filter_route"] = ROUTE_PUBLR
        result.loc[scored, "keep"] = graded[scored]
        result.loc[scored & ~graded & ~clean, "reason"] = "Intra-priming"
        low = scored & clean & (probability < PUBLR_THRESHOLD)
        result.loc[low, "reason"] = (
            "pubLR_prob: " + probability[low].round(3).astype(str)
            + f" < {PUBLR_THRESHOLD}")
        if ml_prob is not None and ml_prob.notna().any():
            cuts = calibrated_cuts(classif, ml_prob, clean)
            for category in PANTRANSCRIPTOME_CATEGORIES:
                logger.info(f"Calibrated cut of {category}: "
                            f"{cuts.get(category, 'none, graded pubLR')}")
            cut = classif["structural_category"].map(cuts).astype(float)
            by_ml = scored & cut.notna() & ml_prob.notna()
            kept = (ml_prob >= cut) & clean
            result.loc[by_ml, "filter_route"] = ROUTE_ML
            result.loc[by_ml, "keep"] = kept[by_ml]
            result.loc[by_ml, "reason"] = ""
            result.loc[by_ml & ~clean, "reason"] = "Intra-priming"
            below = by_ml & clean & ~kept
            result.loc[below, "reason"] = (
                "POS_MLprob: " + ml_prob[below].round(3).astype(str)
                + " < " + cut[below].round(3).astype(str))
        elif ml_prob is not None:
            logger.warning("The random forest was not trained; the "
                           "NIC/NNC are filtered with graded pubLR.")
    result["keep"] = result["keep"].astype(bool)

    return result


def auto_filter(sqanti_class: str, json_file: str, force_multi_exon: bool,
                intrapriming: float, ml_prob: pd.Series | None,
                prefix: str, logger: logging.Logger) -> None:
    """Run the automatic filter and write its results.

    Args:
        sqanti_class: SQANTI3 classification file.
        json_file: JSON rules of the rules filter.
        force_multi_exon: Whether every mono-exon isoform is an
            artifact.
        intrapriming: Percentage of A downstream of the TTS from which
            an isoform is intra-primed.
        ml_prob: POS_MLprob of each isoform id, or None when the random
            forest was not run.
        prefix: Output prefix (directory and name). Writes
            <prefix>_AutoFilter_classification.txt (the classification
            with POS_MLprob, filter_route and filter_result),
            <prefix>_pass_isoforms.txt and
            <prefix>_filtering_reasons.txt.
        logger: Logger for the progress, routes and cuts.
    """
    classif = pd.read_csv(sqanti_class, sep="\t", dtype={"chrom": str})
    if not has_publr(classif):
        logger.warning("The classification has no pubLR probabilities "
                       "(QC was run without --pantranscriptome, or there "
                       "are no multi-exon NIC/NNC): applying the rules "
                       "filter to every isoform.")
    rules_dict = read_json_rules(json_file)
    keep_rules = apply_rules_table(classif, force_multi_exon,
                                   rules_dict) == "Isoform"
    if ml_prob is not None:
        ml_prob = classif["isoform"].map(ml_prob).astype(float)
        classif["POS_MLprob"] = ml_prob
    decisions = auto_decisions(classif, keep_rules, ml_prob, intrapriming,
                               logger)
    classif["filter_route"] = decisions["filter_route"]
    classif["filter_result"] = np.where(decisions["keep"], "Isoform",
                                        "Artifact")
    for route, group in classif.groupby("filter_route"):
        kept = (group["filter_result"] == "Isoform").sum()
        logger.info(f"Route {route}: {len(group)} isoforms, {kept} kept")
    artifacts = classif[classif["filter_result"] == "Artifact"]
    reasons = get_reasons_table(artifacts, force_multi_exon, rules_dict)
    own = decisions.loc[artifacts.index, "reason"]
    reasons.loc[own != "", "filter_reason"] = own[own != ""]
    classif.to_csv(f"{prefix}_AutoFilter_classification.txt", sep="\t",
                   index=False)
    classif.loc[classif["filter_result"] == "Isoform", "isoform"].to_csv(
        f"{prefix}_pass_isoforms.txt", sep="\t", index=False, header=False)
    reasons.to_csv(f"{prefix}_filtering_reasons.txt", sep="\t",
                   index=False)


def read_ml_probability(ml_classification: str) -> pd.Series:
    """Read the random forest probability written by the ML filter.

    Args:
        ml_classification: <prefix>_ML_classification.txt of the ML
            filter.

    Returns:
        POS_MLprob indexed by isoform (NaN where not scored, and
        everywhere when the random forest was not trained).

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not os.path.isfile(ml_classification):
        raise FileNotFoundError(ml_classification)
    table = pd.read_csv(ml_classification, sep="\t",
                        usecols=["isoform", "POS_MLprob"])
    result = pd.to_numeric(table.set_index("isoform")["POS_MLprob"],
                           errors="coerce")

    return result
