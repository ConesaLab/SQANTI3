"""Evidence check of the reference transcripts reintroduced by rescue.

Rescue filters the reference transcriptome without the requisites on long-read
counts (FL, FL.<sample>, prevalence, prevalence_<group>), since reference
transcripts have no counts of their own. Their evidence is the counts of the
artifacts assigned to them, which may have been discarded precisely for lack of
reproducible detection. This module computes, sample by sample, the counts each
reference target would receive from its artifacts with the redistribution of
requantification, recomputes FL, prevalence and prevalence_<group> on them as QC
does, and evaluates the count requisites of the rules with the rules filter
(apply_rules). Targets that fail are not reintroduced:
their artifacts are reassigned to the best supported long-read hit, or left
unrescued so requantification sends their counts to the gene residual.
"""
import json
import sys

import numpy as np
import pandas as pd

from src.config import MIN_EXPRESSION
from src.module_logging import rescue_logger
from src.utils import is_expressed
from src.utilities.rescue.sq_requant import build_artifact_table, redistribute_counts_vectorized
from src.utilities.filter.sqanti3_rules_filter import read_json_rules, apply_rules

# Structural category whose rules apply to reference targets: in the reference
# classification every transcript is a full-splice match of itself.
TARGET_CATEGORY = "full-splice_match"

EVIDENCE_COLUMN = "evidence_check"

# Requisites evaluated on the aggregated counts: the long-read count columns
# written by QC (FL, FL.<sample>, prevalence and prevalence_<group>).
COUNT_RULES = ["FL", "prevalence"]


def is_count_rule(column: str) -> bool:
    """Whether a rule column is computed from long-read counts."""
    return any(column == name or column.startswith((f"{name}.", f"{name}_"))
               for name in COUNT_RULES)


def target_count_rules(json_data: dict, category: str = TARGET_CATEGORY) -> list:
    """Count requisites that apply to reference targets, keeping the rule structure.

    Rules of a category are evaluated with OR, so a rule without count requisites
    makes counts unnecessary for the whole category.

    Args:
        json_data (dict): rules JSON as loaded from the --json_filter file
        category (str): structural category whose rules apply to the targets

    Returns:
        list: one dict of count requisites per rule, or an empty list if the
        targets do not need count evidence
    """
    sc = category if category in json_data else "rest"
    rule_sets = []
    for rule_set in json_data.get(sc, []):
        kept = {col: value for col, value in rule_set.items() if is_count_rule(col)}
        if not kept:
            return []
        rule_sets.append(kept)
    return rule_sets


def count_samples(counts_df: pd.DataFrame) -> list:
    """Sample columns of a count matrix loaded by sq_requant.load_counts()."""
    return [c for c in counts_df.columns if c != "isoform"]


def target_counts(rescue_df: pd.DataFrame, classif_df: pd.DataFrame,
                  counts_df: pd.DataFrame) -> pd.DataFrame:
    """Counts each reference target would receive from its artifacts, per sample.

    Runs the redistribution of requantification on the whole rescue table, before
    any target is removed: each artifact is split among all its targets, long-read
    and reference, in proportion to the counts of each target in the sample (evenly
    if none has counts), and integer counts stay integer. Reference targets have no
    counts of their own, so they get nothing from an artifact tied with a long-read
    isoform that has counts. Removing the targets that fail can only increase the
    counts of the others, so a single pass is enough.

    Args:
        rescue_df (pd.DataFrame): rescue table (artifact, assigned_transcript, origin)
        classif_df (pd.DataFrame): filter classification (isoform, filter_result, associated_gene)
        counts_df (pd.DataFrame): count matrix of requantification (sq_requant.load_counts)

    Returns:
        pd.DataFrame: one row per reference target (index), one column per sample,
        named as the classification (FL.<sample>, or FL with a single sample)
    """
    samples = count_samples(counts_df)
    names = ["FL"] if len(samples) == 1 else [f"FL.{s}" for s in samples]
    refs = rescue_df.loc[rescue_df["origin"] == "reference", "assigned_transcript"].unique()
    if len(refs) == 0:
        return pd.DataFrame(columns=names, dtype=float)

    artifacts_df = build_artifact_table(rescue_df[["artifact", "assigned_transcript"]], classif_df)
    final = redistribute_counts_vectorized(artifacts_df, classif_df, counts_df).set_index("isoform")
    out = final.reindex(refs)[samples].fillna(0)
    out.columns = names
    return out


def add_count_columns(agg: pd.DataFrame, design=None,
                      min_expression: float = MIN_EXPRESSION) -> pd.DataFrame:
    """Add the count columns of QC computed on the aggregated counts.

    FL is the sum over samples, prevalence the number of samples where the target
    is expressed with the --min_expression of QC (is_expressed()) and, with a
    counts design, prevalence_<group> the same over the samples of each group.
    """
    agg = agg.copy()
    sample_cols = list(agg.columns)
    detected = is_expressed(agg[sample_cols], min_expression)
    if sample_cols != ["FL"]:
        agg["FL"] = agg[sample_cols].sum(axis=1)
    agg["prevalence"] = detected.sum(axis=1)
    for group, samples in (design or {}).items():
        agg[f"prevalence_{group}"] = detected[[f"FL.{s}" for s in samples]].sum(axis=1)
    return agg


def failed_targets(rescue_df: pd.DataFrame, classif_df: pd.DataFrame, counts_df: pd.DataFrame,
                   rules_dict: dict, design=None, min_expression: float = MIN_EXPRESSION) -> set:
    """Reference targets whose counts (target_counts()) fail the count requisites.

    Args:
        counts_df (pd.DataFrame): count matrix of requantification (sq_requant.load_counts)
        rules_dict (dict): count requisites for the targets, parsed by read_json_rules()
        design (dict): group -> samples, needed for prevalence_<group> requisites
        min_expression (float): --min_expression of QC (minimum count for expression)

    Returns:
        set: reference targets that do not pass
    """
    agg = target_counts(rescue_df, classif_df, counts_df)
    if agg.empty:
        return set()
    agg = add_count_columns(agg, design, min_expression)
    agg["isoform"] = agg.index
    agg["structural_category"] = TARGET_CATEGORY
    agg["exons"] = np.nan
    result = agg.apply(lambda row: apply_rules(row, False, rules_dict), axis=1)
    return set(agg.index[result == "Artifact"])


def best_lr_fallback(artifacts, hits_df: pd.DataFrame, lr_isoforms) -> pd.DataFrame:
    """Best mapping hit among supported long-read isoforms, per artifact.

    Ties are kept, as in select_best_hits(): requantification splits the counts.

    Returns:
        pd.DataFrame: artifact, assigned_transcript
    """
    cols = ["artifact", "assigned_transcript"]
    hits = hits_df[hits_df["rescue_candidate"].isin(artifacts) &
                   hits_df["mapping_hit"].isin(lr_isoforms)]
    if hits.empty:
        return pd.DataFrame(columns=cols)

    best = hits["alignment_score"] == hits.groupby("rescue_candidate")["alignment_score"].transform("max")
    return (hits.loc[best, ["rescue_candidate", "mapping_hit"]]
            .drop_duplicates()
            .rename(columns={"rescue_candidate": "artifact", "mapping_hit": "assigned_transcript"})
            .reset_index(drop=True))


def write_evidence_rules(rule_sets: list, out_file: str) -> dict:
    """Write the count requisites of the targets as a rules JSON and parse it."""
    with open(out_file, "w") as f:
        json.dump({TARGET_CATEGORY: rule_sets}, f, indent=4)
    return read_json_rules(out_file)


def check_design_for_rules(rule_sets: list, design) -> None:
    """Exit if the rules use prevalence_<group> columns that the design does not define."""
    groups = {col[len("prevalence_"):] for rs in rule_sets for col in rs
              if col.startswith("prevalence_")}
    if not groups:
        return
    if design is None:
        rescue_logger.error(f"The rules use per-group prevalence ({sorted(groups)}): the evidence check "
                            "needs the --counts_design file given to SQANTI3 QC.")
        sys.exit(1)
    missing = sorted(groups - set(design))
    if missing:
        rescue_logger.error(f"Groups used in the rules but not defined in --counts_design: {missing}")
        sys.exit(1)


def apply_evidence_check(rescue_df: pd.DataFrame, inclusion_list, classif_df: pd.DataFrame,
                         counts_df: pd.DataFrame, hits_df: pd.DataFrame, rules_dict: dict, design=None,
                         min_expression: float = MIN_EXPRESSION):
    """Drop reference targets without count evidence and reassign their artifacts.

    Rows of failed targets stay in the rescue table, marked as failed and not
    reintroduced, for traceability. An artifact left with no passing target is
    reassigned to its best supported long-read hit; without one, it is left
    unrescued and requantification sends its counts to the gene residual.

    Returns:
        tuple: (inclusion list without the failed targets, rescue_df with an
        evidence_check column), in the order returned by the other rescue steps
    """
    failed = failed_targets(rescue_df, classif_df, counts_df, rules_dict, design, min_expression)

    df = rescue_df.copy()
    df[EVIDENCE_COLUMN] = np.where(df["origin"] == "reference", "pass", "not_required")

    is_failed = df["origin"].eq("reference") & df["assigned_transcript"].isin(failed)
    df.loc[is_failed, EVIDENCE_COLUMN] = "failed"
    df.loc[is_failed, "reintroduced"] = "no"

    still_assigned = set(df.loc[~is_failed, "artifact"])
    orphans = (df.loc[is_failed & ~df["artifact"].isin(still_assigned), ["artifact", "rescue_mode"]]
               .drop_duplicates("artifact"))

    lr_isoforms = classif_df.loc[classif_df["filter_result"] == "Isoform", "isoform"]
    fallback = best_lr_fallback(orphans["artifact"], hits_df, lr_isoforms).merge(orphans, on="artifact")
    fallback["origin"] = "lr_defined"
    fallback["reintroduced"] = "no"
    fallback[EVIDENCE_COLUMN] = "reassigned"
    if not fallback.empty:
        df = pd.concat([df, fallback[df.columns]], ignore_index=True)

    inclusion = pd.Series(np.asarray(inclusion_list).ravel(), name="isoform")
    inclusion = inclusion[~inclusion.isin(failed)].reset_index(drop=True)

    n_ref = df.loc[df["origin"].eq("reference"), "assigned_transcript"].nunique()
    rescue_logger.info(f"Evidence check: {len(failed)} of {n_ref} reference targets failed the count requisites.")
    rescue_logger.info(f"Evidence check: {fallback['artifact'].nunique()} of {len(orphans)} artifacts of failed "
                       "targets reassigned to long-read isoforms; the rest go to the gene residual.")
    return inclusion, df


def rows_for_requant(rescue_df: pd.DataFrame) -> pd.DataFrame:
    """Rescue table rows whose counts requantification must move.

    Rows of failed targets are kept in the table for traceability only.
    """
    if EVIDENCE_COLUMN not in rescue_df.columns:
        return rescue_df
    return rescue_df[rescue_df[EVIDENCE_COLUMN] != "failed"]
