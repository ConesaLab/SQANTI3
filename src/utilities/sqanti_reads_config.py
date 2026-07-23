"""Configurable thresholds for SQANTI-reads.

Centralizes the numeric cut-offs that were previously hardcoded across the
plotting/report code. All defaults reproduce the historical behaviour, so a run
with no ``--config`` file is unchanged. Users can override any subset via a YAML
file passed to ``sqanti3_reads.py --config my_config.yaml``.
"""
import copy

# Default configuration — every value matches the previous hardcoded behaviour.
DEFAULT_CONFIG = {
    # Per-sample QC flag thresholds (drive the pass/warn/fail panel + qc_summary.json).
    # 'worse': 'high' => higher is worse (warn/fail are lower bounds);
    # 'worse': 'low'  => lower is worse (warn/fail are upper bounds).
    "qc_flags": {
        "perc_reads_intrapriming": {"warn": 20.0, "fail": 40.0, "worse": "high",
                                     "label": "Intra-priming reads (%)"},
        "perc_reads_RTS": {"warn": 15.0, "fail": 30.0, "worse": "high",
                            "label": "RT-switching reads (%)"},
        "perc_reads_non-canonical": {"warn": 10.0, "fail": 20.0, "worse": "high",
                                      "label": "Non-canonical reads (%)"},
        "median_length": {"warn": 500.0, "fail": 250.0, "worse": "low",
                           "label": "Median read length (bp)"},
        # Splice-site imprecision ("fuzziness"): % of donor+acceptor observations
        # that fall near a reference site (within +/- jxn_offset_window) but are
        # off by >0 bp. Higher => noisier splice-site boundaries.
        "perc_sites_imprecise": {"warn": 10.0, "fail": 20.0, "worse": "high",
                                  "label": "Imprecise splice sites (%)"},
        # 5'/3' completeness: % of reads whose end lands within completeness_window
        # bp of the annotated gene's end. A low value means reads are systematically
        # short at that end — consistent with truncation/degradation, but also with
        # alternative TSS/TTS or incomplete annotation, so it is a comparative signal,
        # not a degradation verdict. These absolute thresholds are deliberately lenient;
        # the cohort-relative scorecard (below) is the primary, dataset-agnostic view.
        "perc_5p_within_window": {"warn": 40.0, "fail": 20.0, "worse": "low",
                                   "label": "Reads with complete 5' end (%)"},
        "perc_3p_within_window": {"warn": 40.0, "fail": 20.0, "worse": "low",
                                   "label": "Reads with complete 3' end (%)"},
        # Novel non-canonical junction burden: % of all junctions that are novel
        # AND non-canonical. This class is ENRICHED for alignment / junction-calling
        # artefacts, but not exclusively artefactual (non-canonical splicing occurs
        # biologically, e.g. minor-spliceosome introns). An elevated value relative to
        # peers is worth inspecting; it is a candidate signal, not a definitive verdict.
        "perc_novel_noncanonical_jxn": {"warn": 2.0, "fail": 5.0, "worse": "high",
                                         "label": "Novel non-canonical junctions (%)"},
    },
    # A read is counted as intra-primed when %A downstream of its TTS exceeds this.
    "intrapriming_perc_A_cutoff": 60.0,
    # PCA loadings heatmap keeps the leading PCs explaining at least this fraction.
    "pca_cumulative_variance": 0.85,
    # Splice-site fuzziness: half-width (bp) of the signed-offset histogram and of
    # the precision profile, and the window within which a junction site counts as
    # an observation of a reference site (offsets beyond it are treated as novel,
    # not imprecise reference matches).
    "jxn_offset_window": 15,
    # 5'/3' completeness: a read end is "complete" if it lands within this many bp
    # of the annotated gene start (5') / end (3'). Used by the completeness-profile
    # plots and the perc_5p/3p_within_window flags above.
    "completeness_window": 50,
    # Cohort-relative sample-outlier scorecard. Each metric is turned into a robust
    # z-score against the COHORT's own distribution (median / MAD across samples),
    # so a sample is flagged for diverging from its peers — not against any absolute
    # baseline. This keeps the flag meaningful on any dataset (any organism,
    # protocol, length regime), including the common case where all samples are
    # clean (then every |z| is small and nothing is flagged). Needs >= min_samples
    # samples to be meaningful; disabled below that.
    "sample_scorecard": {
        # |robust z| at/above which a sample-metric cell is warned / failed.
        "z_warn": 2.5,
        "z_fail": 3.5,
        # A sample's OVERALL scorecard flag trips when at least this many of its
        # metrics are individually warned/failed (concordance across independent
        # axes is far more trustworthy than any single metric).
        "min_metrics_warn": 2,
        "min_metrics_fail": 2,
        # Below this many samples, cohort-relative z-scores are noise; skip the
        # scorecard entirely rather than emit misleading flags.
        "min_samples": 4,
        # Metrics entering the scorecard, each with the direction that is "worse".
        # 'high' => a large positive z is bad; 'low' => a large negative z is bad.
        # Any metric absent for a run (e.g. imprecision when hashing was skipped)
        # is dropped automatically.
        "metrics": {
            "median_length": "low",
            "perc_ISM": "high",
            "perc_novel_junctions": "high",
            "perc_5p_within_window": "low",
            "perc_3p_within_window": "low",
            "perc_reads_RTS": "high",
            "perc_reads_intrapriming": "high",
            "perc_sites_imprecise": "high",
            # Fraction of ISM reads that are 3'/5'/internal fragments (vs mono-exon /
            # intron-retention ISMs). A high fragment fraction is NOT inherently bad
            # -- alternative TSS/TTS and genuinely shorter isoforms land here too --
            # so it carries no absolute threshold; it is scorecard-only, where the
            # signal is purely "this sample's ISMs are shaped differently from its
            # peers'", which for replicates of one condition is worth a look.
            "perc_ISM_fragments": "high",
            "perc_novel_noncanonical_jxn": "high",
            # Distinct junctions recovered per 1000 reads; depth-normalised so a low
            # value is low junction complexity (consistent with degradation or a
            # simpler transcriptome) rather than shallow sequencing. Scored
            # cohort-relatively, no absolute threshold.
            "jxn_per_1k_reads": "low",
            # Mean Jensen–Shannon distance of a sample's structural-category
            # composition to the rest of the cohort; scored cohort-relatively, no
            # absolute threshold (a whole-cohort composition is not inherently wrong).
            "composition_drift": "high",
            # B11/B12/B13 orthogonal-support scalars — the share of a sample's
            # transcript 5' ends inside a CAGE peak (perc_within_cage), 3' ends at a
            # polyA site (perc_within_polya), and junctions with short-read coverage
            # (perc_jxn_SR_supported). Cohort-relative only (NO absolute qc_flag): a
            # low value means this sample's ends/junctions are less independently
            # supported than its peers', worth a look but not a fixed pass/fail. Each
            # drops out automatically when the supporting assay wasn't supplied to
            # SQANTI3 (columns all-NaN), so the fixture — which has none — is unaffected.
            "perc_within_cage": "low",
            "perc_within_polya": "low",
            "perc_jxn_SR_supported": "low",
        },
    },
}


def _deep_merge(base, override):
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_merge(base[k], v)
        else:
            base[k] = v
    return base


def load_config(path=None):
    """Return the default config, with a user YAML file merged over it if given."""
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    if path:
        import yaml
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
        _deep_merge(cfg, user)
    return cfg
