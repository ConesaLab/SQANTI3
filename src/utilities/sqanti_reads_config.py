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
        # bp of the annotated gene's end. Lower => more truncation/degradation. These
        # are absolute thresholds; the cohort-relative scorecard (below) is the
        # primary, dataset-agnostic view — leave these lenient by default.
        "perc_5p_within_window": {"warn": 40.0, "fail": 20.0, "worse": "low",
                                   "label": "Reads with complete 5' end (%)"},
        "perc_3p_within_window": {"warn": 40.0, "fail": 20.0, "worse": "low",
                                   "label": "Reads with complete 3' end (%)"},
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
