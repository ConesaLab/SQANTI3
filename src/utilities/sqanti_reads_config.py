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
