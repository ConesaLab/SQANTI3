#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:32:18 2024

@author: nkeil
"""


import sys
import pandas as pd
import argparse
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({
    'figure.figsize': (12, 8),
    'axes.titlesize': 22,
    'axes.labelsize': 16,
    'font.weight': 'bold'
})
import matplotlib
import textwrap
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import jensenshannon
from matplotlib.ticker import FixedLocator
from contextlib import nullcontext
# The interactive HTML report is built by src/utilities/sqanti_reads_report.py
# (Plotly, imported lazily in main()); the old pdf2image/poppler rasterization
# path has been retired.
import logging
from dataclasses import dataclass
from typing import Optional

try:
    from src.module_logging import reads_logger
except ImportError:
    # Add the project root to sys.path if running as a script
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    try:
        from src.module_logging import reads_logger
    except ImportError:
        # Fallback logger if module_logging is not available
        reads_logger = logging.getLogger('reads_plots')
        reads_logger.setLevel(logging.INFO)
        if not reads_logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
            reads_logger.addHandler(handler)


@dataclass
class ReadsPlotArgs:
    """Arguments for the reads plotting functions."""
    inREF: str
    inDESIGN: str
    OUT: str
    PREFIX: str
    inFACTOR: Optional[str] = None
    ANNOTEXP: int = 100
    JXNEXP: int = 10
    PERCCOV: int = 20
    PERCMAXJXN: int = 80
    FACTORLVL: Optional[str] = None
    ALLTABLES: bool = False
    PCATABLES: bool = False
    report: str = 'pdf'
    config: Optional[dict] = None
    jobs: int = 1

    def cfg(self):
        """Return the config dict, loading defaults lazily if none was set."""
        if self.config is None:
            from src.utilities.sqanti_reads_config import load_config
            self.config = load_config()
        return self.config



# --- Shared plotting constants -------------------------------------------------
# Colors/orders for the PDF report (render_report_pdf) and the HTML report
# (sqanti_reads_report.py), kept as module-level names so both resolve the same
# palettes and category orders.
category_color_palette = {
    "FSM": "#6BAED6",
    "ISM": "#FC8D59",
    "NIC": "#78C679",
    "NNC": "#EE6A50",
    "GENIC": "#969696",
    "AS": "#66C2A4",
    "FUS": "#ffc125",
    "INTER": "#e9967a",
    "GI": "#41B6C4",
}

subcat_color_palette = {
    "alternative_5end": '#02314d',
    "alternative_3end": '#0e5a87',
    "alternative_3end5end": '#7ccdfc',
    'reference_match': '#c4e1f2',
    "3prime_fragment": '#c4531d',
    "internal_fragment": '#e37744',
    "5prime_fragment": '#e0936e',
    "combination_of_known_junctions": '#014d02',
    "combination_of_known_splicesites": '#379637',
    "intron_retention": '#81eb82',
    "at_least_one_novel_splicesite": '#6ec091',  ## Changed this from Not comb. of annot. junctions
    "mono-exon_by_intron_retention": '#4aaa72',
    "At least 1 annot. don./accept.": '#32734d',
    "mono-exon": '#cec2d2',
    "multi-exon": '#876a91',
    "mono_in_multi": '#aec6cf',
}

# NNC (novel_not_in_catalog) sub-category palette — reds, matching the NNC
# structural-category color (#EE6A50). Kept separate from subcat_color_palette
# because NNC shares sub-category *names* with NIC (e.g. intron_retention) but
# should read as red, not green. NNC has exactly two sub-categories: a novel
# splice site, or intron retention (which overrides when present).
nnc_subcat_color_palette = {
    "at_least_one_novel_splicesite": '#a50f15',
    "intron_retention": '#fb6a4a',
}

cat_order = ["FSM", "ISM", "NIC", "NNC", "AS", "FUS", "GENIC", "GI", "INTER"]
cat_order_stacked = ["INTER", "GI", "GENIC", "FUS", "AS", "NNC", "NIC", "ISM", "FSM"]

# Aesthetic categorical palette (a green→magenta "good→bad" severity ramp).
aes_palette = {
    "green": "#15918A", "orange": "#F58A53", "yellow": "#FDC659",
    "blue": "#74CDF0", "purple": "#9F7BB8", "pink": "#FDA3D1", "magenta": "#EE446F",
}

# Junction classification (known/novel × canonical/non-canonical). Colors match
# SQANTI3's own canonical palette (myPalette[c(1,7,3,2)] in SQANTI3_report.R).
jxn_palette = {
    "known_canonical": "#6BAED6",       # SQANTI3 myPalette[1]
    "known_non_canonical": "#FFC125",   # myPalette[7] (R goldenrod1)
    "novel_canonical": "#78C679",       # myPalette[3]
    "novel_non_canonical": "#FC8D59",   # myPalette[2]
}

# Three-series plots — artefacts (RTS/intra-priming/non-canonical) and
# donor/acceptor CV (ref_match/cv_0/cv_gt_0), in column order (good→bad).
three_series_palette = [aes_palette["green"], aes_palette["orange"], aes_palette["magenta"]]

# Legend display names for the donor/acceptor position-consistency categories
# (data columns keep their raw names; only the legend labels are prettified).
cv_count_legend_labels = {"ref_match": "ref_match", "cv_0": "CV=0", "cv_gt_0": "CV>0"}
cv_perc_legend_labels = {"perc_ref_match": "ref_match", "perc_cv_0": "CV=0", "perc_cv_gt_0": "CV>0"}

# Read-count bins, most support (best) → least support (worst).
readcount_palette = {
    "100+ reads": aes_palette["green"],
    "50-100 reads": aes_palette["blue"],
    "11-50 reads": aes_palette["yellow"],
    "2-10 reads": aes_palette["orange"],
    "1 read": aes_palette["magenta"],
}

# Read-length bins (count columns + their `_perc` variants), shortest→longest,
# so longer reads (generally better) read green.
_length_bins = ["reads_lt_1kb", "reads_1kb_to_2kb", "reads_2kb_to_3kb", "reads_gt_3kb"]
_length_colors = [aes_palette["magenta"], aes_palette["orange"], aes_palette["yellow"], aes_palette["green"]]
length_palette = {b: c for b, c in zip(_length_bins, _length_colors)}
length_palette.update({b + "_perc": c for b, c in zip(_length_bins, _length_colors)})

# Under-annotation gene categories (good→bad); shared by the PDF and HTML reports.
underannot_palette = {
    "annotated_with_well_covered_FSM": aes_palette["green"],
    "annotated_with_low_coverage_FSM": aes_palette["blue"],
    "underannotated_with_candidate_transcript": aes_palette["orange"],
    "underannotated_no_candidate_transcripts": aes_palette["magenta"],
}

# Per-sample qualitative colors (Plotly "Dark24"), used for scatter/strip/violin
# and PCA so samples are colored the same way in the PDF and HTML.
sample_seq = [
    '#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100',
    '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D',
    '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA',
    '#6C4516', '#0D2A63', '#AF0038',
]

# Full SQANTI3 structural-category names -> the abbreviations used by
# category_color_palette (for UpSet stacked bars, etc.).
SC_ABBR = {
    'full-splice_match': 'FSM', 'incomplete-splice_match': 'ISM',
    'novel_in_catalog': 'NIC', 'novel_not_in_catalog': 'NNC',
    'genic': 'GENIC', 'antisense': 'AS', 'fusion': 'FUS',
    'intergenic': 'INTER', 'genic_intron': 'GI',
}
# ------------------------------------------------------------------------------

def _palette_colors(columns, palette, default="#999999"):
    """Ordered colors for `columns` from a {name: hex} palette.

    Skips columns absent from the palette (e.g. a leftover 'sampleID' column), so
    the returned list matches the numeric series pandas actually plots.
    """
    return [palette[c] for c in columns if c in palette]


def _scatter_labeled(df, x_col, y_col, label_col, title, xlabel, ylabel,
                     factor_col=None, ax=None):
    """Scatter with per-point text labels, matching the HTML report.

    Points are a single color (no factor) or colored by factor level (with a
    factor-only legend); each point is annotated with its `label_col` value.
    This replaces seaborn scatterplots that keyed a busy hue+style legend to
    per-sample color and factor shapes.
    """
    import matplotlib.pyplot as _plt
    ax = ax or _plt.gca()
    faceted = _is_faceted(df, factor_col)
    if faceted:
        for i, lv in enumerate(pd.unique(df[factor_col])):
            d = df[df[factor_col] == lv]
            ax.scatter(d[x_col], d[y_col], s=120, color=sample_seq[i % len(sample_seq)],
                       label=str(lv))
        ax.legend(title=factor_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.scatter(df[x_col], df[y_col], s=120, color="#1f6feb")
    for _, row in df.iterrows():
        ax.annotate(str(row[label_col]), (row[x_col], row[y_col]),
                    textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def _is_faceted(df, factor_col):
    """True when a real (non-synthetic, multi-level) factor should split panels.

    ``temp_factor`` is the single-level placeholder the pipeline adds when the
    user gives no ``--factor``; it never faces.
    """
    return (bool(factor_col) and factor_col in df.columns
            and factor_col != "temp_factor" and df[factor_col].nunique() > 1)


def _stacked_bars_dict(*args, color_palette, exp_factor, **kwargs):
    """FacetGrid ``map_dataframe`` callback: stack every non-id column, colored
    by name from ``color_palette`` (a {column: hex} dict). Columns absent from
    the palette are skipped with a warning; column order sets the stack order."""
    data = kwargs.pop('data')
    kwargs.pop('color', None)
    ax = plt.gca()
    bottom = np.zeros(len(data))
    categories = [col for col in data.columns if col not in ['sampleID', exp_factor]]
    for category in categories:
        values = data[category].values
        if category in color_palette:
            non_zero_indices = values != 0
            ax.bar(data['sampleID'][non_zero_indices], values[non_zero_indices],
                   bottom=bottom[non_zero_indices], color=color_palette[category],
                   label=category, **kwargs)
            bottom += values
        else:
            reads_logger.warning(f"Color for {category} not found in palette.")


def _stacked_bars_indexed(*args, categories, palette, legend_labels=None, **kwargs):
    """FacetGrid ``map_dataframe`` callback: stack a fixed ``categories`` list in
    order, colored by the aligned ``palette`` list. A category missing from a
    facet subset is treated as all-zero so stack order/colors stay stable.
    ``legend_labels`` optionally maps a column name to its legend display name."""
    data = kwargs.pop('data')
    labels = legend_labels or {}
    ax = plt.gca()
    bottom = np.zeros(len(data))
    for idx, category in enumerate(categories):
        values = data[category].values if category in data.columns else np.zeros(len(data))
        non_zero_indices = values != 0
        ax.bar(data['sampleID'][non_zero_indices], values[non_zero_indices],
               bottom=bottom[non_zero_indices], color=palette[idx],
               label=labels.get(category, category))
        bottom += values


def _render_stacked_bar(pdf, df, *, exp_factor, num_factors, title, xlabel, ylabel,
                        palette, categories=None, legend_title=None,
                        height=8, aspect=1.3, sort=True, legend_labels=None):
    """Render one faceted stacked-bar page (one panel per ``exp_factor`` level).

    ``palette`` is either a {column: hex} dict (columns inferred + stacked in
    column order) or an ordered list aligned to ``categories`` (fixed stack).
    ``legend_labels`` optionally maps a category column to its legend display name.
    With the synthetic single-level ``temp_factor`` this collapses to one panel.
    """
    plot_df = df.sort_values(by=[exp_factor, 'sampleID']) if sort else df
    g = sns.FacetGrid(plot_df, col=exp_factor, col_wrap=num_factors, height=height,
                      aspect=aspect, sharex=False, sharey=True)
    if isinstance(palette, dict):
        g.map_dataframe(_stacked_bars_dict, color_palette=palette, exp_factor=exp_factor)
    else:
        g.map_dataframe(_stacked_bars_indexed, categories=categories, palette=palette,
                        legend_labels=legend_labels)
    for ax, (_name, group) in zip(g.axes.flatten(), plot_df.groupby(exp_factor)):
        # Pin tick positions to the sample count before labelling. Using the
        # actual sample positions (not ax.get_xticks()) keeps labels aligned even
        # when a panel drew no bars — e.g. an all-zero facet — which would
        # otherwise mismatch the default auto-ticks and raise a FixedLocator error.
        labels = list(group['sampleID'].unique())
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=90)
    g.set_axis_labels(xlabel, ylabel)
    g.set_titles("" if exp_factor == 'temp_factor' else exp_factor + " = {col_name}")
    title_obj = g.fig.suptitle(title, y=1.02, fontsize=20)
    legend_kwargs = {'bbox_to_anchor': (1.05, 1), 'loc': 'upper left'}
    if legend_title is not None:
        legend_kwargs['title'] = legend_title
    lgd = plt.legend(**legend_kwargs)
    plt.tight_layout()
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.subplots_adjust(top=0.85, right=0.8)
    pdf.savefig(bbox_extra_artists=(lgd, title_obj), bbox_inches='tight')
    plt.close()


def _grouped_bars_indexed(*args, categories, palette, **kwargs):
    """FacetGrid ``map_dataframe`` callback: side-by-side (grouped) bars, one
    group of ``categories`` per sample, colored by the aligned ``palette`` list."""
    data = kwargs.pop('data')
    ax = plt.gca()
    bar_width = 0.35
    num_samples = len(data['sampleID'].unique())
    gap_width = 0.1  # gap between per-sample groups
    total_bar_width = (len(categories) * bar_width) + gap_width
    positions = np.arange(num_samples) * total_bar_width
    for idx, category in enumerate(categories):
        category_positions = positions + idx * bar_width
        values = data[category].values if category in data.columns else np.zeros(len(data))
        non_zero_indices = values != 0
        ax.bar(category_positions[non_zero_indices], values[non_zero_indices],
               width=bar_width, color=palette[idx], label=category)
    central_offset = bar_width * (len(categories) - 1) / 2
    ax.set_xticks(positions + central_offset)
    ax.set_xticklabels(data['sampleID'].unique(), rotation=90)


def _render_grouped_bar(pdf, df, *, exp_factor, num_factors, categories, palette,
                        title, xlabel, ylabel, legend_title=None,
                        height=8, aspect=1.3, sort=True):
    """Render one faceted grouped-bar page (one panel per ``exp_factor`` level)."""
    plot_df = df.sort_values(by=[exp_factor, 'sampleID']) if sort else df
    g = sns.FacetGrid(plot_df, col=exp_factor, col_wrap=num_factors, height=height,
                      aspect=aspect, sharex=False, sharey=True)
    g.map_dataframe(_grouped_bars_indexed, categories=categories, palette=palette)
    for ax, (_name, group) in zip(g.axes.flatten(), plot_df.groupby(exp_factor)):
        ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
    g.set_axis_labels(xlabel, ylabel)
    g.set_titles("" if exp_factor == 'temp_factor' else exp_factor + " = {col_name}")
    title_obj = g.fig.suptitle(title, y=1.02, fontsize=20)
    legend_kwargs = {'bbox_to_anchor': (1.05, 1), 'loc': 'upper left'}
    if legend_title is not None:
        legend_kwargs['title'] = legend_title
    lgd = plt.legend(**legend_kwargs)
    plt.tight_layout()
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.subplots_adjust(top=0.85, right=0.8)
    pdf.savefig(bbox_extra_artists=(lgd, title_obj), bbox_inches='tight')
    plt.close()


def _vectorize_colorbars(fig):
    """Force heatmap/colorbar artists to render as vector paths instead of a
    raster image, so every mark on the page stays editable in vector tools like
    Illustrator. matplotlib rasterizes colorbar solids by default; the two
    heatmap pages (PCA loadings, replicate concordance) are the only report
    pages that would otherwise carry a raster element."""
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(False)
        for im in ax.images:
            im.set_rasterized(False)


# Column index of each flaggable metric in the summary table below.
_SUMMARY_FLAG_COL = {
    'median_length': 3,
    'perc_reads_intrapriming': 6,
    'perc_reads_RTS': 7,
    'perc_reads_non-canonical': 8,
    'perc_5p_within_window': 9,
    'perc_3p_within_window': 10,
    'perc_sites_imprecise': 11,
}
_FLAG_BG = {"warn": "#FFE0B2", "fail": "#FFCDD2"}
_FLAG_FG = {"warn": "#B45309", "fail": "#B71C1C"}
_OVERALL_BG = {"pass": "#C8E6C9", "warn": "#FFE0B2", "fail": "#FFCDD2"}


def _render_summary_table_page(pdf, per_sample, samples, thresholds):
    """Render the per-sample summary-metrics + QC-flag table as one PDF page.

    Mirrors the HTML report's top table: cells that trigger a warn/fail flag are
    bold and tinted with the flag color, and the Overall cell is flag-colored."""
    headers = ["Sample", "Reads", "Genes", "Median len", "% >1kb", "% FSM",
               "% intra-prim", "% RTS", "% non-canon", "% 5' compl", "% 3' compl",
               "% fuzzy sites", "Overall"]
    rows = []
    for s in samples:
        m = per_sample[s]
        rows.append([
            str(s), f'{m["total_reads"]:,}', f'{m["genes_detected"]}',
            f'{m["median_length"]:.0f}', f'{m["perc_reads_gt_1kb"]:.1f}',
            (f'{m["perc_FSM"]:.1f}' if m["perc_FSM"] is not None else "—"),
            f'{m["perc_reads_intrapriming"]:.1f}',
            f'{m["perc_reads_RTS"]:.1f}',
            f'{m["perc_reads_non-canonical"]:.1f}',
            (f'{m["perc_5p_within_window"]:.1f}' if m.get("perc_5p_within_window") is not None else "—"),
            (f'{m["perc_3p_within_window"]:.1f}' if m.get("perc_3p_within_window") is not None else "—"),
            (f'{m["perc_sites_imprecise"]:.1f}' if m.get("perc_sites_imprecise") is not None else "—"),
            m["overall_flag"].upper(),
        ])

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.axis('off')
    ax.set_title("Summary metrics & QC flags", fontsize=20, pad=24)
    tbl = ax.table(cellText=rows, colLabels=headers, cellLoc='center', loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1, 1.8)

    for j in range(len(headers)):  # header row
        cell = tbl[0, j]
        cell.set_text_props(weight='bold', color='white')
        cell.set_facecolor(aes_palette["green"])

    for i, s in enumerate(samples, start=1):
        m = per_sample[s]
        for metric, col in _SUMMARY_FLAG_COL.items():
            fl = m["flags"].get(metric, "pass")
            if fl in ("warn", "fail"):
                cell = tbl[i, col]
                cell.set_text_props(weight='bold', color=_FLAG_FG[fl])
                cell.set_facecolor(_FLAG_BG[fl])
        ov = tbl[i, len(headers) - 1]
        ov.set_text_props(weight='bold')
        ov.set_facecolor(_OVERALL_BG.get(m["overall_flag"], "#ffffff"))

    matplotlib.rcParams['pdf.fonttype'] = 42
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def compute_ujc_metrics(ujc_count_DF, factor_col=None, n_depths=25):
    """Derive UJC-level QC metrics shared by the PDF and HTML reports.

    Returns a dict with:
      - 'samples': sample order
      - 'sample_factor': {sampleID: factor level} (None-safe)
      - 'saturation': long DataFrame [sampleID, depth, unique_ujcs] (rarefaction)
      - 'concordance': square DataFrame of per-sample UJC read-count correlations
      - 'upset': DataFrame indexed by jxnHash with one boolean column per sample
                 plus a 'structural_category' column (for the UpSet plot)
    """
    from scipy.special import gammaln

    fac = factor_col if (factor_col and factor_col in ujc_count_DF.columns) else None
    df = ujc_count_DF[["jxnHash", "read_count", "structural_category", "sampleID"]].copy()
    samples = list(pd.unique(df["sampleID"]))
    if fac:
        sample_factor = (ujc_count_DF.drop_duplicates("sampleID")
                         .set_index("sampleID")[fac].astype(str).to_dict())
    else:
        sample_factor = {s: None for s in samples}

    # --- Saturation / rarefaction: expected unique UJCs vs subsampling depth ---
    # `counts` is the read-count vector of the UJCs of interest; `total` is the
    # sample's TOTAL read count. For a category subset, summing over that subset
    # with the sample total gives the category's contribution to the overall
    # curve, on the same (total-depth) x-axis — so the categories decompose it.
    def _rarefy(counts, depths, total):
        counts = np.asarray(counts, dtype=float)
        S = len(counts)
        if total <= 0 or S == 0:
            return [0.0] * len(depths)
        out = []
        for n in depths:
            if n <= 0:
                out.append(0.0)
            elif n >= total:
                out.append(float(S))
            else:
                log_c_N_n = gammaln(total + 1) - gammaln(n + 1) - gammaln(total - n + 1)
                nc = total - counts
                ratio = np.zeros(S)
                ok = nc >= n
                ratio[ok] = np.exp((gammaln(nc[ok] + 1) - gammaln(n + 1)
                                    - gammaln(nc[ok] - n + 1)) - log_c_N_n)
                out.append(float((1.0 - ratio).sum()))
        return out

    sat_rows = []
    sat_cat_rows = []
    for s in samples:
        sdf = df[df["sampleID"] == s]
        counts = sdf["read_count"].values
        total = int(counts.sum())
        depths = np.unique(np.linspace(0, total, n_depths).astype(int))
        for d, r in zip(depths, _rarefy(counts, depths, total)):
            sat_rows.append({"sampleID": s, "depth": int(d), "unique_ujcs": r})
        # per-structural-category saturation (same total-depth x-axis)
        for cat, cdf in sdf.groupby("structural_category"):
            abbr = SC_ABBR.get(cat, str(cat))
            cc = cdf["read_count"].values
            for d, r in zip(depths, _rarefy(cc, depths, total)):
                sat_cat_rows.append({"sampleID": s, "structural_category": abbr,
                                     "depth": int(d), "unique_ujcs": r})
    saturation = pd.DataFrame(sat_rows)
    saturation_by_category = pd.DataFrame(sat_cat_rows)

    # --- Per-sample UJC read-count matrix (jxnHash x sample) ---
    mat = df.pivot_table(index="jxnHash", columns="sampleID", values="read_count",
                         aggfunc="sum", fill_value=0).reindex(columns=samples)

    # --- Replicate concordance: Pearson correlation of the count vectors ---
    concordance = mat.corr(method="pearson")

    # --- UpSet membership matrix + structural category per UJC ---
    upset = (mat > 0)
    sc = df.groupby("jxnHash")["structural_category"].agg(
        lambda x: x.value_counts().index[0])
    upset = upset.assign(structural_category=sc.reindex(upset.index))

    return {"samples": samples, "sample_factor": sample_factor,
            "saturation": saturation, "saturation_by_category": saturation_by_category,
            "concordance": concordance, "upset": upset}


def compute_ujc_overlap(ujc_metrics):
    """A6 — pairwise UJC overlap index between samples.

    From the UJC presence matrix (``ujc_metrics['upset']``, one boolean column
    per sample), build an asymmetric matrix whose entry (row A, col B) is
    ``|A ∩ B| / |A|``: the fraction of sample A's unique junction chains that are
    also detected in sample B. The diagonal is 1. Reading across row A shows how
    much of A's UJC repertoire each other sample recovers; the matrix is
    asymmetric (A→B ≠ B→A) whenever the two repertoires differ in size, which is
    exactly what distinguishes a small deeply-shared sample from a large one.
    This is a between-sample comparative view (like the read-count concordance
    and UpSet), not a per-sample quality score, so it feeds no scorecard metric.
    Empty-safe: returns ``{'samples': [], 'matrix': None}`` with fewer than two
    samples that carry UJCs.
    """
    if not ujc_metrics:
        return {"samples": [], "matrix": None}
    upset = ujc_metrics.get("upset")
    samples = [s for s in (ujc_metrics.get("samples") or [])
               if upset is not None and s in upset.columns]
    if upset is None or len(samples) < 2:
        return {"samples": [], "matrix": None}
    present = {s: upset[s].to_numpy(dtype=bool) for s in samples}
    mat = pd.DataFrame(index=samples, columns=samples, dtype=float)
    for a in samples:
        na = int(present[a].sum())
        for b in samples:
            inter = int(np.logical_and(present[a], present[b]).sum())
            mat.loc[a, b] = (inter / na) if na > 0 else 0.0
    return {"samples": samples, "matrix": mat}


def compute_upset_intersections(upset_DF, samples, max_intersections=20):
    """Compute UpSet intersections of shared UJCs, broken down by structural category.

    Returns {'intersections': [{'combo': (samples...), 'total': n,
                                'sc_counts': {SC: n}} ...] sorted by size,
             'set_sizes': {sample: {SC: n}},
             'sc_order': [SC ...]}  — shared by the PDF and HTML renderers.
    """
    if upset_DF is None or upset_DF.empty or len(samples) < 2:
        return None
    df = upset_DF.copy()
    df["SC"] = df["structural_category"].map(lambda x: SC_ABBR.get(x, "other")).fillna("other")
    memb = df[list(samples)].astype(bool)
    df["_combo"] = memb.apply(lambda r: tuple(s for s in samples if r[s]), axis=1)
    df = df[df["_combo"].map(len) > 0]
    if df.empty:
        return None
    inter = [{"combo": combo, "total": int(len(g)),
              "sc_counts": g["SC"].value_counts().to_dict()}
             for combo, g in df.groupby("_combo")]
    inter.sort(key=lambda d: -d["total"])
    inter = inter[:max_intersections]
    set_sizes = {s: df[df["_combo"].map(lambda c: s in c)]["SC"].value_counts().to_dict()
                 for s in samples}
    present = set()
    for d in inter:
        present |= set(d["sc_counts"])
    sc_order = [c for c in cat_order if c in present] + sorted(present - set(cat_order))
    return {"intersections": inter, "set_sizes": set_sizes, "sc_order": sc_order}


def _plot_upset_pdf(pdf, upset_DF, samples):
    """Append a (custom) UpSet plot of shared UJCs, stacked by structural category."""
    import matplotlib.gridspec as gridspec
    up = compute_upset_intersections(upset_DF, samples)
    if not up:
        return
    inter, set_sizes, sc_order = up["intersections"], up["set_sizes"], up["sc_order"]
    colors = {c: category_color_palette.get(c, "#969696") for c in sc_order}
    n = len(inter)
    ns = len(samples)

    fig = plt.figure(figsize=(max(11, n * 0.7 + 4), 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2.2, 10], height_ratios=[3, 1.6],
                           hspace=0.06, wspace=0.06)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0])

    x = np.arange(n)
    # Intersection bars, stacked by structural category
    bottom = np.zeros(n)
    for c in sc_order:
        vals = np.array([d["sc_counts"].get(c, 0) for d in inter], dtype=float)
        ax_bar.bar(x, vals, bottom=bottom, color=colors[c], label=c, width=0.8)
        bottom += vals
    for xi, d in enumerate(inter):
        ax_bar.text(xi, d["total"], str(d["total"]), ha="center", va="bottom", fontsize=8)
    ax_bar.set_ylabel("# UJCs")
    ax_bar.set_title("Shared UJCs across samples (stacked by structural category)")
    ax_bar.legend(title="Structural category", bbox_to_anchor=(1.01, 1),
                  loc="upper left", fontsize=8)
    plt.setp(ax_bar.get_xticklabels(), visible=False)

    # Membership matrix (dots + connecting line)
    idx = {s: i for i, s in enumerate(samples)}
    for j, d in enumerate(inter):
        for si in range(ns):
            ax_mat.plot(j, si, "o", color="#dddddd", markersize=13, zorder=1)
        present = sorted(idx[s] for s in d["combo"])
        for si in present:
            ax_mat.plot(j, si, "o", color="#333333", markersize=13, zorder=3)
        if len(present) > 1:
            ax_mat.plot([j, j], [present[0], present[-1]], color="#333333", lw=2, zorder=2)
    ax_mat.set_xticks([])
    ax_mat.set_yticks(range(ns))
    ax_mat.set_yticklabels(samples)
    ax_mat.set_ylim(ns - 0.5, -0.5)  # sample 0 at top (no invert_yaxis)

    # Set-size bars (per sample), stacked by structural category
    for i, s in enumerate(samples):
        left = 0.0
        for c in sc_order:
            v = set_sizes.get(s, {}).get(c, 0)
            ax_set.barh(i, v, left=left, color=colors[c], height=0.6)
            left += v
    ax_set.set_xlabel("Set size")
    ax_set.set_ylim(ns - 0.5, -0.5)
    ax_set.set_yticks(range(ns))
    ax_set.set_yticklabels(samples)
    ax_set.invert_xaxis()

    matplotlib.rcParams['pdf.fonttype'] = 42
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def compute_jxn_offset_metrics(jxn_offset_DF, window=15):
    """Derive splice-site fuzziness views from the collapsed offset table.

    Input is the long [sampleID, site, offset, canonical, count] frame from
    ``_summarize_jxn_offsets`` (already restricted to +/- ``window`` bp). Returns a
    dict shared by the PDF and HTML renderers:

      - 'samples'   : sample order
      - 'window'    : the +/- bp half-width used
      - 'spectrum'  : [sampleID, site, offset, count] signed-offset histogram
                      (exact matches, offset==0, kept — the plot excludes them)
      - 'profile'   : [sampleID, k, perc_within] cumulative % of sites within +/-k
      - 'per_sample': [sampleID, n_sites, perc_exact, perc_imprecise,
                       median_abs_offset, perc_downstream] one row per sample
      - 'by_class'  : [sampleID, canonical, perc_imprecise, n] canonical split

    ``perc_imprecise`` (share of near-reference site observations off by >0 bp) is
    the scalar wired into the QC-flag panel / qc_summary.json.
    """
    empty = {"samples": [], "window": window,
             "spectrum": pd.DataFrame(columns=["sampleID", "site", "offset", "count"]),
             "profile": pd.DataFrame(columns=["sampleID", "k", "perc_within"]),
             "per_sample": pd.DataFrame(columns=["sampleID", "n_sites", "perc_exact",
                                                 "perc_imprecise", "median_abs_offset",
                                                 "perc_downstream"]),
             "by_class": pd.DataFrame(columns=["sampleID", "canonical", "perc_imprecise", "n"])}
    if jxn_offset_DF is None or jxn_offset_DF.empty:
        return empty

    df = jxn_offset_DF.copy()
    df["offset"] = df["offset"].astype(int)
    df["count"] = df["count"].astype(int)
    df["abs"] = df["offset"].abs()
    samples = list(pd.unique(df["sampleID"]))

    # Signed-offset spectrum per (sample, site).
    spectrum = (df.groupby(["sampleID", "site", "offset"], observed=True)["count"]
                  .sum().reset_index())

    # Cumulative precision profile: % of a sample's site observations within +/-k.
    prof_rows = []
    ks = range(0, window + 1)
    for s in samples:
        sd = df[df["sampleID"] == s]
        tot = int(sd["count"].sum())
        for k in ks:
            within = int(sd.loc[sd["abs"] <= k, "count"].sum())
            prof_rows.append({"sampleID": s, "k": k,
                              "perc_within": (within / tot * 100) if tot else 0.0})
    profile = pd.DataFrame(prof_rows)

    # Per-sample scalars.
    ps_rows = []
    for s in samples:
        sd = df[df["sampleID"] == s]
        tot = int(sd["count"].sum())
        exact = int(sd.loc[sd["offset"] == 0, "count"].sum())
        fuzzy = sd[sd["offset"] != 0]
        n_fuzzy = int(fuzzy["count"].sum())
        # weighted median of |offset| over fuzzy observations
        med = 0.0
        if n_fuzzy:
            f = fuzzy.sort_values("abs")
            cum = f["count"].cumsum()
            med = float(f.loc[cum >= n_fuzzy / 2, "abs"].iloc[0])
        downstream = int(fuzzy.loc[fuzzy["offset"] > 0, "count"].sum())
        ps_rows.append({
            "sampleID": s,
            "n_sites": tot,
            "perc_exact": (exact / tot * 100) if tot else 0.0,
            "perc_imprecise": (n_fuzzy / tot * 100) if tot else 0.0,
            "median_abs_offset": med,
            "perc_downstream": (downstream / n_fuzzy * 100) if n_fuzzy else 0.0,
        })
    per_sample = pd.DataFrame(ps_rows)

    # Canonical vs non-canonical fuzziness.
    bc_rows = []
    if "canonical" in df.columns:
        for (s, can), g in df.groupby(["sampleID", "canonical"], observed=True):
            tot = int(g["count"].sum())
            n_fuzzy = int(g.loc[g["offset"] != 0, "count"].sum())
            bc_rows.append({"sampleID": s, "canonical": str(can),
                            "perc_imprecise": (n_fuzzy / tot * 100) if tot else 0.0,
                            "n": tot})
    by_class = pd.DataFrame(bc_rows) if bc_rows else empty["by_class"]

    return {"samples": samples, "window": window, "spectrum": spectrum,
            "profile": profile, "per_sample": per_sample, "by_class": by_class}


def compute_completeness_metrics(completeness_DF, window=50):
    """Derive 5'/3' completeness views from the collapsed |distance| digest.

    Input is the long [sampleID, end, abs_dist, count] frame from
    ``_summarize_completeness``. Returns a dict shared by both renderers:

      - 'samples'   : sample order
      - 'window'    : the "complete within +/- window bp" cut used for scalars
      - 'profile'   : [sampleID, end, k, perc_within] cumulative % of reads whose
                      end lands within k bp, for a set of k values (the ECDF)
      - 'per_sample': [sampleID, perc_5p_within_window, perc_3p_within_window,
                       median_abs_5p, median_abs_3p]

    Empty-safe: an empty input yields empty frames and an empty sample list, so
    the plot/scorecard code no-ops (e.g. when diff_to_gene_* were not in the
    SQANTI3 output).
    """
    empty = {"samples": [], "window": window,
             "profile": pd.DataFrame(columns=["sampleID", "end", "k", "perc_within"]),
             "per_sample": pd.DataFrame(columns=["sampleID", "perc_5p_within_window",
                                                 "perc_3p_within_window",
                                                 "median_abs_5p", "median_abs_3p"])}
    if completeness_DF is None or completeness_DF.empty:
        return empty

    df = completeness_DF.copy()
    df["abs_dist"] = df["abs_dist"].astype(int)
    df["count"] = df["count"].astype(int)
    samples = list(pd.unique(df["sampleID"]))

    # k grid for the ECDF: dense near 0, sparser out to 10*window.
    ceiling = int(10 * window)
    ks = sorted(set(list(range(0, min(200, ceiling) + 1, 5)) + [window, ceiling]))

    def _wmedian(sub):
        tot = int(sub["count"].sum())
        if not tot:
            return float("nan")
        s = sub.sort_values("abs_dist")
        cum = s["count"].cumsum()
        return float(s.loc[cum >= tot / 2, "abs_dist"].iloc[0])

    prof_rows, ps_rows = [], []
    end_key = {"5prime": "perc_5p_within_window", "3prime": "perc_3p_within_window"}
    med_key = {"5prime": "median_abs_5p", "3prime": "median_abs_3p"}
    for s in samples:
        sd = df[df["sampleID"] == s]
        row = {"sampleID": s}
        for end in ("5prime", "3prime"):
            ed = sd[sd["end"] == end]
            tot = int(ed["count"].sum())
            for k in ks:
                within = int(ed.loc[ed["abs_dist"] <= k, "count"].sum())
                prof_rows.append({"sampleID": s, "end": end, "k": k,
                                  "perc_within": (within / tot * 100) if tot else 0.0})
            within_w = int(ed.loc[ed["abs_dist"] <= window, "count"].sum())
            row[end_key[end]] = (within_w / tot * 100) if tot else float("nan")
            row[med_key[end]] = _wmedian(ed)
        ps_rows.append(row)

    return {"samples": samples, "window": window,
            "profile": pd.DataFrame(prof_rows), "per_sample": pd.DataFrame(ps_rows)}


def compute_jxn_yield_vs_depth(jxn_counts, read_counts):
    """Depth-normalised junction yield per sample (A4).

    ``jxn_per_1k_reads = n_jxn / n_reads * 1000`` — junctions recovered per 1000
    reads, so a low value is genuinely lower junction complexity (consistent with
    degradation or a simpler transcriptome) rather than merely shallower
    sequencing. Cohort-relative; no absolute meaning. Empty-safe on missing counts.
    """
    empty = {"samples": [],
             "per_sample": pd.DataFrame(columns=["sampleID", "n_jxn", "n_reads",
                                                 "jxn_per_1k_reads"])}
    if not jxn_counts or not read_counts:
        return empty
    rows = []
    for s in jxn_counts:
        n_reads = read_counts.get(s)
        if n_reads is None or n_reads <= 0:
            continue
        n_jxn = jxn_counts[s]
        rows.append({"sampleID": str(s), "n_jxn": int(n_jxn), "n_reads": int(n_reads),
                     "jxn_per_1k_reads": n_jxn / n_reads * 1000.0})
    if not rows:
        return empty
    return {"samples": [r["sampleID"] for r in rows],
            "per_sample": pd.DataFrame(rows)}


def compute_composition_drift(all_gene_percs_pivot_DF):
    """Between-sample composition drift via pairwise Jensen–Shannon distance (A3).

    Each sample's structural-category percentages become a probability vector over
    the shared category set; pairwise Jensen–Shannon distance (base-2, symmetric,
    bounded [0,1]) measures how different two samples' category mixes are, reduced
    to each sample's mean distance to the rest of the cohort. High = an atypical
    composition relative to peers — worth checking whether it is a genuine
    biological difference (e.g. a distinct condition) or a technical one; the value
    itself is not a pass/fail. Empty-safe (<2 samples or no categories).
    """
    empty = {"samples": [], "matrix": pd.DataFrame(),
             "per_sample": pd.DataFrame(columns=["sampleID", "mean_jsd"])}
    df = all_gene_percs_pivot_DF
    if df is None or "sampleID" not in getattr(df, "columns", []):
        return empty
    cat_cols = [c for c in cat_order if c in df.columns]
    if len(df) < 2 or not cat_cols:
        return empty
    samples = df["sampleID"].astype(str).tolist()
    M = df[cat_cols].fillna(0).to_numpy(dtype=float)
    # Drop categories that are zero across the whole cohort, then row-normalise.
    keep = M.sum(axis=0) > 0
    M = M[:, keep]
    if M.shape[1] == 0:
        return empty
    row_tot = M.sum(axis=1, keepdims=True)
    P = np.divide(M, row_tot, out=np.zeros_like(M), where=row_tot > 0)
    n = len(samples)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = jensenshannon(P[i], P[j], base=2)
            d = 0.0 if not np.isfinite(d) else float(d)
            D[i, j] = D[j, i] = d
    matrix = pd.DataFrame(D, index=samples, columns=samples)
    mean_jsd = D.sum(axis=1) / max(n - 1, 1)
    per_sample = pd.DataFrame({"sampleID": samples, "mean_jsd": mean_jsd})
    return {"samples": samples, "matrix": matrix, "per_sample": per_sample}


def compute_tandem_sites(jxn_offset_metrics, tandem_offsets=(3, -3, 4, -4, 6, -6)):
    """Tandem splice-site (NAGNAG) detector (F6).

    Excess mass at ±3 bp (and ±4/±6) in the signed splice-site offset spectrum is
    the signature of tandem donor/acceptor usage — reads legitimately choosing a
    nearby alternative site, a real biological phenomenon distinct from random
    boundary imprecision. Reports, per sample, the % of imprecise (offset!=0)
    observations sitting at the tandem offsets. Descriptive (no scorecard, no
    threshold). Empty-safe. Reads the existing offset spectrum, not new columns.
    """
    empty = {"samples": [],
             "per_sample": pd.DataFrame(columns=["sampleID", "perc_tandem", "perc_tandem_3bp"]),
             "by_offset": pd.DataFrame(columns=["sampleID", "offset", "count"])}
    m = jxn_offset_metrics
    if not m or not m.get("samples") or m.get("spectrum") is None or m["spectrum"].empty:
        return empty
    sp = m["spectrum"]
    sp = sp[sp["offset"] != 0]   # imprecise observations only
    if sp.empty:
        return empty
    tand3 = {3, -3}
    tandall = set(tandem_offsets)
    rows = []
    for s, g in sp.groupby("sampleID"):
        tot = float(g["count"].sum())
        c3 = float(g.loc[g["offset"].isin(tand3), "count"].sum())
        call = float(g.loc[g["offset"].isin(tandall), "count"].sum())
        rows.append({"sampleID": str(s),
                     "perc_tandem_3bp": (c3 / tot * 100) if tot else 0.0,
                     "perc_tandem": (call / tot * 100) if tot else 0.0})
    near = (sp[sp["offset"].abs() <= 10]
            .groupby(["sampleID", "offset"], as_index=False)["count"].sum())
    return {"samples": [r["sampleID"] for r in rows],
            "per_sample": pd.DataFrame(rows), "by_offset": near}


_A5_AXES = ["category_composition", "length_profile", "imprecision"]
# A8 within-group robust-spread (MAD) companion columns, aligned to _A5_AXES.
_A5_MAD_COLS = ["mad_category_composition", "mad_length_profile", "mad_imprecision"]


def _rel_agreement(x, y):
    """1 - mean relative difference between two vectors, in [0, 1] (1 = identical)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    denom = np.abs(x) + np.abs(y)
    rel = np.divide(np.abs(x - y), denom, out=np.zeros_like(denom), where=denom > 0)
    return float(1.0 - np.mean(rel))


def _mad(values):
    """Median absolute deviation (robust spread): median(|x − median(x)|).

    0 when the replicates are identical on this axis; larger = more within-group
    spread. NaN for fewer than two finite values. Not scaled by 1.4826 — reported
    as a plain typical-deviation in the axis's own units."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if v.size < 2:
        return np.nan
    return float(np.median(np.abs(v - np.median(v))))


def compute_replicate_concordance(all_gene_percs_pivot_DF, length_DF,
                                  jxn_offset_metrics, factor_map):
    """Multi-axis within-group replicate concordance (A5).

    For samples sharing a design-factor level (replicates), how closely each
    agrees with its group-mates on: structural-category composition (1 − mean
    Jensen–Shannon distance), read-length profile (q25/median/q75), and
    splice-site imprecision. Each axis is a 0–1 agreement (1 = matches its
    replicates). A replicate low on one axis diverges from its peers on that
    specific property — more sensitive than UJC overlap alone. Only meaningful
    with a factor and ≥2 replicates in a level; otherwise disabled with a reason.
    """
    def disabled(reason):
        return {"enabled": False, "reason": reason, "samples": [], "groups": {},
                "axes": _A5_AXES, "mad_axes": _A5_MAD_COLS,
                "per_sample": pd.DataFrame(
                    columns=["sampleID", "group"] + _A5_AXES + _A5_MAD_COLS)}

    if not factor_map:
        return disabled("no design factor supplied")
    groups = {}
    for s, lv in factor_map.items():
        groups.setdefault(str(lv), []).append(str(s))
    multi = {lv: ss for lv, ss in groups.items() if len(ss) >= 2}
    if not multi:
        return disabled("no factor level has ≥2 replicates")

    # Per-sample vectors from tables already computed.
    comp = {}
    if all_gene_percs_pivot_DF is not None and "sampleID" in all_gene_percs_pivot_DF.columns:
        cat_cols = [c for c in cat_order if c in all_gene_percs_pivot_DF.columns]
        piv = all_gene_percs_pivot_DF.set_index("sampleID")
        for s in piv.index:
            v = piv.loc[s, cat_cols].fillna(0).to_numpy(dtype=float) if cat_cols else np.array([])
            tot = v.sum()
            comp[str(s)] = v / tot if tot > 0 else v
    lenv = {}
    lcols = ["q25_length", "median_length", "q75_length"]
    if length_DF is not None and all(c in length_DF.columns for c in lcols):
        li = length_DF.set_index("sampleID")
        for s in li.index:
            lenv[str(s)] = li.loc[s, lcols].to_numpy(dtype=float)
    impr = {}
    if jxn_offset_metrics and not jxn_offset_metrics.get("per_sample", pd.DataFrame()).empty:
        impr = {str(k): float(v) for k, v in
                jxn_offset_metrics["per_sample"].set_index("sampleID")["perc_imprecise"].items()}

    # A8 — per-group robust within-group spread (MAD) on each axis, a companion to
    # the agreement scores above: the agreement says how close each replicate is to
    # its group-mates, the MAD says how tightly the whole group clusters (and, unlike
    # pairwise agreement, is non-degenerate for exactly two replicates). Each axis is
    # reduced to one scalar per sample, then MAD'd across the group (same value for
    # every member; reported in the axis's own units):
    #   composition -> Jensen–Shannon distance of the sample's category mix to the
    #                  group-mean mix (JSD, 0–1); length -> median read length (bp);
    #                  imprecision -> % imprecise splice sites (percentage points).
    group_mad = {}
    for lv, members in multi.items():
        comp_scalars = []
        if all(s in comp and comp[s].size for s in members):
            stack = np.vstack([comp[s] for s in members])
            centroid = stack.mean(axis=0)
            csum = centroid.sum()
            if csum > 0:
                centroid = centroid / csum
                for s in members:
                    d = jensenshannon(comp[s], centroid, base=2)
                    comp_scalars.append(0.0 if not np.isfinite(d) else float(d))
        len_scalars = [float(lenv[s][1]) for s in members if s in lenv and len(lenv[s]) >= 2]
        impr_scalars = [float(impr[s]) for s in members if s in impr]
        group_mad[lv] = {
            "mad_category_composition": _mad(comp_scalars) if len(comp_scalars) == len(members) else np.nan,
            "mad_length_profile": _mad(len_scalars) if len(len_scalars) == len(members) else np.nan,
            "mad_imprecision": _mad(impr_scalars) if len(impr_scalars) == len(members) else np.nan,
        }

    rows = []
    for lv, members in multi.items():
        for s in members:
            others = [o for o in members if o != s]
            comp_ag = np.nan
            if s in comp and comp[s].size and all(o in comp and comp[o].size for o in others):
                jsds = []
                for o in others:
                    d = jensenshannon(comp[s], comp[o], base=2)
                    jsds.append(0.0 if not np.isfinite(d) else d)
                comp_ag = 1.0 - float(np.mean(jsds))
            len_ag = (float(np.mean([_rel_agreement(lenv[s], lenv[o]) for o in others]))
                      if s in lenv and all(o in lenv for o in others) else np.nan)
            impr_ag = (float(np.mean([_rel_agreement([impr[s]], [impr[o]]) for o in others]))
                       if s in impr and all(o in impr for o in others) else np.nan)
            rows.append({"sampleID": s, "group": lv,
                         "category_composition": comp_ag,
                         "length_profile": len_ag, "imprecision": impr_ag,
                         **group_mad[lv]})
    return {"enabled": True, "reason": "", "samples": [r["sampleID"] for r in rows],
            "groups": groups, "axes": _A5_AXES, "mad_axes": _A5_MAD_COLS,
            "per_sample": pd.DataFrame(rows)}


_F7_OFFSET_TOL = 3.0  # bp; per-site offsets within this are treated as agreeing


def compute_fuzziness_concordance(site_offset_DF, factor_map):
    """Within-group concordance of per-reference-site splice precision (F7).

    For replicate samples, whether the same reference splice sites are placed at
    the same sub-bp offsets. Per factor group, restrict to reference sites detected
    in ≥2 replicates and score each sample's agreement (1 − mean |Δoffset| capped
    at ±{tol} bp) with its group-mates on those shared sites. High agreement means
    the imprecision pattern is reproducible; a replicate that disagrees has
    sample-specific boundary placement. Disabled with a reason when there is no
    factor or no ≥2-replicate group.
    """
    def disabled(reason):
        return {"enabled": False, "reason": reason, "samples": [], "groups": {},
                "per_sample": pd.DataFrame(columns=["sampleID", "group",
                                                    "site_precision_agreement"])}
    if not factor_map:
        return disabled("no design factor supplied")
    if site_offset_DF is None or site_offset_DF.empty:
        return disabled("no splice-site offset data")
    groups = {}
    for s, lv in factor_map.items():
        groups.setdefault(str(lv), []).append(str(s))
    multi = {lv: ss for lv, ss in groups.items() if len(ss) >= 2}
    if not multi:
        return disabled("no factor level has ≥2 replicates")

    df = site_offset_DF.copy()
    df["sampleID"] = df["sampleID"].astype(str)
    df["site"] = df["ref_key"].astype(str) + "|" + df["site_type"].astype(str)
    rows = []
    for lv, members in multi.items():
        piv = df[df["sampleID"].isin(members)].pivot_table(
            index="site", columns="sampleID", values="median_offset", aggfunc="median")
        piv = piv.reindex(columns=members)
        piv = piv[piv.notna().sum(axis=1) >= 2]   # sites in >=2 replicates
        for s in members:
            ags = []
            for o in members:
                if o == s:
                    continue
                pair = piv[[s, o]].dropna()
                if len(pair):
                    d = np.abs(pair[s].to_numpy(dtype=float) - pair[o].to_numpy(dtype=float))
                    ags.append(float(np.mean(np.clip(1.0 - d / _F7_OFFSET_TOL, 0.0, 1.0))))
            rows.append({"sampleID": s, "group": lv,
                         "site_precision_agreement": float(np.mean(ags)) if ags else np.nan})
    return {"enabled": True, "reason": "", "samples": [r["sampleID"] for r in rows],
            "groups": groups, "per_sample": pd.DataFrame(rows)}


def compute_fuzz_depth_metrics(fuzz_depth_DF):
    """Per-sample splice-site imprecision vs junction depth-bin (F4).

    Imprecision concentrated at low-depth junctions is consistent with limited
    read support; imprecision that persists at high depth points to a systematic
    placement offset. Empty-safe (the fixture has no per-junction coverage, so
    this returns no samples and the page/section are skipped).
    """
    empty = {"samples": [],
             "profile": pd.DataFrame(columns=["sampleID", "depth_bin", "perc_imprecise"])}
    if fuzz_depth_DF is None or fuzz_depth_DF.empty:
        return empty
    g = (fuzz_depth_DF.groupby(["sampleID", "depth_bin"], observed=True)
         .agg(n_sites=("n_sites", "sum"), n_imprecise=("n_imprecise", "sum")).reset_index())
    g["perc_imprecise"] = np.where(g["n_sites"] > 0, g["n_imprecise"] / g["n_sites"] * 100, np.nan)
    return {"samples": sorted(g["sampleID"].astype(str).unique().tolist()),
            "profile": g[["sampleID", "depth_bin", "perc_imprecise"]]}


def _robust_z(values):
    """Robust z-scores of a 1-D sequence: (x - median) / (1.4826 * MAD).

    1.4826 makes the MAD a consistent estimator of the SD under normality, so the
    z threshold reads on a familiar scale. Falls back to the standard-deviation
    z if MAD is 0 (all-but-one identical), and returns all-zeros if that is 0 too
    (every value identical => no outliers). NaNs propagate as 0 contribution.
    """
    v = np.asarray(values, dtype=float)
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med))
    if mad > 0:
        z = (v - med) / (1.4826 * mad)
    else:
        sd = np.nanstd(v)
        z = (v - med) / sd if sd > 0 else np.zeros_like(v)
    return np.nan_to_num(z, nan=0.0)


def compute_sample_scorecard(sample_metrics, cfg):
    """Cohort-relative sample-outlier scorecard.

    ``sample_metrics`` is a list of per-sample dicts (one per sample) carrying the
    raw metric values named in cfg['sample_scorecard']['metrics']; missing keys
    are simply not scored. Each metric is turned into a robust z-score across the
    COHORT (median/MAD over samples), oriented so a positive ``signed_z`` always
    means "worse" (large positive z for 'high'-is-worse metrics, sign-flipped for
    'low'-is-worse). A sample-metric cell is flagged warn/fail on |z| thresholds;
    a sample's overall flag trips when enough cells are flagged. This is
    deliberately relative and threshold-free at the metric level, so it is
    dataset-agnostic and stays quiet when all samples agree.

    Returns a dict:
      - 'enabled'   : bool (False if too few samples or no usable metrics)
      - 'reason'    : str when disabled
      - 'samples'   : sample order
      - 'metrics'   : metric names actually scored
      - 'z'         : {sampleID: {metric: signed_z}}
      - 'raw'       : {sampleID: {metric: raw_value}}
      - 'cell_flags': {sampleID: {metric: 'pass'|'warn'|'fail'}}
      - 'n_flagged' : {sampleID: int}   (warn+fail cells)
      - 'overall'   : {sampleID: 'pass'|'warn'|'fail'}
    """
    sc = (cfg or {}).get("sample_scorecard", {})
    samples = [m["sampleID"] for m in sample_metrics]
    metric_dirs = sc.get("metrics", {})
    min_samples = sc.get("min_samples", 4)

    result = {"enabled": False, "samples": samples, "metrics": [],
              "z": {}, "raw": {}, "cell_flags": {}, "n_flagged": {}, "overall": {}}

    if len(samples) < min_samples:
        result["reason"] = (f"cohort has {len(samples)} sample(s); the relative "
                            f"scorecard needs >= {min_samples} to be meaningful")
        return result

    # Keep only metrics present (non-NaN) for at least half the cohort.
    usable = []
    for metric in metric_dirs:
        vals = [m.get(metric) for m in sample_metrics]
        n_ok = sum(1 for v in vals if v is not None and not (isinstance(v, float) and np.isnan(v)))
        if n_ok >= max(3, len(samples) // 2):
            usable.append(metric)
    if not usable:
        result["reason"] = "no metric available for enough samples to score"
        return result

    z_warn, z_fail = sc.get("z_warn", 2.5), sc.get("z_fail", 3.5)
    n_warn = sc.get("min_metrics_warn", 2)
    n_fail = sc.get("min_metrics_fail", 2)

    z_by_metric = {}
    for metric in usable:
        raw = np.array([(m.get(metric) if m.get(metric) is not None else np.nan)
                        for m in sample_metrics], dtype=float)
        z = _robust_z(raw)
        if metric_dirs[metric] == "low":   # lower is worse -> flip so +z == worse
            z = -z
        z_by_metric[metric] = z

    for i, s in enumerate(samples):
        zrow, rawrow, cellrow = {}, {}, {}
        warn_hits = fail_hits = 0
        for metric in usable:
            zz = float(z_by_metric[metric][i])
            zrow[metric] = zz
            rawrow[metric] = sample_metrics[i].get(metric)
            # Only a "worse" deviation (positive signed z) is a flag; a sample far
            # on the GOOD side is never penalised.
            if zz >= z_fail:
                cellrow[metric] = "fail"; fail_hits += 1
            elif zz >= z_warn:
                cellrow[metric] = "warn"; warn_hits += 1
            else:
                cellrow[metric] = "pass"
        result["z"][s] = zrow
        result["raw"][s] = rawrow
        result["cell_flags"][s] = cellrow
        result["n_flagged"][s] = warn_hits + fail_hits
        if fail_hits >= n_fail:
            result["overall"][s] = "fail"
        elif (warn_hits + fail_hits) >= n_warn:
            result["overall"][s] = "warn"
        else:
            result["overall"][s] = "pass"

    result["enabled"] = True
    result["metrics"] = usable
    return result


def _ism_fragment_pct(ISM_DF):
    """{sampleID: % of ISM reads that are 3'/5'/internal fragments}.

    Splits the ISM bucket into its fragment sub-types (3'/5'/internal) versus the
    rest (mono-exon, intron retention). A high fragment fraction is a shape signal,
    NOT a quality verdict: it is consistent with truncation/degradation but equally
    with genuine alternative isoforms (alternative TSS/TTS, shorter isoforms). It is
    interpreted cohort-relatively (does a sample's ISM shape diverge from its
    peers?), never against an absolute cut-off. Empty dict when ISM_DF is
    absent/empty."""
    if ISM_DF is None or "sampleID" not in getattr(ISM_DF, "columns", []):
        return {}
    frag_cols = ["3prime_fragment", "5prime_fragment", "internal_fragment"]
    subcols = [c for c in ISM_DF.columns if c != "sampleID"
               and pd.api.types.is_numeric_dtype(ISM_DF[c])]
    if not subcols:
        return {}
    idf = ISM_DF.set_index("sampleID")[subcols].astype(float)
    tot = idf.sum(axis=1)
    present = [c for c in frag_cols if c in idf.columns]
    frag = idf[present].sum(axis=1) if present else tot * 0
    return (frag / tot * 100).replace([np.inf, -np.inf], np.nan).to_dict()


def _novel_noncanonical_pct(nov_can_DF):
    """{sampleID: % of all junctions that are novel AND non-canonical} — the
    junction class most enriched for alignment/calling artefacts (though not
    exclusively artefactual; non-canonical splicing also occurs biologically).
    Empty dict when the four known/novel x canonical columns aren't present."""
    if nov_can_DF is None or "sampleID" not in getattr(nov_can_DF, "columns", []):
        return {}
    four = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
    if not all(c in nov_can_DF.columns for c in four):
        return {}
    nj = nov_can_DF.set_index("sampleID")[four].astype(float)
    tot = nj.sum(axis=1)
    return (nj['novel_non_canonical'] / tot * 100
            ).replace([np.inf, -np.inf], np.nan).to_dict()


def assemble_scorecard_metrics(samples, length_DF=None, err_DF=None,
                               all_gene_percs_pivot_DF=None, nov_can_DF=None,
                               completeness_metrics=None, jxn_offset_metrics=None,
                               ISM_DF=None, yield_metrics=None, drift_metrics=None,
                               support_DF=None):
    """Gather the raw per-sample metric values the scorecard scores.

    Pulls each metric from the table that already computed it (no recomputation),
    returning a list of per-sample dicts keyed by sampleID + the scorecard metric
    names. Any source missing simply leaves that metric out (the scorecard drops
    metrics that aren't present for enough samples), so this is safe across runs
    with different SQANTI3 feature sets.
    """
    def _idx(df, col):
        if df is None or col not in getattr(df, "columns", []):
            return {}
        try:
            return df.set_index("sampleID")[col].astype(float).to_dict()
        except Exception:
            return {}

    med_len = _idx(length_DF, "median_length")
    rts = _idx(err_DF, "perc_reads_RTS")
    ip = _idx(err_DF, "perc_reads_intrapriming")
    ism = _idx(all_gene_percs_pivot_DF, "ISM")

    # % novel junctions from the known/novel × canonical junction counts.
    novj = {}
    if nov_can_DF is not None and "sampleID" in getattr(nov_can_DF, "columns", []):
        four = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
        if all(c in nov_can_DF.columns for c in four):
            nj = nov_can_DF.set_index("sampleID")[four].astype(float)
            tot = nj.sum(axis=1)
            novj = ((nj['novel_canonical'] + nj['novel_non_canonical']) / tot * 100
                    ).replace([np.inf, -np.inf], np.nan).to_dict()

    comp = {}
    if completeness_metrics and not completeness_metrics.get("per_sample", pd.DataFrame()).empty:
        cp = completeness_metrics["per_sample"].set_index("sampleID")
        comp = {"perc_5p_within_window": cp["perc_5p_within_window"].to_dict(),
                "perc_3p_within_window": cp["perc_3p_within_window"].to_dict()}

    impr = {}
    if jxn_offset_metrics and not jxn_offset_metrics.get("per_sample", pd.DataFrame()).empty:
        impr = jxn_offset_metrics["per_sample"].set_index("sampleID")["perc_imprecise"].to_dict()

    # Derived quality scalars (single source of truth: the module helpers).
    ism_frag = _ism_fragment_pct(ISM_DF)
    novnc = _novel_noncanonical_pct(nov_can_DF)

    # Between-sample comparison scalars (A4 depth-normalised yield; A3 drift).
    yld = _idx(yield_metrics.get("per_sample") if yield_metrics else None, "jxn_per_1k_reads")
    drift = _idx(drift_metrics.get("per_sample") if drift_metrics else None, "mean_jsd")

    # B11/B12/B13 orthogonal-support scalars (absent unless the assay was supplied).
    cage = _idx(support_DF, "perc_within_cage")
    polya = _idx(support_DF, "perc_within_polya")
    sr = _idx(support_DF, "perc_jxn_SR_supported")

    out = []
    for s in samples:
        out.append({
            "sampleID": s,
            "median_length": med_len.get(s),
            "perc_ISM": ism.get(s),
            "perc_novel_junctions": novj.get(s),
            "perc_5p_within_window": comp.get("perc_5p_within_window", {}).get(s),
            "perc_3p_within_window": comp.get("perc_3p_within_window", {}).get(s),
            "perc_reads_RTS": rts.get(s),
            "perc_reads_intrapriming": ip.get(s),
            "perc_sites_imprecise": impr.get(s),
            "perc_ISM_fragments": ism_frag.get(s),
            "perc_novel_noncanonical_jxn": novnc.get(s),
            "jxn_per_1k_reads": yld.get(s),
            "composition_drift": drift.get(s),
            "perc_within_cage": cage.get(s),
            "perc_within_polya": polya.get(s),
            "perc_jxn_SR_supported": sr.get(s),
        })
    return out


def plot_completeness_pages(pdf, completeness_metrics):
    """Append the 5'/3' completeness profile page (|distance|-to-gene-end ECDFs)."""
    m = completeness_metrics
    if not m or not m.get("samples"):
        return
    samples = m["samples"]
    window = m["window"]
    sample_palette = {s: sample_seq[i % len(sample_seq)] for i, s in enumerate(samples)}
    matplotlib.rcParams['pdf.fonttype'] = 42

    prof = m["profile"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, end, title in zip(axes, ["5prime", "3prime"],
                              ["5' completeness", "3' completeness"]):
        ed = prof[prof["end"] == end]
        for s in samples:
            d = ed[ed["sampleID"] == s].sort_values("k")
            if not d.empty:
                ax.plot(d["k"], d["perc_within"], marker="o", markersize=2,
                        color=sample_palette[s], label=s)
        ax.axvline(window, color="#999999", lw=0.8, ls=":")
        ax.set_title(title)
        ax.set_xlabel(f"|distance to annotated gene {'5' if end=='5prime' else '3'}' end| (bp)")
        ax.set_xlim(0, min(200, 10 * window))
    axes[0].set_ylabel("Cumulative % of reads")
    axes[1].legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.suptitle("Read-end completeness profile", y=1.0, fontsize=16)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_scorecard_page(pdf, scorecard):
    """Append the cohort-relative sample-outlier scorecard as a heatmap page.

    Cells are signed robust z-scores (red = worse than cohort, blue = better);
    warn/fail cells are outlined. A sample's overall flag is shown as a colored
    strip beside its row."""
    sc = scorecard
    if not sc or not sc.get("enabled"):
        return
    samples = sc["samples"]
    metrics = sc["metrics"]
    matplotlib.rcParams['pdf.fonttype'] = 42

    Z = np.array([[sc["z"][s][mt] for mt in metrics] for s in samples], dtype=float)
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(metrics) + 3),
                                    max(4, 0.6 * len(samples) + 2)))
    vlim = max(3.5, float(np.abs(Z).max()) if Z.size else 3.5)
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([mt.replace("perc_", "%").replace("_", " ") for mt in metrics],
                       rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=8)
    # annotate z + outline flagged cells
    flag_edge = {"warn": "#FF9800", "fail": "#F44336"}
    for i, s in enumerate(samples):
        for j, mt in enumerate(metrics):
            ax.text(j, i, f"{Z[i, j]:.1f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(Z[i, j]) > vlim * 0.6 else "black")
            fl = sc["cell_flags"][s][mt]
            if fl in flag_edge:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor=flag_edge[fl], lw=2.5))
    # overall-flag strip on the right
    overall_color = {"pass": "#4CAF50", "warn": "#FF9800", "fail": "#F44336"}
    for i, s in enumerate(samples):
        ax.add_patch(plt.Rectangle((len(metrics) - 0.4, i - 0.5), 0.5, 1,
                                   transform=ax.transData, clip_on=False,
                                   color=overall_color[sc["overall"][s]]))
    ax.set_title("Sample-outlier scorecard (cohort-relative robust z; red = worse than peers)",
                 fontsize=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.12)
    cbar.set_label("robust z vs cohort (+ = worse)")
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_ujc_overlap_page(pdf, overlap):
    """A6 — pairwise UJC overlap-index heatmap (asymmetric ``|A∩B|/|A|``, Viridis
    0–1, cells annotated). No-op when fewer than two samples carry UJCs."""
    if not overlap or overlap.get("matrix") is None:
        return
    mat = overlap["matrix"]
    samples = overlap["samples"]
    M = mat.to_numpy(dtype=float)
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(samples) + 3),
                                    max(5, 0.7 * len(samples) + 2)))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=8)
    for i in range(len(samples)):
        for j in range(len(samples)):
            # Viridis is dark in its lower half — use white text on low cells.
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if M[i, j] < 0.55 else "black")
    ax.set_xlabel("fraction also detected in this sample (B)")
    ax.set_ylabel("fraction of this sample's UJCs (A)")
    ax.set_title("Pairwise UJC overlap index (row A ∩ col B / row A)", fontsize=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.08)
    cbar.set_label("overlap fraction (1 = all of A's UJCs seen in B)")
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_support_pages(pdf, cage_metrics, polya_metrics, sr_metrics, scorecard=None):
    """B11/B12/B13 orthogonal-support pages — one per-sample cohort bar per view.

    CAGE 5' support (B11), polyA 3' support (B12), and short-read support (B13:
    TSS ratio + junction coverage). Each view no-ops when its evidence column was
    not supplied to SQANTI3 (empty metrics). Cohort-relative only: a low value
    means this sample's ends/junctions are less independently supported than its
    peers', not a fixed pass/fail — so no warn/fail threshold lines are drawn."""
    if cage_metrics and cage_metrics.get("samples"):
        _plot_metric_cohort_page(
            pdf, cage_metrics["per_sample"], "perc_within_cage", scorecard=scorecard,
            title="CAGE 5'-end support: reads with a TSS inside a CAGE peak",
            ylabel="% of reads within a CAGE peak", color=aes_palette["blue"])
    if polya_metrics and polya_metrics.get("samples"):
        _plot_metric_cohort_page(
            pdf, polya_metrics["per_sample"], "perc_within_polya", scorecard=scorecard,
            title="PolyA 3'-end support: reads with a TTS at a polyA site",
            ylabel="% of reads at a polyA site", color=aes_palette["magenta"])
    if sr_metrics and sr_metrics.get("samples"):
        _plot_metric_cohort_page(
            pdf, sr_metrics["per_sample"], "median_ratio_TSS",
            title="Short-read 5' support: median TSS ratio",
            ylabel="median short-read TSS ratio", color=aes_palette["green"])
        _plot_metric_cohort_page(
            pdf, sr_metrics["per_sample"], "perc_jxn_SR_supported", scorecard=scorecard,
            title="Short-read junction support: junctions with short-read coverage",
            ylabel="% of junctions with short-read coverage", color=aes_palette["green"])


def _plot_metric_cohort_page(pdf, metric_DF, metric, *, title, ylabel,
                             threshold=None, scorecard=None, color=None):
    """One per-sample bar page that puts a single QC rate in cohort context.

    Draws each sample's value for ``metric`` as a bar, overlays the cohort median
    (dashed grey) and, when configured, the warn/fail threshold lines from the
    config, and marks any sample the scorecard flagged *on this metric* with a
    colored dot + border. This is the per-metric companion to the scorecard
    heatmap: the heatmap says "which samples diverge and by how much"; this page
    shows the actual rate behind that verdict and where the peers sit.

    ``metric_DF`` needs a ``sampleID`` column and ``metric``; robust to a missing
    metric (no page) and to no scorecard (bars + median only, no flag markers).
    """
    if metric_DF is None or metric not in getattr(metric_DF, "columns", []):
        return
    d = metric_DF[["sampleID", metric]].dropna()
    if d.empty:
        return
    samples = d["sampleID"].astype(str).tolist()
    vals = d[metric].astype(float).tolist()
    color = color or aes_palette["blue"]
    matplotlib.rcParams['pdf.fonttype'] = 42

    fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(samples) + 3), 6))
    bars = ax.bar(samples, vals, color=color, edgecolor="none", zorder=2)

    med = float(np.median(vals))
    ax.axhline(med, color="#777777", ls="--", lw=1.2, zorder=1,
               label=f"cohort median ({med:.1f})")

    # Config warn/fail lines (absolute reference; the scorecard itself is relative).
    if threshold:
        for lvl, ls in (("warn", ":"), ("fail", "-.")):
            if lvl in threshold:
                ax.axhline(float(threshold[lvl]), color=_FLAG_FG.get(lvl, "#B71C1C"),
                           ls=ls, lw=1.1, zorder=1,
                           label=f"{lvl} threshold ({threshold[lvl]:.0f})")

    # Mark samples the scorecard flagged on THIS metric.
    flag_color = {"warn": "#FF9800", "fail": "#F44336"}
    if scorecard and scorecard.get("enabled") and metric in scorecard.get("metrics", []):
        for bar, s in zip(bars, samples):
            fl = scorecard["cell_flags"].get(s, {}).get(metric, "pass")
            if fl in flag_color:
                bar.set_edgecolor(flag_color[fl])
                bar.set_linewidth(2.5)
                ax.plot(bar.get_x() + bar.get_width() / 2,
                        bar.get_height(), marker="o", markersize=8,
                        color=flag_color[fl], zorder=3)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Sample ID")
    ax.set_ylabel(ylabel)
    ax.margins(y=0.15)
    if len(samples) > 6:
        plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_artefact_metric_pages(pdf, err_DF, qc_flags, scorecard=None):
    """Per-sample cohort-context pages for the read-artefact QC rates.

    One page each for RT-switching and intra-priming (the two artefact rates that
    feed the scorecard), each showing the per-sample percentage against the cohort
    median, the config warn/fail thresholds, and the scorecard flag markers.
    No-ops on missing columns."""
    if err_DF is None:
        return
    specs = [
        ("perc_reads_RTS", "RT-switching reads per sample",
         "% reads with RT-switching evidence", aes_palette["orange"]),
        ("perc_reads_intrapriming", "Intra-priming reads per sample",
         "% reads with intra-priming evidence", aes_palette["magenta"]),
    ]
    for metric, title, ylabel, color in specs:
        _plot_metric_cohort_page(
            pdf, err_DF, metric, title=title, ylabel=ylabel,
            threshold=(qc_flags or {}).get(metric), scorecard=scorecard, color=color)


def plot_quality_metric_pages(pdf, ISM_DF, nov_can_DF, qc_flags, scorecard=None):
    """Per-sample cohort-context pages for two derived read-quality scalars that
    feed the scorecard:

      * ISM fragment fraction  -- what share of a sample's ISM reads are 3'/5'/
        internal fragments rather than mono-exon / intron-retention ISMs. This is a
        SHAPE signal, not a quality verdict: fragment-shaped ISMs arise from
        truncation/degradation but equally from genuine alternative isoforms, so it
        is read cohort-relatively (does this sample's ISM shape diverge from its
        peers?) and carries no absolute threshold;
      * Novel non-canonical junction burden -- the junction class most enriched for
        alignment/calling artefacts (not exclusively artefactual), separated from
        novel-canonical.

    Each rendered as bars against the cohort median, any config warn/fail lines,
    and scorecard flag markers. No-ops when the source table is missing."""
    frag = _ism_fragment_pct(ISM_DF)
    if frag:
        d = pd.DataFrame({"sampleID": list(frag.keys()),
                          "perc_ISM_fragments": list(frag.values())})
        _plot_metric_cohort_page(
            pdf, d, "perc_ISM_fragments",
            title="ISM fragment fraction per sample (shape, not quality)",
            ylabel="% of ISM reads that are 3'/5'/internal fragments",
            threshold=(qc_flags or {}).get("perc_ISM_fragments"),
            scorecard=scorecard, color=subcat_color_palette.get("3prime_fragment", aes_palette["orange"]))

    novnc = _novel_noncanonical_pct(nov_can_DF)
    if novnc:
        d = pd.DataFrame({"sampleID": list(novnc.keys()),
                          "perc_novel_noncanonical_jxn": list(novnc.values())})
        _plot_metric_cohort_page(
            pdf, d, "perc_novel_noncanonical_jxn",
            title="Novel non-canonical junction burden, per sample",
            ylabel="% of junctions that are novel & non-canonical",
            threshold=(qc_flags or {}).get("perc_novel_noncanonical_jxn"),
            scorecard=scorecard, color=jxn_palette.get("novel_non_canonical", aes_palette["magenta"]))


def plot_jxn_yield_page(pdf, yield_metrics, scorecard=None):
    """Per-sample depth-normalised junction yield (A4), in cohort context.

    A sample low here after depth-normalisation has genuinely lower junction
    complexity, not just fewer reads. No-op on empty."""
    m = yield_metrics
    if not m or not m.get("samples"):
        return
    _plot_metric_cohort_page(
        pdf, m["per_sample"], "jxn_per_1k_reads",
        title="Junction yield per 1000 reads (depth-normalised)",
        ylabel="distinct junctions / 1k reads",
        scorecard=scorecard, color=aes_palette["blue"])


def plot_composition_drift_page(pdf, drift):
    """Pairwise composition-drift heatmap (A3): Jensen–Shannon distance between
    samples' structural-category composition vectors. Sequential scale (distance
    is >= 0). No-op on empty."""
    m = drift
    if not m or not m.get("samples") or m.get("matrix") is None or m["matrix"].empty:
        return
    mat = m["matrix"]
    samples = list(mat.index)
    Z = mat.to_numpy(dtype=float)
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(samples) + 3),
                                    max(5, 0.6 * len(samples) + 2)))
    im = ax.imshow(Z, cmap="viridis", vmin=0, vmax=max(float(Z.max()), 1e-6), aspect="auto")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=8)
    for i in range(len(samples)):
        for j in range(len(samples)):
            ax.text(j, i, f"{Z[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if Z[i, j] > Z.max() * 0.6 else "black")
    ax.set_title("Composition drift — pairwise Jensen–Shannon distance\n"
                 "(structural-category composition; larger = more different)", fontsize=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.08)
    cbar.set_label("Jensen–Shannon distance")
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_tandem_sites_page(pdf, tandem):
    """Tandem-site (NAGNAG) page (F6): near-zero offset histogram (summed over
    samples) with the ±3/±4/±6 tandem bars highlighted, plus each sample's tandem
    fraction at ±3 bp. Descriptive — a higher tandem fraction is a property of the
    sample, not a defect. No-op on empty."""
    m = tandem
    if not m or not m.get("samples") or m.get("by_offset") is None or m["by_offset"].empty:
        return
    agg = (m["by_offset"].groupby("offset", as_index=False)["count"].sum()
           .sort_values("offset"))
    offs = agg["offset"].tolist()
    counts = agg["count"].tolist()
    tand = {3, -3, 4, -4, 6, -6}
    colors = [aes_palette["magenta"] if o in tand else "#BBBBBB" for o in offs]
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(offs, counts, color=colors, width=0.8, zorder=2)
    ax.set_xlabel("signed splice-site offset (bp)")
    ax.set_ylabel("imprecise observations (all samples)")
    ax.set_title("Tandem splice sites (NAGNAG): excess mass at ±3/±4/±6 bp "
                 "(highlighted)", fontsize=12)
    ps = m["per_sample"]
    note = "   ".join(f"{r.sampleID}: {r.perc_tandem_3bp:.1f}%" for r in ps.itertuples())
    ax.text(0.5, -0.22, "tandem fraction at ±3 bp — " + note,
            transform=ax.transAxes, ha="center", fontsize=8, color="#555555")
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_replicate_concordance_page(pdf, concordance):
    """Multi-axis replicate-concordance page (A5): per-sample within-group
    agreement on composition / length / imprecision, grouped by factor level.
    No-op when disabled (no factor or no ≥2-replicate level)."""
    m = concordance
    if not m or not m.get("enabled") or m["per_sample"].empty:
        return
    ps = m["per_sample"].reset_index(drop=True)
    axes = m["axes"]
    samples = ps["sampleID"].tolist()
    x = np.arange(len(samples))
    w = 0.8 / len(axes)
    colors = [aes_palette["green"], aes_palette["orange"], aes_palette["blue"]]
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(max(8, 1.0 * len(samples) + 3), 5.5))
    for i, axis in enumerate(axes):
        ax.bar(x + i * w, ps[axis].fillna(0).to_numpy(dtype=float), width=w,
               label=axis.replace("_", " "), color=colors[i % len(colors)], zorder=2)
    ax.set_xticks(x + w * (len(axes) - 1) / 2)
    ax.set_xticklabels([f"{s}\n({g})" for s, g in zip(samples, ps["group"])], fontsize=8)
    ax.set_ylabel("within-group agreement (1 = matches replicates)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Replicate concordance (multi-axis): agreement with group-mates\n"
                 "(only factor levels with ≥2 replicates shown)", fontsize=12)
    ax.legend(fontsize=8, frameon=False, ncol=len(axes))
    # A8 companion — per-group robust within-group spread (MAD) in each axis's own
    # units. Complements the agreement bars: how tightly the whole group clusters
    # (non-degenerate even for exactly two replicates, unlike pairwise agreement).
    mad_cols = m.get("mad_axes") or []
    if mad_cols and all(c in ps.columns for c in mad_cols):
        lines = []
        for g, gdf in ps.groupby("group"):
            r = gdf.iloc[0]
            def _fmt(v):
                return "n/a" if not np.isfinite(v) else f"{v:.3g}"
            lines.append(f"{g}: composition {_fmt(r['mad_category_composition'])} JSD · "
                         f"length {_fmt(r['mad_length_profile'])} bp · "
                         f"imprecision {_fmt(r['mad_imprecision'])} pp")
        fig.text(0.5, -0.02,
                 "Within-group spread (robust MAD; 0 = replicates identical) — "
                 + "   |   ".join(lines),
                 ha="center", va="top", fontsize=7, color="#555555", wrap=True)
    # In a 2-replicate group each sample's only group-mate is the other, and all
    # distances are symmetric, so the pair's bars are identical by construction.
    if (ps.groupby("group")["sampleID"].transform("size") == 2).any():
        fig.text(0.5, -0.08,
                 "Note: groups of exactly 2 replicates show identical bars for both samples "
                 "(the score is the pair's mutual agreement); per-sample differences appear "
                 "only with ≥3 replicates in a group.",
                 ha="center", va="top", fontsize=7, style="italic", wrap=True)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_fuzziness_concordance_page(pdf, concordance):
    """Replicate concordance of splice-site precision (F7): per-sample within-group
    agreement on per-reference-site median offsets, colored by factor level. No-op
    when disabled (no factor or no ≥2-replicate level)."""
    m = concordance
    if not m or not m.get("enabled") or m["per_sample"].empty:
        return
    ps = m["per_sample"].reset_index(drop=True)
    samples = ps["sampleID"].tolist()
    grps = ps["group"].tolist()
    uniq = list(dict.fromkeys(grps))
    gcolor = {g: sample_seq[i % len(sample_seq)] for i, g in enumerate(uniq)}
    x = np.arange(len(samples))
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(max(8, 0.9 * len(samples) + 3), 5.5))
    ax.bar(x, ps["site_precision_agreement"].fillna(0).to_numpy(dtype=float),
           color=[gcolor[g] for g in grps], zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n({g})" for s, g in zip(samples, grps)], fontsize=8)
    ax.set_ylabel("site-precision agreement (1 = matches replicates)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Replicate concordance of splice-site precision\n"
                 "(per reference site; only ≥2-replicate levels shown)", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_fuzz_depth_page(pdf, metrics):
    """Splice-site imprecision vs junction read-depth bin (F4): one line per
    sample. No-op on empty (the fixture has no per-junction coverage → 0 pages)."""
    m = metrics
    if not m or not m.get("samples") or m["profile"].empty:
        return
    prof = m["profile"].copy()
    order = _fuzz_depth_labels(_FUZZ_DEPTH_BINS)
    prof["depth_bin"] = prof["depth_bin"].astype(str)
    xpos = {lab: i for i, lab in enumerate(order)}
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, (s, g) in enumerate(prof.groupby("sampleID")):
        g = g[g["depth_bin"].isin(order)].sort_values(
            "depth_bin", key=lambda col: col.map(xpos))
        ax.plot([xpos[b] for b in g["depth_bin"]], g["perc_imprecise"],
                marker="o", label=str(s), color=sample_seq[i % len(sample_seq)])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=9)
    ax.set_xlabel("junction read depth (total_coverage_unique)")
    ax.set_ylabel("% imprecise splice-site observations")
    ax.set_title("Splice-site imprecision vs junction depth", fontsize=12)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_jxn_offset_pages(pdf, jxn_offset_metrics):
    """Append the splice-site fuzziness pages (spectrum + precision profile +
    directionality + canonical split) to `pdf`."""
    m = jxn_offset_metrics
    if not m or not m.get("samples"):
        return
    samples = m["samples"]
    window = m["window"]
    sample_palette = {s: sample_seq[i % len(sample_seq)] for i, s in enumerate(samples)}
    matplotlib.rcParams['pdf.fonttype'] = 42

    # Title page for the section.
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.5, 0.5, "Splice-site fuzziness analysis", ha="center", va="center",
             fontsize=24, fontweight="bold")
    pdf.savefig(fig)
    plt.close(fig)

    # 1. Signed-offset spectrum (exact matches excluded), donor + acceptor panels.
    spec = m["spectrum"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, site in zip(axes, ["donor", "acceptor"]):
        ss = spec[(spec["site"] == site) & (spec["offset"] != 0)]
        for s in samples:
            d = ss[ss["sampleID"] == s].sort_values("offset")
            if not d.empty:
                ax.plot(d["offset"], d["count"], color=sample_palette[s], lw=1.4, label=s)
        ax.axvline(0, color="#999999", lw=0.8, ls=":")
        ax.set_title(site.capitalize())
        ax.set_xlabel(f"Offset from reference {site} (bp)")
        ax.set_xlim(-window, window)
    axes[0].set_ylabel("Site observations (exact matches excluded)")
    axes[1].legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.suptitle("Splice-site offset spectrum", y=1.0, fontsize=16)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # 2. Precision profile: cumulative % within +/-k bp.
    prof = m["profile"]
    fig, ax = plt.subplots(figsize=(11, 7))
    for s in samples:
        d = prof[prof["sampleID"] == s].sort_values("k")
        ax.plot(d["k"], d["perc_within"], marker="o", markersize=3,
                color=sample_palette[s], label=s)
    ax.set_xlabel("Window half-width k (bp)")
    ax.set_ylabel("% of site observations within +/- k of reference")
    ax.set_title("Splice-site precision profile")
    ax.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # 3. Canonical vs non-canonical imprecision (grouped bars).
    bc = m["by_class"]
    if bc is not None and not bc.empty:
        classes = [c for c in ["canonical", "non_canonical"] if c in set(bc["canonical"])]
        classes += [c for c in sorted(set(bc["canonical"])) if c not in classes]
        class_colors = {"canonical": aes_palette["green"],
                        "non_canonical": aes_palette["magenta"]}
        x = np.arange(len(samples))
        w = 0.8 / max(1, len(classes))
        fig, ax = plt.subplots(figsize=(11, 7))
        for i, can in enumerate(classes):
            vals = [bc[(bc["sampleID"] == s) & (bc["canonical"] == can)]["perc_imprecise"]
                    for s in samples]
            vals = [float(v.iloc[0]) if len(v) else 0.0 for v in vals]
            ax.bar(x + (i - (len(classes) - 1) / 2) * w, vals, w,
                   color=class_colors.get(can, "#999999"), label=can.replace("_", " "))
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=90)
        ax.set_ylabel("% site observations imprecise (offset != 0)")
        ax.set_title("Splice-site fuzziness by canonical class")
        ax.legend(title="Splice-site class")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def plot_ujc_metrics_pages(pdf, ujc_metrics, factor_col=None):
    """Append the UJC saturation, replicate-concordance and UpSet pages to `pdf`."""
    if not ujc_metrics:
        return
    samples = ujc_metrics["samples"]
    sample_palette = {s: sample_seq[i % len(sample_seq)] for i, s in enumerate(samples)}

    # 1. Saturation / rarefaction curves
    sat = ujc_metrics["saturation"]
    plt.figure(figsize=(14, 10))
    for s in samples:
        d = sat[sat["sampleID"] == s]
        plt.plot(d["depth"], d["unique_ujcs"], marker="o", markersize=3,
                 label=s, color=sample_palette[s])
    plt.title("UJC saturation (rarefaction)")
    plt.xlabel("Reads sampled")
    plt.ylabel("Expected unique UJCs")
    plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    matplotlib.rcParams['pdf.fonttype'] = 42
    pdf.savefig()
    plt.close()

    # 1b. Saturation per structural category (same total-depth x-axis)
    sat_c = ujc_metrics.get("saturation_by_category")
    if sat_c is not None and not sat_c.empty:
        cats = [c for c in cat_order if c in set(sat_c["structural_category"])]
        cats += [c for c in sorted(set(sat_c["structural_category"])) if c not in cats]
        ncols = min(3, len(cats))
        nrows = int(np.ceil(len(cats) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                                 squeeze=False, sharex=True)
        for i, cat in enumerate(cats):
            ax = axes[i // ncols][i % ncols]
            cd = sat_c[sat_c["structural_category"] == cat]
            for s in samples:
                d = cd[cd["sampleID"] == s]
                ax.plot(d["depth"], d["unique_ujcs"], marker="o", markersize=2,
                        label=s, color=sample_palette[s])
            ax.set_title(cat)
            ax.set_xlabel("Reads sampled")
            ax.set_ylabel("Expected unique UJCs")
        # blank any unused axes
        for j in range(len(cats), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        handles, labels = axes[0][0].get_legend_handles_labels()
        fig.legend(handles, labels, title="Sample", loc="upper right")
        fig.suptitle("UJC saturation by structural category", y=1.0, fontsize=16)
        plt.tight_layout()
        matplotlib.rcParams['pdf.fonttype'] = 42
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    # 2. Replicate concordance heatmap
    conc = ujc_metrics["concordance"]
    plt.figure(figsize=(10, 8))
    sns.heatmap(conc, annot=True, fmt=".2f", cmap="viridis", vmin=0, vmax=1,
                square=True, cbar_kws={"label": "Pearson r"})
    plt.title("Replicate concordance (per-UJC read counts)")
    plt.tight_layout()
    matplotlib.rcParams['pdf.fonttype'] = 42
    _vectorize_colorbars(plt.gcf())
    pdf.savefig()
    plt.close()

    # 3. UpSet plot of shared UJCs, stacked by structural category
    _plot_upset_pdf(pdf, ujc_metrics["upset"], samples)


# CLI arguments are defined once in src/reads_argparse.py; this module is driven
# via run_reads_plots(). (A previous standalone getOptions()/__main__ entry point
# duplicated --report and the thresholds and has been removed.)

def load_sqanti_file(file,col_Lst, dtype_Dct, tolerant=False):
    """Load a tab-separated SQANTI3 table, requesting only ``col_Lst``.

    With ``tolerant=True`` the header is read first and only the requested columns
    that actually exist are loaded — so optional columns (e.g. the CAGE/polyA/
    short-read support columns, absent in older or evidence-free SQANTI3 runs) do
    not raise. Downstream code guards on ``_column_is_populated`` before using any
    optional column, so a missing column simply makes its view no-op."""
    if tolerant:
        header = pd.read_csv(file, sep="\t", nrows=0).columns
        present = [c for c in col_Lst if c in header]
        dtype_Dct = {k: v for k, v in dtype_Dct.items() if k in present}
        return pd.read_csv(file, sep="\t", usecols=present, dtype=dtype_Dct, low_memory=True)
    return pd.read_csv(file, sep="\t", usecols=col_Lst, dtype=dtype_Dct, low_memory=True)


def _column_is_populated(df, col):
    """True when ``col`` exists in ``df`` and has at least one non-null value.

    The empty-safe gate for every optional support column: an absent or all-NaN
    column (no CAGE/polyA/short-read evidence supplied to SQANTI3) reads as
    unpopulated, so the dependent metric/page/section degrades to a no-op."""
    return (df is not None and col in getattr(df, "columns", [])
            and bool(df[col].notna().any()))


def merge_dfs(df1, df2, column1, column2, how='outer'):
    """
    Merge two pandas DataFrames on specified columns with different IDs.
    Returns:
    - Merged DataFrame.
    """
    # Rename the column in df2 to match df1 for the merge
    df2_renamed = df2.rename(columns={column2: column1})
    
    # Perform the merge
    merged_df = pd.merge(df1, df2_renamed, on=column1, how=how)
    
    return merged_df

def merge_three_dfs(df1, df2, df3, col1, col2, col3):
    """
    Merge three DataFrames using an outer join based on one specified column from each DataFrame,
    ensuring these columns contain the same unique set of values.
    
    Parameters:
    - df1, df2, df3: the DataFrames to merge.
    - col1: the column name in df1 to merge on.
    - col2: the column name in df2 to merge on.
    - col3: the column name in df3 to merge on.
    
    Returns:
    - Merged DataFrame with all rows and columns from df1, df2, and df3.
    """
    
    # Check if the merge columns are unique in each DataFrame
    if not df1[col1].is_unique or not df2[col2].is_unique or not df3[col3].is_unique:
        raise ValueError("Merge columns must be unique in each DataFrame.")
    
    # Merge df1 and df2 on their specified columns, then rename the column in df2 to match df1 for simplicity
    merged_df1_df2 = pd.merge(df1, df2.rename(columns={col2: col1}), how='outer', on=col1)
    
    # Merge the result with df3, renaming the column in df3 to match df1
    final_merged_df = pd.merge(merged_df1_df2, df3.rename(columns={col3: col1}), how='outer', on=col1)
    
    return final_merged_df

def fix_sample_order(df, design):

    if design is not None and 'sampleID' in df.columns:
        sample_order = design['sampleID'].unique().tolist()
        df['sampleID'] = pd.Categorical(df['sampleID'], categories=sample_order, ordered=True)
    return df

def categorize_by_readcount(read_count):
    if read_count == 1:
        return '1 read'
    elif read_count <= 10:
        return '2-10 reads'
    elif read_count <= 50:
        return '11-50 reads'
    elif read_count <= 100:
        return '50-100 reads'
    else:
        return '100+ reads'
    
def cv(x):
    if np.mean(x) == 0:
        return np.nan
    else:
        return np.std(x) / np.mean(x)
    
def calc_jxn_cv(jxnDF, classDF, refDF, dropFlag=True):
    """
    Calculate cv of reference donors and acceptors from sqanti jxn file
    
    Parameters:
    - jxnDF, dataframe of sqanti junction file
    - classDF: dataframe of sqanti classification file
    - refDF: dataframe with column 'gene_id' that contains all reference genes
    - dropFlag: flag to drop jxns in novel genes or jxns associated with multiple genes
    
    Returns:
    - Dataframe with cv of each reference donor and acceptor 
    """
    ##Drop all junctions where diff to nearest ref junction cannot be calculated
    filter_jxnDF=jxnDF.dropna(subset=['diff_to_Ref_start_site','diff_to_Ref_end_site'])
    
    filter_jxnDF=pd.merge(filter_jxnDF, classDF[['isoform', 'associated_gene']], on='isoform', how='left')
    
    
    filter_jxnDF['ref_junction_start']=(filter_jxnDF['genomic_start_coord']+filter_jxnDF['diff_to_Ref_start_site']).astype(int)
    filter_jxnDF['ref_junction_end']=(filter_jxnDF['genomic_end_coord']+filter_jxnDF['diff_to_Ref_end_site']).astype(int)
    filter_jxnDF['abs_diff_to_start']=abs(filter_jxnDF['diff_to_Ref_start_site'])
    filter_jxnDF['abs_diff_to_end']=abs(filter_jxnDF['diff_to_Ref_end_site'])
    

    
    cv_startDF =filter_jxnDF.groupby(['chrom', 'strand', 'ref_junction_start']).agg({'abs_diff_to_start': ['mean', 'std', cv, 'size'], 'associated_gene': [lambda x: '|'.join(set(x)), lambda x: len(set(x)) ]}).reset_index()
    cv_startDF.columns = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff','cv', 'count','associated_gene','gene_count']
    
    cv_startDF['flag_multi_gene'] = cv_startDF['gene_count'].apply(lambda x: 1 if x > 1 else 0)
    cv_startDF['flag_annotated_gene'] = np.where(cv_startDF['associated_gene'].isin(refDF['gene_id'].unique()), 1, 0)
    
    cv_startDF['flag_donor'] = cv_startDF['strand'].apply(lambda x: 1 if x == '+' else 0)
    cv_startDF['flag_acceptor'] = cv_startDF['strand'].apply(lambda x: 1 if x == '-' else 0)
    #cv_startDF['flag_single'] = cv_startDF['count'].apply(lambda x: 1 if x == 1 else 0)
    cv_startDF['flag_mean_0'] = cv_startDF['mean_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    #cv_startDF['flag_std_0'] = cv_startDF['std_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    cv_startDF['flag_std_0'] = cv_startDF['std_abs_diff'].apply(lambda x: 0 if pd.isna(x) else (1 if x == 0 else 0))
    #cv_startDF['flag_cv_0'] = cv_startDF['cv'].apply(lambda x: 1 if x == 0 else 0)
    #cv_startDF['flag_cv_0'] = cv_startDF['cv'].apply(lambda x: 0 if pd.isna(x) else (1 if x == 0 else 0))
    
    cv_endDF =filter_jxnDF.groupby(['chrom', 'strand', 'ref_junction_end']).agg({'abs_diff_to_end': ['mean', 'std', cv, 'size'], 'associated_gene': [lambda x: '|'.join(set(x)), lambda x: len(set(x)) ]}).reset_index()
    cv_endDF.columns = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff','cv', 'count','associated_gene','gene_count']
    
   
   
    cv_endDF['flag_multi_gene'] = cv_endDF['gene_count'].apply(lambda x: 1 if x > 1 else 0)
    cv_endDF['flag_annotated_gene'] =  np.where(cv_endDF['associated_gene'].isin(refDF['gene_id'].unique()), 1, 0)
    
    
    cv_endDF['flag_donor'] = cv_endDF['strand'].apply(lambda x: 1 if x == '-' else 0)
    cv_endDF['flag_acceptor'] = cv_endDF['strand'].apply(lambda x: 1 if x == '+' else 0)
    #cv_endDF['flag_single'] = cv_endDF['count'].apply(lambda x: 1 if x == 1 else 0)
    cv_endDF['flag_mean_0'] = cv_endDF['mean_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    #cv_endDF['flag_std_0'] = cv_endDF['std_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    cv_endDF['flag_std_0'] = cv_endDF['std_abs_diff'].apply(lambda x: 0 if pd.isna(x) else (1 if x == 0 else 0))
    #cv_endDF['flag_cv_0'] = cv_endDF['cv'].apply(lambda x: 1 if x == 0 else 0)
    #cv_endDF['flag_cv_0'] = cv_endDF['cv'].apply(lambda x: 0 if pd.isna(x) else (1 if x == 0 else 0))
    
    
    cvDF=pd.concat([cv_startDF,cv_endDF])
    cvDF['flag_ref_match'] = cvDF['flag_mean_0'].apply(lambda x: 1 if x == 1 else 0)
    #cvDF['flag_cv_0'] = cvDF.apply(lambda row: 1 if row['cv'] == 0 and row['flag_mean_0'] != 1 else 0, axis=1)
    #cvDF['flag_cv_0'] = cvDF.apply(lambda row: 1 if row['cv'] == 0 and row['flag_mean_0'] != 1 else 0, axis=1)
    #cvDF['flag_cv_0'] = cvDF.apply(lambda row: 1 if not np.isnan(row['cv']) and row['cv'] == 0 and row['flag_mean_0'] != 1 else 0, axis=1)
    #cvDF['flag_cv_gt_0'] = cvDF['cv'].apply(lambda x: 1 if x > 0 else 0)
    #cvDF['flag_cv_gt_0'] = cvDF['cv'].apply(lambda x: 1 if not np.isnan(x) and x > 0 else 0)
    cvDF['flag_cv_0'] = cvDF.apply(lambda row: 1 if pd.notna(row['cv']) and row['cv'] == 0 and row['flag_mean_0'] != 1 else 0, axis=1)
    cvDF['flag_cv_gt_0'] = cvDF['cv'].apply(lambda x: 1 if pd.notna(x) and x > 0 else 0)
    
    if dropFlag:
         
        #count=((cvDF['flag_annotated_gene'] != 1) | (cvDF['flag_multi_gene'] == 1)).sum()
        #print(str(count)+"donors and acceptors assigned to multiple and/or novel genes")
        
        drop=(cvDF['flag_annotated_gene'] != 1) | (cvDF['flag_multi_gene'] == 1)
        cvDF_filtered=cvDF[~drop]
        
        return cvDF_filtered
    else :
        return cvDF
    
def flag_ref_monoexon(inRef):
    
    """
    Creates a dataframe with all reference genes and identifies genes with monoexon transcri[td]
    
    Parameters:
    - inRef, path to reference gtf
    
    Returns:
    - Dataframe with all reference gene ids and flag if it has a monoexon transcript or not
    """
    # Get exon counts from reference GTF (comment='#' skips GTF header lines)
    refgtf = pd.read_csv(inRef, names=['chr','source','feature','start','end','score',
                                       'strand','frame','attribute'], sep="\t",
                         comment='#', low_memory=False)
    refexon = refgtf[refgtf["feature"] == "exon"].copy()

    # Vectorized extraction of gene_id / transcript_id from the GTF attribute
    # column (format: key "value"; ...). Replaces a fragile per-row parser and
    # yields the same unquoted IDs it produced.
    refexon["gene_id"] = refexon["attribute"].str.extract(r'gene_id "([^"]+)"')
    refexon["transcript_id"] = refexon["attribute"].str.extract(r'transcript_id "([^"]+)"')
    missing = refexon["transcript_id"].isna()
    if missing.any():
        reads_logger.warning(f"transcript_id not found on {int(missing.sum())} exon line(s) of {inRef}")

    # Get min exons per transcript for each gene to flag genes with at least one monoexon transcript
    refxcrpt = refexon.groupby(["gene_id", "transcript_id"])["feature"].count().reset_index().rename(columns={"feature": "num_exon"})
    refgene = refxcrpt.groupby("gene_id")["num_exon"].min().reset_index().rename(columns={"num_exon": "min_exon"})
    refgene["flag_ref_monoexon"] = np.where(
        refgene["min_exon"] == 1,
        1,
        0
    )
    return refgene

def _run_parallel(func, items, jobs):
    """Run func(item) over items, in parallel with a thread pool when jobs>1.

    Threads (not processes) because the per-sample work is dominated by file I/O
    (pandas read_csv releases the GIL during parsing); each sample writes only
    its own dict keys, which is thread-safe in CPython.
    """
    jobs = max(1, int(jobs or 1))
    items = list(items)
    if jobs <= 1 or len(items) <= 1:
        for it in items:
            func(it)
        return
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=min(jobs, len(items))) as ex:
        list(ex.map(func, items))


# --- Per-sample summary builders --------------------------------------------
# Each takes one sample's loaded/merged classification (and junction) frame and
# returns one summary table. They are pure (no shared state) so proc_samples'
# inner loop reads as an orchestration and each metric can be unit-tested.

_JXN_CLASSES = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']


def _summarize_junctions(jxn_DF):
    """Per-isoform counts of the four junction classes (known/novel × canonical)."""
    count_jxns_DF = jxn_DF.pivot_table(index='isoform',
                                       columns=['junction_category', 'canonical'],
                                       aggfunc='size',
                                       fill_value=0).reset_index()
    # Flatten the MultiIndex columns by joining with an underscore, except 'isoform'.
    count_jxns_DF.columns = ['isoform'] + ['_'.join([str(c) for c in col]).strip()
                                           for col in count_jxns_DF.columns[1:]]
    for col in _JXN_CLASSES:
        if col not in count_jxns_DF.columns:
            count_jxns_DF[col] = 0
    return count_jxns_DF[['isoform'] + _JXN_CLASSES]


def _summarize_gene_counts(class_DF, sampleID, exp_factor, exp_factor_val):
    """Per-gene read counts by structural category + unique-UJC count."""
    gene_category_count_DF = class_DF.pivot_table(index='associated_gene',
                                                  columns='structural_category',
                                                  aggfunc='count',
                                                  fill_value=0)['isoform'].reset_index()
    gene_category_count_DF['total_read_count'] = gene_category_count_DF.iloc[:, 1:].sum(axis=1)
    gene_UJC_count_DF = class_DF.groupby('associated_gene')['jxnHash'].nunique().reset_index(name='unique_jxnHash_counts')
    gene_count_DF = pd.merge(gene_category_count_DF, gene_UJC_count_DF, how='outer', on='associated_gene')
    gene_count_DF['sampleID'] = sampleID
    gene_count_DF[exp_factor] = exp_factor_val
    return gene_count_DF


def _summarize_subcategory(class_DF, categories, wanted, sampleID, exp_factor, exp_factor_val):
    """Per-subcategory read counts for the structural categories in ``wanted``.

    Returns an empty frame when none of ``wanted`` occurs in this sample (matches
    the historical FSM/ISM/NIC-NNC guards)."""
    df = pd.DataFrame()
    if any(w in categories for w in wanted):
        df = class_DF[class_DF['structural_category'].isin(wanted)].copy()
        df['sampleID'] = sampleID
        df = df.pivot_table(index='sampleID', columns='subcategory',
                            aggfunc='count', fill_value=0)['isoform'].reset_index()
        df[exp_factor] = exp_factor_val
    return df


def _summarize_ujc(class_DF, sampleID, exp_factor, exp_factor_val):
    """Per-UJC (jxnHash × gene × structural category) read counts + MEI flag."""
    ujc_group_cols = ['jxnHash', 'associated_gene', 'structural_category']
    ujc_count_DF = class_DF.groupby(ujc_group_cols).agg({
        'isoform': 'nunique'  # Count unique isoforms for this group
    }).reset_index()

    for col in _JXN_CLASSES:
        if col in class_DF.columns:
            first_vals = class_DF.groupby(ujc_group_cols)[col].first().reset_index()
            ujc_count_DF = pd.merge(ujc_count_DF, first_vals, on=ujc_group_cols, how='left')
        else:
            ujc_count_DF[col] = 0
    for col in ujc_count_DF.columns:
        if ujc_count_DF[col].dtype == 'int64':
            ujc_count_DF[col] = ujc_count_DF[col].fillna(0)
        elif ujc_count_DF[col].dtype == 'object' or ujc_count_DF[col].dtype.name == 'string':
            ujc_count_DF[col] = ujc_count_DF[col].fillna('0')
    ujc_count_DF.rename(columns={'isoform': 'read_count'}, inplace=True)
    ujc_count_DF['flag_MEI'] = ujc_count_DF.groupby('associated_gene')['read_count'] \
        .transform(lambda s: (s == s.max()).astype(int))
    ujc_count_DF['sampleID'] = sampleID
    ujc_count_DF[exp_factor] = exp_factor_val

    # Reorder columns to original layout
    desired_cols = [
        'jxnHash', 'read_count', 'associated_gene',
        'known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical',
        'structural_category', 'flag_MEI', 'sampleID', 'flag_annotated_gene',
    ]
    existing_desired = [c for c in desired_cols if c in ujc_count_DF.columns]
    remaining_cols = [c for c in ujc_count_DF.columns if c not in existing_desired]
    return ujc_count_DF[existing_desired + remaining_cols]


def _summarize_length(class_DF, sampleID, exp_factor, exp_factor_val):
    """Read-length summary stats (counts/percentages over 1/2/3 kb + quantiles)."""
    total_reads = len(class_DF)
    reads_gt_1kb = (class_DF['length'] > 1000).sum()
    reads_gt_2kb = (class_DF['length'] > 2000).sum()
    reads_gt_3kb = (class_DF['length'] > 3000).sum()
    length_summary_DF = pd.DataFrame({
        'total_reads': [total_reads],
        'reads_gt_1kb': [reads_gt_1kb],
        'reads_gt_2kb': [reads_gt_2kb],
        'reads_gt_3kb': [reads_gt_3kb],
        'perc_reads_gt_1kb': [(reads_gt_1kb / total_reads) * 100],
        'perc_reads_gt_2kb': [(reads_gt_2kb / total_reads) * 100],
        'perc_reads_gt_3kb': [(reads_gt_3kb / total_reads) * 100],
        'average_length': [class_DF['length'].mean()],
        'median_length': [class_DF['length'].median()],
        'min_length': [class_DF['length'].min()],
        'max_length': [class_DF['length'].max()],
        'q25_length': [class_DF['length'].quantile(0.25)],
        'q75_length': [class_DF['length'].quantile(0.75)],
        'sampleID': [sampleID],
    })
    length_summary_DF[exp_factor] = exp_factor_val
    return length_summary_DF


def _summarize_errors(class_DF, ip_cutoff, sampleID, exp_factor, exp_factor_val):
    """Artefact-read counts/percentages: RT-switching, intra-priming, non-canonical."""
    total_reads = len(class_DF)
    num_reads_RTS = (class_DF['RTS_stage'] == True).sum()
    num_reads_intrapriming = (class_DF['perc_A_downstream_TTS'] > ip_cutoff).sum()
    num_reads_non_can = (class_DF['all_canonical'] == 'non_canonical').sum()
    err_DF = pd.DataFrame({
        'num_reads_RTS': [num_reads_RTS],
        'perc_reads_RTS': [(num_reads_RTS / total_reads) * 100],
        'num_reads_intrapriming': [num_reads_intrapriming],
        'perc_reads_intrapriming': [(num_reads_intrapriming / total_reads) * 100],
        'num_reads_non-canonical': [num_reads_non_can],
        'perc_reads_non-canonical': [(num_reads_non_can / total_reads) * 100],
    })
    err_DF['sampleID'] = sampleID
    err_DF[exp_factor] = exp_factor_val
    return err_DF


def _summarize_jxn_offsets(jxn_DF, sampleID, exp_factor, exp_factor_val, window=15):
    """Compact per-sample histogram of signed splice-site offsets to the reference.

    Each junction contributes two *site observations* — a donor and an acceptor —
    whose signed distance to the nearest annotated site is ``diff_to_Ref_*``. The
    assignment is strand-aware: on ``+`` the junction start is the donor and the
    end is the acceptor; on ``-`` this is reversed. Offsets are collapsed to a
    per-(site, offset, canonical) count within +/- ``window`` bp (the near-reference
    "fuzziness" zone); observations farther than ``window`` (or with no computable
    diff) are genuinely novel sites and are excluded here — they are already
    covered by the known/novel junction plots.

    Returns a long DataFrame [sampleID, exp_factor, site, offset, canonical, count]
    from which the offset spectrum, precision profile, directionality and the
    per-sample imprecision scalar are all derived downstream. Collapsing to counts
    keeps memory bounded (<=~4*(2*window+1) rows/sample) regardless of read depth.
    """
    cols = ["sampleID", exp_factor, "site", "offset", "canonical", "count"]
    if jxn_DF is None or jxn_DF.empty:
        return pd.DataFrame(columns=cols)

    df = jxn_DF.copy()
    for c in ("diff_to_Ref_start_site", "diff_to_Ref_end_site"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    can = df["canonical"] if "canonical" in df.columns else pd.Series("NA", index=df.index)

    plus, minus = df["strand"] == "+", df["strand"] == "-"
    parts = [
        pd.DataFrame({"site": "donor",    "offset": df["diff_to_Ref_start_site"].where(plus),  "canonical": can}),
        pd.DataFrame({"site": "acceptor", "offset": df["diff_to_Ref_end_site"].where(plus),    "canonical": can}),
        pd.DataFrame({"site": "acceptor", "offset": df["diff_to_Ref_start_site"].where(minus), "canonical": can}),
        pd.DataFrame({"site": "donor",    "offset": df["diff_to_Ref_end_site"].where(minus),   "canonical": can}),
    ]
    long = pd.concat(parts, ignore_index=True).dropna(subset=["offset"])
    long = long[long["offset"].abs() <= window]
    if long.empty:
        return pd.DataFrame(columns=cols)
    long["offset"] = long["offset"].astype(int)
    long["canonical"] = long["canonical"].fillna("NA").astype(str)

    out = (long.groupby(["site", "offset", "canonical"], observed=True)
                .size().reset_index(name="count"))
    out["sampleID"] = sampleID
    out[exp_factor] = exp_factor_val
    return out[cols]


def _summarize_site_offsets(jxn_DF, sampleID, exp_factor, exp_factor_val, window=15):
    """Per reference splice site, this sample's median signed offset (F7).

    Keeps the reference-site identity that ``_summarize_jxn_offsets`` collapses
    away, so replicates can be compared site-by-site. The reference coordinate is
    reconstructed exactly as ``calc_jxn_cv`` does (``genomic_*_coord +
    diff_to_Ref_*``); the donor/acceptor label is strand-aware (matching
    ``_summarize_jxn_offsets``). Only near-reference observations (|offset| <=
    window) are kept. Typed-empty safe.
    Returns [sampleID, exp_factor, ref_key, site_type, median_offset, n].
    """
    cols = ["sampleID", exp_factor, "ref_key", "site_type", "median_offset", "n"]
    need = ["chrom", "strand", "genomic_start_coord", "genomic_end_coord",
            "diff_to_Ref_start_site", "diff_to_Ref_end_site"]
    if jxn_DF is None or jxn_DF.empty or not all(c in jxn_DF.columns for c in need):
        return pd.DataFrame(columns=cols)
    df = jxn_DF.copy()
    for c in ("genomic_start_coord", "genomic_end_coord",
              "diff_to_Ref_start_site", "diff_to_Ref_end_site"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    plus, minus = df["strand"] == "+", df["strand"] == "-"

    def _part(site_type, mask, coord_col, diff_col):
        d = df[mask]
        return pd.DataFrame({"chrom": d["chrom"], "strand": d["strand"],
                             "ref_coord": d[coord_col] + d[diff_col],
                             "offset": d[diff_col], "site_type": site_type})

    parts = [
        _part("donor", plus, "genomic_start_coord", "diff_to_Ref_start_site"),
        _part("acceptor", plus, "genomic_end_coord", "diff_to_Ref_end_site"),
        _part("acceptor", minus, "genomic_start_coord", "diff_to_Ref_start_site"),
        _part("donor", minus, "genomic_end_coord", "diff_to_Ref_end_site"),
    ]
    long = pd.concat(parts, ignore_index=True).dropna(subset=["ref_coord", "offset"])
    long = long[long["offset"].abs() <= window]
    if long.empty:
        return pd.DataFrame(columns=cols)
    long["ref_key"] = (long["chrom"].astype(str) + ":" + long["strand"].astype(str)
                       + ":" + long["ref_coord"].astype("int64").astype(str))
    out = (long.groupby(["ref_key", "site_type"], observed=True)["offset"]
                .agg(median_offset="median", n="size").reset_index())
    out["sampleID"] = sampleID
    out[exp_factor] = exp_factor_val
    return out[cols]


_FUZZ_DEPTH_BINS = (1, 5, 10, 25, 50, 100)


def _fuzz_depth_labels(bins):
    labels = [f"{bins[i]}-{bins[i + 1] - 1}" for i in range(len(bins) - 1)]
    return labels + [f"{bins[-1]}+"]


def _summarize_fuzz_by_depth(jxn_DF, sampleID, exp_factor, exp_factor_val,
                             window=15, depth_bins=_FUZZ_DEPTH_BINS):
    """Splice-site imprecision as a function of junction read depth (F4).

    Each donor/acceptor observation (strand-aware, within ±window) is bucketed by
    its junction's ``total_coverage_unique`` into depth bins, counting imprecise
    (offset≠0) vs total. Requires per-junction coverage; when that column is
    absent or all-NaN (the chr22 fixture reality), returns a typed-empty frame so
    F4 no-ops. Returns [sampleID, exp_factor, depth_bin, n_sites, n_imprecise].
    """
    cols = ["sampleID", exp_factor, "depth_bin", "n_sites", "n_imprecise"]
    if (jxn_DF is None or jxn_DF.empty
            or "total_coverage_unique" not in jxn_DF.columns):
        return pd.DataFrame(columns=cols)
    cov = pd.to_numeric(jxn_DF["total_coverage_unique"], errors="coerce")
    if cov.notna().sum() == 0:   # all-NaN -> nothing to bin (fixture case)
        return pd.DataFrame(columns=cols)
    df = jxn_DF.copy()
    df["_cov"] = cov
    for c in ("diff_to_Ref_start_site", "diff_to_Ref_end_site"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    plus, minus = df["strand"] == "+", df["strand"] == "-"
    parts = [
        pd.DataFrame({"offset": df["diff_to_Ref_start_site"].where(plus),  "_cov": df["_cov"]}),
        pd.DataFrame({"offset": df["diff_to_Ref_end_site"].where(plus),    "_cov": df["_cov"]}),
        pd.DataFrame({"offset": df["diff_to_Ref_start_site"].where(minus), "_cov": df["_cov"]}),
        pd.DataFrame({"offset": df["diff_to_Ref_end_site"].where(minus),   "_cov": df["_cov"]}),
    ]
    long = pd.concat(parts, ignore_index=True).dropna(subset=["offset", "_cov"])
    long = long[long["offset"].abs() <= window]
    if long.empty:
        return pd.DataFrame(columns=cols)
    edges = list(depth_bins) + [np.inf]
    labels = _fuzz_depth_labels(depth_bins)
    long["depth_bin"] = pd.cut(long["_cov"], bins=edges, right=False, labels=labels)
    long = long.dropna(subset=["depth_bin"])
    long["imprecise"] = (long["offset"] != 0).astype(int)
    out = (long.groupby("depth_bin", observed=True)
                .agg(n_sites=("offset", "size"), n_imprecise=("imprecise", "sum"))
                .reset_index())
    out["sampleID"] = sampleID
    out[exp_factor] = exp_factor_val
    return out[cols]


def _summarize_support(class_DF, jxn_DF, sampleID, exp_factor, exp_factor_val):
    """B11/B12/B13 — one row of orthogonal-support scalars for a sample.

    From the classification table's optional evidence columns (present only when
    CAGE / polyA / short-read data were supplied to SQANTI3):
      perc_within_cage   B11 % of reads whose 5' end falls inside a CAGE peak
      median_dist_cage   B11 median |distance| (bp) from the 5' end to the nearest CAGE peak
      perc_within_polya  B12 % of reads whose 3' end sits at an annotated polyA site
      perc_polya_motif   B12 % of reads with a polyA motif found near the 3' end
      median_dist_polya  B12 median |distance| (bp) to the nearest polyA site
      median_ratio_TSS   B13 median short-read TSS ratio (5'-end coverage support)
    and from the junction table:
      perc_jxn_SR_supported  B13 % of junctions with short-read coverage (>0 unique reads)
    Every scalar is NaN when its source column is absent or all-NaN, so a sample
    with no evidence contributes an all-NaN row and the downstream views no-op."""
    def _truthy_frac(df, col):
        if not _column_is_populated(df, col):
            return np.nan
        s = df[col].dropna().astype(str).str.upper()
        return float((s == "TRUE").mean() * 100) if len(s) else np.nan

    def _median_abs(df, col):
        if not _column_is_populated(df, col):
            return np.nan
        v = pd.to_numeric(df[col], errors="coerce").dropna().abs()
        return float(v.median()) if len(v) else np.nan

    def _median(df, col):
        if not _column_is_populated(df, col):
            return np.nan
        v = pd.to_numeric(df[col], errors="coerce").dropna()
        return float(v.median()) if len(v) else np.nan

    perc_jxn_sr = np.nan
    if _column_is_populated(jxn_DF, "total_coverage_unique"):
        cov = pd.to_numeric(jxn_DF["total_coverage_unique"], errors="coerce").dropna()
        perc_jxn_sr = float((cov > 0).mean() * 100) if len(cov) else np.nan

    return pd.DataFrame([{
        "sampleID": sampleID, exp_factor: exp_factor_val,
        "perc_within_cage": _truthy_frac(class_DF, "within_CAGE_peak"),
        "median_dist_cage": _median_abs(class_DF, "dist_to_CAGE_peak"),
        "perc_within_polya": _truthy_frac(class_DF, "within_polyA_site"),
        "perc_polya_motif": _truthy_frac(class_DF, "polyA_motif_found"),
        "median_dist_polya": _median_abs(class_DF, "dist_to_polyA_site"),
        "median_ratio_TSS": _median(class_DF, "ratio_TSS"),
        "perc_jxn_SR_supported": perc_jxn_sr,
    }])


_SUPPORT_COLS = ["perc_within_cage", "median_dist_cage", "perc_within_polya",
                 "perc_polya_motif", "median_dist_polya", "median_ratio_TSS",
                 "perc_jxn_SR_supported"]


def _support_metrics(support_DF, cols):
    """Shared empty-safe reducer for the B views: returns ``{samples, per_sample}``
    where ``per_sample`` is ``support_DF`` restricted to sampleID + ``cols`` with at
    least one populated column, or ``{'samples': []}`` when none are populated."""
    if support_DF is None or support_DF.empty:
        return {"samples": [], "per_sample": pd.DataFrame(columns=["sampleID"] + cols)}
    have = [c for c in cols if _column_is_populated(support_DF, c)]
    if not have:
        return {"samples": [], "per_sample": pd.DataFrame(columns=["sampleID"] + cols)}
    keep = ["sampleID"] + [c for c in cols if c in support_DF.columns]
    ps = support_DF[keep].copy()
    ps = ps[ps[have].notna().any(axis=1)]
    return {"samples": ps["sampleID"].astype(str).tolist(), "per_sample": ps.reset_index(drop=True)}


def compute_cage_metrics(support_DF):
    """B11 — CAGE 5'-end support per sample (no-op when no CAGE evidence)."""
    return _support_metrics(support_DF, ["perc_within_cage", "median_dist_cage"])


def compute_polya_metrics(support_DF):
    """B12 — polyA 3'-end support per sample (no-op when no polyA evidence)."""
    return _support_metrics(support_DF, ["perc_within_polya", "perc_polya_motif",
                                         "median_dist_polya"])


def compute_sr_metrics(support_DF):
    """B13 — short-read support per sample: TSS ratio + junction coverage fraction
    (no-op when no short-read evidence)."""
    return _support_metrics(support_DF, ["median_ratio_TSS", "perc_jxn_SR_supported"])


def missing_support_types(cage_metrics, polya_metrics, sr_metrics):
    """Part C — the orthogonal-support evidence types NOT present in this run (their
    B metrics have no populated samples), each paired with the SQANTI3 QC flag(s)
    that would supply it. Returns a list of ``(label, flags)`` tuples; empty when
    every support type is already present."""
    missing = []
    if not (cage_metrics and cage_metrics.get("samples")):
        missing.append(("CAGE 5'-end support (transcription start sites)", "--CAGE_peak"))
    if not (polya_metrics and polya_metrics.get("samples")):
        missing.append(("polyA 3'-end support (cleavage sites)",
                        "--polyA_peak / --polyA_motif_list"))
    if not (sr_metrics and sr_metrics.get("samples")):
        missing.append(("short-read support (junction coverage / TSS ratio)",
                        "--short_reads / --coverage"))
    return missing


def plot_optional_support_note_page(pdf, cage_metrics, polya_metrics, sr_metrics):
    """Part C — a compact end-of-report note naming the orthogonal-support types
    that were not supplied and the (optional) SQANTI3 flags that would add them.

    Drawn only when at least one type is missing. Deliberately non-judgemental:
    this evidence is optional and the report is complete without it — the note
    exists so users who *have* CAGE/polyA/short-read data know they can enrich a
    future run, never to imply the current report is deficient."""
    missing = missing_support_types(cage_metrics, polya_metrics, sr_metrics)
    if not missing:
        return
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.5, 0.9, "Optional orthogonal support", ha="center",
             fontsize=16, weight="bold", color="#333333")
    lines = ["This report is complete as it stands. The independent end / junction",
             "evidence below was not supplied to SQANTI3 and is entirely optional —",
             "providing it would add extra cross-checks of transcript 5'/3' ends and",
             "splice junctions in a future run, but it is not required for any metric",
             "in this report:", ""]
    for label, flags in missing:
        lines.append(f"   •  {label}")
        lines.append(f"         add with:  {flags}")
        lines.append("")
    lines.append("None of these are needed to interpret the QC results above.")
    fig.text(0.12, 0.8, "\n".join(lines), ha="left", va="top",
             fontsize=12, color="#555555")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _summarize_completeness(class_DF, sampleID, exp_factor, exp_factor_val, window=50):
    """Per-sample 5'/3' completeness digest from read-to-gene end distances.

    ``diff_to_gene_TSS``/``diff_to_gene_TTS`` give each read's signed distance to
    the annotated gene 5'/3' end. We keep a compact digest sufficient to redraw a
    cumulative completeness profile (|distance| ECDF) and to compute the
    per-sample "% complete within +/- window bp" scalars — collapsed to a bounded
    histogram of |distance| so memory does not scale with read count.

    Returns [sampleID, exp_factor, end, abs_dist, count], where ``end`` is
    '5prime' or '3prime' and ``abs_dist`` is clipped at a generous ceiling
    (10*window) so a few extreme reads do not blow up the table; the ceiling bin
    is the ">= ceiling" tail. Reads with no computable distance are dropped.
    Robust to the columns being absent (older SQANTI3 output): returns an empty
    typed frame, and the downstream metric/plot code no-ops.
    """
    cols = ["sampleID", exp_factor, "end", "abs_dist", "count"]
    have = all(c in class_DF.columns for c in ("diff_to_gene_TSS", "diff_to_gene_TTS"))
    if class_DF is None or class_DF.empty or not have:
        return pd.DataFrame(columns=cols)

    ceiling = int(10 * window)
    parts = []
    for end, col in (("5prime", "diff_to_gene_TSS"), ("3prime", "diff_to_gene_TTS")):
        d = pd.to_numeric(class_DF[col], errors="coerce").abs().dropna()
        if d.empty:
            continue
        d = d.clip(upper=ceiling).astype(int)
        h = d.value_counts().reset_index()
        h.columns = ["abs_dist", "count"]
        h["end"] = end
        parts.append(h)
    if not parts:
        return pd.DataFrame(columns=cols)
    out = pd.concat(parts, ignore_index=True)
    out["sampleID"] = sampleID
    out[exp_factor] = exp_factor_val
    return out[cols]


def proc_samples(args, design_file, ref):
    # Read design file
    design_DF = pd.read_csv(design_file, sep=",")
    

    
    # Make dictionaries to store each file type
    gene_count_dfs = {}
    ujc_count_dfs = {}
    length_dfs = {}
    err_dfs = {}
    cv_dfs = {}
    jxn_offset_dfs = {}
    completeness_dfs = {}
    site_offset_dfs = {}       # per reference-site median offset (F7 fuzziness concordance)
    fuzz_depth_dfs = {}        # imprecision vs junction depth bin (F4)
    jxn_count_by_sample = {}   # per-sample junction-record count (A4 depth-normalised yield)
    support_dfs = {}           # per-sample orthogonal-support scalars (B11/B12/B13)
    fsm_dfs = {}
    ism_dfs = {}
    nic_dfs = {}
    nnc_dfs = {}
    nov_can_dfs = {}
    length_Dct ={}
    
    
    if args.inFACTOR is None:
        exp_factor = 'temp_factor'
    else:
        exp_factor = args.inFACTOR
    
    
    ##Flag ref genes with at least one mono exonic transcript
    ref_DF = flag_ref_monoexon(ref)
    
    ## CREATE SUMMARY FILES TO MAKE PLOTS FROM
    jxn_cols = ['isoform','chrom','strand','junction_number','genomic_start_coord','genomic_end_coord','junction_category',
                'diff_to_Ref_start_site','diff_to_Ref_end_site','canonical','total_coverage_unique']

    jxn_dtypes = {'isoform':'string', 'chrom': 'string', 'strand': 'string', 'junction_number': 'string',
                  'genomic_start_coord':'Int64', 'genomic_end_coord': 'Int64', 'junction_category':'string',
                  'diff_to_Ref_start_site': 'Int64', 'diff_to_Ref_end_site': 'Int64', 'canonical': 'string',
                  'total_coverage_unique': 'Int64'}   # per-junction read depth (F4); NaN-heavy
    
    class_cols = ['isoform','chrom','strand','exons','associated_gene','associated_transcript','structural_category','subcategory',
                    'length', 'RTS_stage','perc_A_downstream_TTS','ref_length','ref_exons','all_canonical',
                    'diff_to_gene_TSS','diff_to_gene_TTS', "jxn_string", "jxnHash",
                    # B11/B12/B13 orthogonal-support columns (present only when the
                    # matching assay was supplied to SQANTI3; tolerant load skips any
                    # that are absent, and the views no-op when they are all-NaN).
                    'within_CAGE_peak', 'dist_to_CAGE_peak',
                    'within_polyA_site', 'dist_to_polyA_site', 'polyA_motif_found',
                    'ratio_TSS']

    class_dtypes = {'isoform': 'string', 'chrom': 'string', 'strand': 'string', 'exons': 'Int64', 'associated_gene': 'string','associated_transcript': 'string',
                    'structural_category': 'string', 'subcategory': 'string','length': 'Int64', 'RTS_stage': 'boolean', 'perc_A_downstream_TTS': float,
                    'ref_length': 'Int64','ref_exons': 'Int64', 'all_canonical': 'string',
                    'diff_to_gene_TSS': 'Int64', 'diff_to_gene_TTS': 'Int64', 'jxn_string':'string', "jxnHash":'string',
                    'within_CAGE_peak': 'string', 'dist_to_CAGE_peak': 'Float64',
                    'within_polyA_site': 'string', 'dist_to_polyA_site': 'Float64',
                    'polyA_motif_found': 'string', 'ratio_TSS': 'Float64'}
    
    # Per-sample processing is independent (each writes only its own dict keys),
    # so it runs in parallel when --jobs>1 (thread pool; file loading dominates).
    def _one(row):
        # gtf = row['gtf_file']
        sampleID = row['sampleID']
        class_file = row['classification_file']
        jxn_file = row['junction_file']
        
        if exp_factor == 'temp_factor' :
            exp_factor_val = 0
        else:
            exp_factor_val = row[exp_factor]
        
        reads_logger.info("Loading junction file: "+ sampleID)
        jxn_DF = load_sqanti_file(jxn_file, jxn_cols, jxn_dtypes)
        #jxn_DF = load_sqanti_file(jxn_file, jxn_cols)
        
        reads_logger.info("Loading classification file: "+ sampleID)
        # Tolerant: the optional CAGE/polyA/short-read support columns are absent in
        # evidence-free or older SQANTI3 output, and must not break the load.
        class_DF = load_sqanti_file(class_file, class_cols, class_dtypes, tolerant=True)
        #class_DF = load_sqanti_file(class_file, class_cols)
    
        ##Merge in ref DF
        class_DF = merge_dfs(class_DF, ref_DF, 'associated_gene', 'gene_id', 'left')
    
        ##Get number of novel, known canonical and non-canonical junctions
        count_jxns_DF = _summarize_junctions(jxn_DF)

        # Total canonical/novel junctions for this sample
        nov_can_DF = count_jxns_DF.drop('isoform', axis=1).sum().to_frame().transpose()
        nov_can_DF['sampleID'] = sampleID
        nov_can_DF[exp_factor] = exp_factor_val

        ##Merge classification DF and per-isoform junction counts into one
        class_DF = merge_dfs(class_DF, count_jxns_DF, 'isoform', 'isoform')
        categories = list(class_DF['structural_category'].unique())

        gene_count_DF = _summarize_gene_counts(class_DF, sampleID, exp_factor, exp_factor_val)
        ISM_DF = _summarize_subcategory(class_DF, categories, ['incomplete-splice_match'], sampleID, exp_factor, exp_factor_val)
        FSM_DF = _summarize_subcategory(class_DF, categories, ['full-splice_match'], sampleID, exp_factor, exp_factor_val)
        NIC_DF = _summarize_subcategory(class_DF, categories, ['novel_in_catalog'], sampleID, exp_factor, exp_factor_val)
        NNC_DF = _summarize_subcategory(class_DF, categories, ['novel_not_in_catalog'], sampleID, exp_factor, exp_factor_val)
        ujc_count_DF = _summarize_ujc(class_DF, sampleID, exp_factor, exp_factor_val)
        length_summary_DF = _summarize_length(class_DF, sampleID, exp_factor, exp_factor_val)

        ##Length arrays for violin plot
        length_Dct[sampleID] = np.array(class_DF['length'])

        err_DF = _summarize_errors(class_DF, args.cfg()['intrapriming_perc_A_cutoff'], sampleID, exp_factor, exp_factor_val)

        ##Calculate junction cv for each of the samples
        cv_DF = calc_jxn_cv(jxn_DF, class_DF, ref_DF, dropFlag=True)
        cv_DF['sampleID'] = sampleID
        cv_DF[exp_factor] = exp_factor_val

        ##Signed splice-site offset spectrum (fuzziness) for each sample
        jxn_offset_DF = _summarize_jxn_offsets(
            jxn_DF, sampleID, exp_factor, exp_factor_val,
            window=args.cfg().get('jxn_offset_window', 15))

        ##5'/3' completeness digest for each sample
        completeness_DF = _summarize_completeness(
            class_DF, sampleID, exp_factor, exp_factor_val,
            window=args.cfg().get('completeness_window', 50))

        ##Per reference-site median offset (F7 fuzziness concordance)
        site_offset_DF = _summarize_site_offsets(
            jxn_DF, sampleID, exp_factor, exp_factor_val,
            window=args.cfg().get('jxn_offset_window', 15))

        ##Imprecision vs junction read depth (F4); no-op when coverage is absent
        fuzz_depth_DF = _summarize_fuzz_by_depth(
            jxn_DF, sampleID, exp_factor, exp_factor_val,
            window=args.cfg().get('jxn_offset_window', 15))

        ##Orthogonal-support scalars (B11/B12/B13); all-NaN when no evidence supplied
        support_DF = _summarize_support(class_DF, jxn_DF, sampleID, exp_factor, exp_factor_val)

        ##Store all summary dataframes in dictionaries
        gene_count_dfs[sampleID] = gene_count_DF
        ujc_count_dfs[sampleID] = ujc_count_DF
        length_dfs[sampleID] = length_summary_DF
        cv_dfs[sampleID] = cv_DF
        jxn_offset_dfs[sampleID] = jxn_offset_DF
        completeness_dfs[sampleID] = completeness_DF
        site_offset_dfs[sampleID] = site_offset_DF
        fuzz_depth_dfs[sampleID] = fuzz_depth_DF
        jxn_count_by_sample[sampleID] = int(len(jxn_DF))
        support_dfs[sampleID] = support_DF
        err_dfs[sampleID] = err_DF
        fsm_dfs[sampleID] = FSM_DF
        ism_dfs[sampleID] = ISM_DF
        nic_dfs[sampleID] = NIC_DF
        nnc_dfs[sampleID] = NNC_DF
        nov_can_dfs[sampleID] = nov_can_DF
    
        reads_logger.info(sampleID + " done processing")
        del(cv_DF)
        del(class_DF)
        del(jxn_DF)
        reads_logger.debug("Memory cleared for next sample")

    _run_parallel(_one, [row for _, row in design_DF.iterrows()], getattr(args, 'jobs', 1))

    # Restore design order: under --jobs>1 the dicts are filled in completion
    # order, so reorder them by the design so downstream output is deterministic.
    _order = design_DF['sampleID'].tolist()
    _reord = lambda d: {s: d[s] for s in _order if s in d}
    gene_count_dfs, ujc_count_dfs, length_dfs = _reord(gene_count_dfs), _reord(ujc_count_dfs), _reord(length_dfs)
    cv_dfs, err_dfs, fsm_dfs = _reord(cv_dfs), _reord(err_dfs), _reord(fsm_dfs)
    ism_dfs, nov_can_dfs = _reord(ism_dfs), _reord(nov_can_dfs)
    nic_dfs, nnc_dfs = _reord(nic_dfs), _reord(nnc_dfs)
    length_Dct = _reord(length_Dct)
    jxn_offset_dfs = _reord(jxn_offset_dfs)
    completeness_dfs = _reord(completeness_dfs)
    site_offset_dfs = _reord(site_offset_dfs)
    fuzz_depth_dfs = _reord(fuzz_depth_dfs)
    jxn_count_by_sample = _reord(jxn_count_by_sample)
    support_dfs = _reord(support_dfs)

    return(ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_dfs, nnc_dfs, nov_can_dfs,length_Dct, jxn_offset_dfs, completeness_dfs, site_offset_dfs, fuzz_depth_dfs, jxn_count_by_sample, support_dfs )

def _export_reads_tables(args, exp_factor, core_tables, alltables_tables):
    """Write the per-run summary CSVs. The synthetic ``temp_factor`` column is
    dropped in no-factor runs; the ``*_counts`` extras are written only with
    ``--all-tables``. ``*_tables`` are (DataFrame, filename-suffix) pairs."""
    drop_factor = args.inFACTOR is None

    def _write(df, suffix):
        out = df.drop(columns=[exp_factor]) if drop_factor else df
        out.to_csv(os.path.join(args.OUT, args.PREFIX + suffix), index=False)

    for df, suffix in core_tables:
        _write(df, suffix)
    if args.ALLTABLES:
        for df, suffix in alltables_tables:
            _write(df, suffix)


def prep_tables(args, ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_dfs, nnc_dfs, nov_can_dfs,length_Dct, jxn_offset_dfs, completeness_dfs, site_offset_dfs, fuzz_depth_dfs, support_dfs ):
    
    if args.inFACTOR is None:
        exp_factor = 'temp_factor'
    else:
        exp_factor = args.inFACTOR
    
    
    #Combine datframes for all samples
    #Cat together all gene count dfs
    gene_count_DF = pd.concat(gene_count_dfs.values(), sort=False)
    
    # Add flags to gene count DF
  
    gene_count_DF['flag_annotated_gene'] = np.where(gene_count_DF['associated_gene'].isin(ref_DF['gene_id'].unique()), 1, 0)
    
    #Cat together all ujc dfs
    ujc_count_DF = pd.concat(ujc_count_dfs.values())
    
    ujc_count_DF['flag_annotated_gene'] = np.where(ujc_count_DF['associated_gene'].isin(ref_DF['gene_id'].unique()), 1, 0)
    
    ##Cat together all length DFs 
    length_DF = pd.concat(length_dfs.values(), sort=False)
    cols = ['sampleID', exp_factor] + [col for col in length_DF.columns if col not in ['sampleID', exp_factor]]
    length_DF = length_DF[cols]
    
    #Cat together all err DFs
    
    err_DF = pd.concat(err_dfs.values(), sort=False)
    cols = ['sampleID', exp_factor] + [col for col in err_DF.columns if col not in ['sampleID', exp_factor]]
    err_DF = err_DF[cols]
    ##Cat together all cv DFs
    
    cv_DF = pd.concat(cv_dfs.values(), sort=False)

    ##Cat together all junction-offset (fuzziness) DFs
    _offset_frames = [d for d in jxn_offset_dfs.values() if d is not None and not d.empty]
    if _offset_frames:
        jxn_offset_DF = pd.concat(_offset_frames, sort=False, ignore_index=True)
    else:
        jxn_offset_DF = pd.DataFrame(columns=['sampleID', exp_factor, 'site', 'offset', 'canonical', 'count'])

    ##Cat together all 5'/3' completeness DFs
    _comp_frames = [d for d in completeness_dfs.values() if d is not None and not d.empty]
    if _comp_frames:
        completeness_DF = pd.concat(_comp_frames, sort=False, ignore_index=True)
    else:
        completeness_DF = pd.DataFrame(columns=['sampleID', exp_factor, 'end', 'abs_dist', 'count'])

    _site_frames = [d for d in site_offset_dfs.values() if d is not None and not d.empty]
    if _site_frames:
        site_offset_DF = pd.concat(_site_frames, sort=False, ignore_index=True)
    else:
        site_offset_DF = pd.DataFrame(columns=['sampleID', exp_factor, 'ref_key',
                                               'site_type', 'median_offset', 'n'])

    _fuzz_frames = [d for d in fuzz_depth_dfs.values() if d is not None and not d.empty]
    if _fuzz_frames:
        fuzz_depth_DF = pd.concat(_fuzz_frames, sort=False, ignore_index=True)
    else:
        fuzz_depth_DF = pd.DataFrame(columns=['sampleID', exp_factor, 'depth_bin',
                                              'n_sites', 'n_imprecise'])

    ##Cat together the per-sample orthogonal-support rows (B11/B12/B13)
    _support_frames = [d for d in support_dfs.values() if d is not None and not d.empty]
    if _support_frames:
        support_DF = pd.concat(_support_frames, sort=False, ignore_index=True)
    else:
        support_DF = pd.DataFrame(columns=['sampleID', exp_factor] + _SUPPORT_COLS)

    #Cat subcategory DFs
    
    FSM_DF = pd.concat(fsm_dfs.values(), sort=False)
    FSM_DF.fillna(0, inplace=True)
    
    ISM_DF = pd.concat(ism_dfs.values(), sort=False)
    ISM_DF.fillna(0, inplace=True)
    
    NIC_DF = pd.concat(nic_dfs.values(), sort=False)
    NIC_DF.fillna(0, inplace=True)

    NNC_DF = pd.concat(nnc_dfs.values(), sort=False)
    NNC_DF.fillna(0, inplace=True)

    #Cat nov_can_dfs
    nov_can_DF = pd.concat(nov_can_dfs.values(), sort=False)
    nov_can_DF.fillna(0, inplace=True)
        
    ##Export tables
    _export_reads_tables(
        args, exp_factor,
        core_tables=[
            (gene_count_DF, '_gene_counts.csv'),
            (ujc_count_DF, '_ujc_counts.csv'),
            (length_DF, '_length_summary.csv'),
            (cv_DF, '_cv.csv'),
            (jxn_offset_DF, '_jxn_offsets.csv'),
            (completeness_DF, '_completeness.csv'),
            (site_offset_DF, '_site_offsets.csv'),
            (fuzz_depth_DF, '_fuzz_by_depth.csv'),
            (support_DF, '_support.csv'),
        ],
        alltables_tables=[
            (err_DF, '_err_counts.csv'),
            (FSM_DF, '_FSM_counts.csv'),
            (ISM_DF, '_ISM_counts.csv'),
            (NIC_DF, '_NIC_counts.csv'),
            (NNC_DF, '_NNC_counts.csv'),
            (nov_can_DF, '_jxn_counts.csv'),
        ],
    )

    return (gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_DF, NNC_DF, nov_can_DF, length_Dct, jxn_offset_DF, completeness_DF, site_offset_DF, fuzz_depth_DF, support_DF)




def identify_cand_underannot(args, ujc_count_DF, factor_level=None, pdf=None, out_path=None, plot=True):
    """Compute under-annotation tables (writes CSVs) and append the
    under-annotation plot section to the report.

    If ``pdf`` (an already-open ``PdfPages``) is given, pages are appended to
    that shared report; otherwise a standalone PDF is written to ``out_path``.
    When ``plot`` is False, only the CSV tables are written (no plotting) — used
    for HTML-only runs that render the section interactively instead.
    """

    exp_factor = args.inFACTOR
    
    gene_cov_thresh = args.ANNOTEXP
    #gene_cov_thresh = 100
    ujc_perc_cov_thresh = args.PERCCOV
    #ujc_perc_cov_thresh = 20
    max_jxn_thresh = args.PERCMAXJXN
    #max_jxn_thresh = 80
    #factor_level = 'PacBio'
    
    flag_cov_col = f'flag_cov_gt_{ujc_perc_cov_thresh}_perc'
    flag_jxn_col = f'flag_gt_{max_jxn_thresh}_perc_maxjxns'
    
    flag_ujc_in_gene_col = f'flag_ujc_gt_{ujc_perc_cov_thresh}_perc_in_gene'
    flag_fsm_in_gene_col = f'flag_FSM_gt_{ujc_perc_cov_thresh}_perc_in_gene'

    
    def flag_gene_annotated_ujc(group):
        if any(group['structural_category'] == 'full-splice_match'):
            return 1
        else:
            return 0

    def gene_has_well_covered_transcript(group):
        if any(group[flag_cov_col] == 1):
            return 1
        else:
            return 0

    def gene_has_annotated_well_covered_transcript(group):
        if any((group[flag_cov_col] == 1) & (group['structural_category'] == 'full-splice_match')):
            return 1
        else:
            return 0
        
    def categorize_gene(group):
        # Extract the flags
        flag_FSM_in_gene = group['flag_FSM_in_gene'].max()  # max to check if any row has flag_FSM_in_gene == 1
       # flag_ujc_gt_x_perc = group[flag_ujc_in_gene_col].max()
        flag_FSM_gt_x_perc = group[flag_fsm_in_gene_col].max()
        
        well_covered_ujc_novel = group.loc[
        (group[flag_ujc_in_gene_col] == 1) &
        (group['structural_category'].isin(['novel_in_catalog', 'novel_not_in_catalog']))
            ].shape[0] > 0
    
        # Apply the conditions
        if flag_FSM_in_gene == 0 and well_covered_ujc_novel:
            return 'underannotated_with_candidate_transcript'
        elif flag_FSM_in_gene == 0 and not well_covered_ujc_novel:
            return 'underannotated_no_candidate_transcripts'
        elif flag_FSM_in_gene == 1 and flag_FSM_gt_x_perc == 1:
            return 'annotated_with_well_covered_FSM'
        elif flag_FSM_in_gene == 1 and flag_FSM_gt_x_perc == 0:
            return 'annotated_with_low_coverage_FSM'
        else:
            return 'unclassified'
    
    # Calculate total_jxns for each row
    ujc_count_DF['total_jxns'] = ujc_count_DF[['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']].sum(axis=1)
    
    if factor_level is not None and args.inFACTOR is not None:
        ujc_DF = ujc_count_DF[ujc_count_DF[exp_factor] == factor_level]
    else:
        ujc_DF = ujc_count_DF
    
    #Keep only annotated genes
    ujc_DF = ujc_DF[ujc_DF['flag_annotated_gene'] == 1]
    
    #Drop monoexons
    ujc_DF.loc[:,'total_jxns'] = ujc_DF[['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']].sum(axis=1)
    ujc_DF = ujc_DF[ujc_DF['total_jxns'] > 0]
    
    #Get total number of counts across all samples for each UJC
    read_count_sum = ujc_DF.groupby(['jxnHash', 'structural_category','associated_gene', 'total_jxns'])['read_count'].sum().reset_index()
    read_count_sum.columns = ['jxnHash', 'structural_category', 'associated_gene', 'total_jxns', 'total_read_count']

    # Calculate the total coverage of each gene (sum of read_count for all jxnHashes in the same gene)
    gene_coverage_sum = ujc_DF.groupby('associated_gene')['read_count'].sum().reset_index()
    gene_coverage_sum.columns = ['associated_gene', 'total_gene_coverage']
    
    # Get the maximum number of total_jxns for each gene
    max_jxns = ujc_DF.groupby('associated_gene')['total_jxns'].max().reset_index()
    max_jxns.columns = ['associated_gene', 'gene_max_jxns']

    #Merge gene coverage, read counts per UJC and max junctions per gene
    merged_df = read_count_sum.merge(gene_coverage_sum, on='associated_gene', how='left')
    merged_df  = merged_df.merge(max_jxns, on='associated_gene', how='left')
    #Calculate the perecentage of the maxjxns and the perc of the total gene coverage
    merged_df['perc_max_jxns'] = (merged_df['total_jxns'] / merged_df['gene_max_jxns']) * 100
    merged_df['perc_gene_coverage'] = (merged_df['total_read_count'] / merged_df['total_gene_coverage']) * 100

    #Only keep FSMs, NICs and NNCs in genes with at least 100 reads
    merged_df = merged_df[merged_df['total_gene_coverage'] >= gene_cov_thresh]
    merged_df = merged_df[merged_df['structural_category'] != 'incomplete-splice_match'] 
    merged_df = merged_df[merged_df['structural_category'] != 'genic']
    

    merged_df['flag_FSM_in_gene'] = merged_df.groupby('associated_gene')['structural_category'].transform(lambda x: flag_gene_annotated_ujc(merged_df.loc[x.index]))
    #merged_df['flag_underannotated_gene'] = merged_df['flag_FSM_in_gene'].apply(lambda x: 1 if x == 0 else 0)
    
    merged_df[flag_cov_col] = merged_df['perc_gene_coverage'].apply(lambda x: 1 if x > ujc_perc_cov_thresh else 0)
    
    merged_df[flag_ujc_in_gene_col] = merged_df.groupby('associated_gene')[flag_cov_col].transform(lambda x:  gene_has_well_covered_transcript(merged_df.loc[x.index]))
    merged_df[flag_fsm_in_gene_col] = merged_df.groupby('associated_gene')[flag_cov_col].transform(lambda x:   gene_has_annotated_well_covered_transcript(merged_df.loc[x.index]))
    
    merged_df[flag_jxn_col] = merged_df['perc_max_jxns'].apply(lambda x: 1 if x > max_jxn_thresh else 0)
    merged_df['flag_putative_novel_transcript'] = merged_df.apply(
                                                lambda row: 1 if row[flag_cov_col] == 1 and row[flag_jxn_col] == 1 else 0,axis=1)
    
    ## Categorize genes based on gene categories
    if merged_df.empty:
        # No gene passed the under-annotation filters (commonly because the
        # gene_expression cutoff is higher than any gene's coverage). Rather
        # than aborting the whole report, warn, emit empty tables, add a short
        # note page to the report, and let the main QC plots proceed.
        reads_logger.warning(
            f"No genes met the under-annotation criteria (gene_expression={gene_cov_thresh}). "
            f"Skipping under-annotation plots and writing empty tables. "
            f"Lower --gene_expression to include lower-coverage genes."
        )
        empty_summary = pd.DataFrame(columns=['associated_gene', 'gene_category'])
        if 'gene_category' not in merged_df.columns:
            merged_df['gene_category'] = pd.Series(dtype='object')
        if factor_level is not None and args.inFACTOR is not None:
            merged_df.to_csv(os.path.join(args.OUT, args.PREFIX + '_' + factor_level + '_putative_underannotated_transcripts.csv'), index=False)
            empty_summary.to_csv(os.path.join(args.OUT, args.PREFIX + '_' + factor_level + '_gene_classfication.csv'), index=False)
        else:
            merged_df.to_csv(os.path.join(args.OUT, args.PREFIX + '_putative_underannotation.csv'), index=False)
            empty_summary.to_csv(os.path.join(args.OUT, args.PREFIX + '_gene_classfication.csv'), index=False)
        if not plot:
            return
        with (PdfPages(out_path) if pdf is None else nullcontext(pdf)) as pdf:
            divider = plt.figure(figsize=(14, 10))
            divider.text(0.5, 0.5, "Under-annotation analysis", ha='center', va='center', fontsize=26)
            pdf.savefig(divider)
            plt.close(divider)
            fig = plt.figure(figsize=(14, 10))
            fig.text(0.5, 0.5, "No genes met the under-annotation criteria",
                     ha='center', va='center', fontsize=20)
            pdf.savefig(fig)
            plt.close(fig)
        return
    # Group the original DataFrame by associated_gene
    grouped = merged_df.groupby('associated_gene')
    # Apply the categorization function to each group
    gene_categories = grouped.apply(categorize_gene,include_groups=False)

    # Convert the result to a DataFrame
    summary_df = pd.DataFrame({
        'associated_gene': gene_categories.index,
        'gene_category': gene_categories.values
            }).reset_index(drop=True)

    merged_df = merged_df[merged_df['structural_category'] != 'full-splice_match']
    merged_df = merged_df[(merged_df['flag_putative_novel_transcript'] == 1) | (merged_df['flag_FSM_in_gene'] == 0)]
    merged_df = pd.merge(merged_df, summary_df[['associated_gene', 'gene_category']], on='associated_gene', how='left')
    
    if factor_level is not None and args.inFACTOR is not None:
        merged_df.to_csv(os.path.join(args.OUT, args.PREFIX + '_' + factor_level +'_putative_underannotated_transcripts.csv'), index=False)
        summary_df.to_csv(os.path.join(args.OUT, args.PREFIX + '_' + factor_level +'_gene_classfication.csv'), index=False)
    else:
        merged_df.to_csv(os.path.join(args.OUT, args.PREFIX + '_putative_underannotation.csv'), index=False)
        summary_df.to_csv(os.path.join(args.OUT, args.PREFIX + '_gene_classfication.csv'), index=False)

    if not plot:
        return

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['pdf.fonttype'] = 42
    
    with (PdfPages(out_path) if pdf is None else nullcontext(pdf)) as pdf:
        # Section divider page for the under-annotation section of the report
        title_fig = plt.figure(figsize=(14, 10))
        title_fig.text(0.5, 0.5, "Under-annotation analysis", ha='center', va='center', fontsize=26)
        pdf.savefig(title_fig)
        plt.close(title_fig)
        
        # Gene annotation plot
        color_mapping = underannot_palette

        # Count the occurrences of each gene_category
        category_counts = summary_df['gene_category'].value_counts()

        # Assign colors based on the category
        colors = [color_mapping[category] for category in category_counts.index]

        # Create a barplot with specified colors
        plt.figure(figsize=(14, 10))
        bars = category_counts.plot(kind='bar', color=colors)
        #plt.title('Number of Genes in each annotation category')
        plt.xlabel('Gene Category')
        plt.ylabel('Number of Genes')
        plt.xticks(rotation=90, ha="right")
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        # Add the counts on top of each bar
        for bar in bars.patches:
            bars.annotate(format(bar.get_height(), '.0f'),
                          (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                          ha='center', va='center', size=12, xytext=(0, 8),
                          textcoords='offset points')
        plt.tight_layout()
        matplotlib.rcParams['pdf.fonttype'] = 42
        pdf.savefig()
        plt.close()

        #UJC scatterplots
        for gene_category in merged_df['gene_category'].unique():
            # Filter the DataFrame for the current gene category
            df = merged_df[merged_df['gene_category'] == gene_category].copy() # copy is used to avoid SettingWithCopyWarning
            
            # Create a new column to indicate the color based on the thresholds
            df.loc[:,'Putative Unannotated'] = df.apply(
                lambda row: 'Yes' if row['perc_gene_coverage'] > ujc_perc_cov_thresh and row['perc_max_jxns'] > max_jxn_thresh else 'No',
                axis=1
            )
            
            # Create the scatter plot
            plt.figure(figsize=(14, 10))
            scatter_plot = sns.scatterplot(
                data=df,
                x='perc_gene_coverage',
                y='perc_max_jxns',
                hue='Putative Unannotated',
                palette={'Yes': 'green', 'No': 'grey'},
                markers=True
            )
        
            # Set plot title and labels
            plot_title = gene_category.replace('_', ' ')
            scatter_plot.set_title(plot_title)
            scatter_plot.set_xlabel('Percent of Total Gene Coverage')
            scatter_plot.set_ylabel('Percentage of Max Junctions')
        
            # Show the plot
            plt.xlim(0,100)
            plt.ylim(0,100)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            matplotlib.rcParams['pdf.fonttype'] = 42
            # Save the plot to file
            pdf.savefig()
            plt.close()

def _cv_category_summary(cv_df, exp_factor, sample_factor_DF):
    """Per-sample counts of reference-matching / zero-CV / positive-CV junctions.

    Robust to an empty input (no donor/acceptor junctions pass the -je threshold,
    or a dataset with only monoexons): returns one zero-filled numeric row per
    sample so the counts/percentage columns always exist with the right dtype,
    instead of KeyError-ing or breaking the downstream stacked-bar arithmetic.
    """
    cols = ['ref_match', 'cv_0', 'cv_gt_0']
    if cv_df.empty:
        base = sample_factor_DF[['sampleID', exp_factor]].drop_duplicates().reset_index(drop=True)
        for col in cols:
            base[col] = 0
        return base
    summary = cv_df.groupby(['sampleID', exp_factor]).apply(
        lambda x: pd.Series({
            'ref_match': (x['mean_abs_diff'] == 0).sum(),
            'cv_0': ((x['cv'] == 0) & (x['mean_abs_diff'] != 0)).sum(),
            'cv_gt_0': (x['cv'] > 0).sum(),
        }),
        include_groups=False,
    ).reset_index()
    for col in cols:
        if col not in summary.columns:
            summary[col] = 0
    summary.fillna(0, inplace=True)
    return summary


def prep_data_4_plots(args, gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_DF, NNC_DF, nov_can_DF, length_Dct):
    
    abbr_mapping = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
    "antisense": "AS",
    "fusion": "FUS",
    "genic": "GENIC",
    "genic_intron": "GI",
    "intergenic": "INTER"
    }
    
    if args.inFACTOR is None:
        exp_factor = 'temp_factor'
    else:
        exp_factor= args.inFACTOR
    
    exp_factor_DF = gene_count_DF[['sampleID', exp_factor]].drop_duplicates()
    
    ##Prep data for plot 0/1
    total_reads_per_sample = gene_count_DF.groupby('sampleID')['total_read_count'].sum()
    all_gene_percs = {}

    all_categories = [col for col in gene_count_DF.columns if col in abbr_mapping.keys()]
    
    # Loop through each category
    for category in all_categories:
        
        # Sum the category values for annotated genes grouped by sampleID
        all_category_sum_per_sample = gene_count_DF.groupby('sampleID')[category].sum()
        
        # Calculate the percentage of total reads per sample
        all_gene_percs[f'percent_{category}_annotated'] = (all_category_sum_per_sample / total_reads_per_sample) * 100
    
    # Convert the percentages dictionary to a DataFrame 
    all_gene_percs_DF = pd.DataFrame(all_gene_percs)
    all_gene_percs_DF.reset_index(inplace=True)
    
    all_gene_percs_DF.columns = ['sampleID'] + [abbr_mapping[cat] for cat in all_categories]
    all_gene_percs_long_DF = all_gene_percs_DF.melt(id_vars='sampleID', var_name='Category', value_name='Percentage')
    all_gene_percs_long_DF = merge_dfs(all_gene_percs_long_DF,exp_factor_DF, 'sampleID', 'sampleID')
       
    all_gene_percs_pivot_DF = all_gene_percs_long_DF.pivot_table(index=['sampleID', exp_factor], columns='Category', values='Percentage', fill_value=0).reset_index()
   
    #Prep data for PLOT 1 from the gene count DF - Boxplot %reads in detected genes by structural category
    # Identify annotated genes
    annotated_gene_DF = gene_count_DF[gene_count_DF['flag_annotated_gene'] == 1].copy()
    
    # Aggregate total reads per sample
    #total_reads_per_sample = gene_count_DF.groupby('sampleID')['total_read_count'].sum()
    total_reads_in_annotated_gene_DF =  annotated_gene_DF.groupby('sampleID')['total_read_count'].sum()
    
    #Initialize gene perc dict
    annot_gene_percs = {}
    
    # Categories to calculate percentages for, in annotated genes
    annot_categories = ['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog','genic','genic_intron']
    
    # Loop through each category
    for category in annot_categories:
    
        if category in annotated_gene_DF.columns:
            # Sum the category values for annotated genes grouped by sampleID
            annot_category_sum_per_sample = annotated_gene_DF.groupby('sampleID')[category].sum()
            
            # Calculate the percentage of total reads per sample
            annot_gene_percs[f'{category}'] = (annot_category_sum_per_sample / total_reads_in_annotated_gene_DF) * 100
        
    # Convert the percentages dictionary to a DataFrame 
    annot_gene_percs_DF = pd.DataFrame(annot_gene_percs)
    annot_gene_percs_DF.reset_index(inplace=True)
    
    annot_gene_percs_DF.columns = ['sampleID'] + [abbr_mapping[cat] for cat in annot_categories if cat in annotated_gene_DF.columns] 
    
    #prep data for Plot 2 - Vertical scatter plot %reads in annotated genes by structural category
    
    annot_gene_percs_long_DF = annot_gene_percs_DF.melt(id_vars='sampleID', var_name='Category', value_name='Percentage')
    annot_gene_percs_long_DF = merge_dfs(annot_gene_percs_long_DF,exp_factor_DF, 'sampleID', 'sampleID')
    annot_gene_percs_pivot_DF = annot_gene_percs_long_DF.pivot_table(index=['sampleID', exp_factor], columns='Category', values='Percentage', fill_value=0).reset_index()

    #Prep data for plot 3 - Barplot - number of genes per sample, colured by number of reads assigned to gene
    annotated_gene_DF['read_category'] = annotated_gene_DF['total_read_count'].apply(categorize_by_readcount)
    
    gene_agg_DF = annotated_gene_DF.groupby(['sampleID', 'read_category'])['associated_gene'].nunique().unstack(fill_value=0)
    gene_agg_DF = merge_dfs(gene_agg_DF,exp_factor_DF, 'sampleID', 'sampleID')
    
    
    ##Prep data for plot 4 - Barplot - % genes per sample coloured by number of reads in the gene
    
    # Calculate total counts by sampleID to use for percentage calculation
    gene_counts_by_sample = annotated_gene_DF.groupby(['sampleID', 'read_category'])['associated_gene'].nunique()
    total_gene_counts = annotated_gene_DF.groupby('sampleID')['associated_gene'].nunique()
    
    # Calculate percentages
    gene_percs_read_cat_DF = gene_counts_by_sample.div(total_gene_counts, level='sampleID') * 100  # level='sampleID' ensures division is done within each sampleID group
    # Unstack for plotting
    gene_percs_unstacked = gene_percs_read_cat_DF.unstack(fill_value=0)
    gene_percs_unstacked = merge_dfs(gene_percs_unstacked,exp_factor_DF, 'sampleID', 'sampleID')
    
    ##Prep data for plot 5 -  Boxplots Distribution of % structural category (FSM ISM NIC NNC) 
    #Foucsing on annotated genes
    for category in ['full-splice_match','incomplete-splice_match','novel_in_catalog', 'novel_not_in_catalog']:
        annotated_gene_DF[f'percent_{category}'] = (annotated_gene_DF[category] / annotated_gene_DF['total_read_count']) * 100

    
    melted_annotated_gene_DF=  annotated_gene_DF.melt(id_vars=['sampleID'], 
                                         value_vars=['percent_full-splice_match','percent_incomplete-splice_match','percent_novel_in_catalog', 'percent_novel_not_in_catalog'],
                        var_name='category', value_name='percentage')
    
    melted_annotated_gene_DF = merge_dfs(melted_annotated_gene_DF,exp_factor_DF, 'sampleID', 'sampleID')
    
    
    ###UJCs
    
    ##Prep data for  plots 6 -9 - #Barplot - numberand % UJCs per sample, colured by stackby
     # Add flags to ujc count DF
    
    ujc_count_DF['read_category'] = ujc_count_DF['read_count'].apply(categorize_by_readcount)
    annot_ujc_count_DF = ujc_count_DF[ujc_count_DF['flag_annotated_gene'] == 1]
    
    ujc_cnts_dct = {}
    ujc_percs_dct = {}
    
    for stack_by in ['read_category', 'structural_category']:
    # Assume annot_ujc_count_DF and exp_factor_DF are defined and pre-processed as before

        # Group and count unique 'jxnHash' within each 'sampleID' and 'stack_by' category
        ujc_agg_DF = annot_ujc_count_DF.groupby(['sampleID', stack_by])['jxnHash'].nunique().unstack(fill_value=0).reset_index()
        if stack_by == 'structural_category':
           ujc_agg_DF.columns = ['sampleID'] + [abbr_mapping[cat] for cat in ujc_agg_DF.columns[1:] if cat in abbr_mapping]
        
        # Merge with exp_factor_DF if there's additional information needed
        ujc_agg_DF = ujc_agg_DF.merge(exp_factor_DF, on='sampleID', how='left')

        # Calculate total unique 'jxnHash' counts for each 'sampleID'
        total_ujc_counts_DF = annot_ujc_count_DF.groupby('sampleID')['jxnHash'].nunique().reset_index()

        # Create an empty DataFrame for percentages
        ujc_percs_DF = pd.DataFrame()
        
       
    # Iterate over columns to calculate percentages, excluding 'sampleID' and any columns from 'exp_factor_DF'
        for column in ujc_agg_DF.columns:
            if column not in ['sampleID'] + list(exp_factor_DF.columns):
                # Align and divide by total counts using a mapping
                total_counts_map = total_ujc_counts_DF.set_index('sampleID')['jxnHash'].to_dict()
                ujc_agg_DF['total_counts'] = ujc_agg_DF['sampleID'].map(total_counts_map)
                ujc_percs_DF[column] = (ujc_agg_DF[column] / ujc_agg_DF['total_counts']) * 100
                
        ujc_percs_DF['sampleID'] = ujc_agg_DF['sampleID']
        
        
        ujc_total_DF = ujc_agg_DF[['total_counts','sampleID']]
        ujc_agg_DF.drop('total_counts', axis=1, inplace=True)
        ujc_percs_DF = merge_dfs(ujc_percs_DF,exp_factor_DF, 'sampleID', 'sampleID')
       
        # Store results in dictionaries, assuming they were defined earlier
        ujc_cnts_dct[stack_by] = ujc_agg_DF
        ujc_percs_dct[stack_by] = ujc_percs_DF

    ##Lengths
    ##Prep data for plot 11 - Barplot - number of reads caoloured by read count category
    # Calculate counts for each category
    length_DF['reads_lt_1kb'] = length_DF['total_reads'] - length_DF['reads_gt_1kb']
    length_DF['reads_1kb_to_2kb'] = length_DF['reads_gt_1kb'] - length_DF['reads_gt_2kb']
    length_DF['reads_2kb_to_3kb'] = length_DF['reads_gt_2kb'] - length_DF['reads_gt_3kb']
    length_DF['reads_gt_3kb'] = length_DF['reads_gt_3kb']  # This line is not necessary but added for clarity
    
    # Calculate percentages for each category
    length_cols = ['reads_lt_1kb', 'reads_1kb_to_2kb', 'reads_2kb_to_3kb', 'reads_gt_3kb']
    for category in length_cols:
       length_DF[f'{category}_perc'] = (length_DF[category] / length_DF['total_reads']) * 100
    # Aggregating the counts into a new DataFrame suitable for plotting
    length_cnts_agg = length_DF.set_index('sampleID')[length_cols]
    length_cnts_agg = merge_dfs(length_cnts_agg,exp_factor_DF, 'sampleID', 'sampleID')
    
    ##Prep data for plot 12 - Barplot % reads by read count category
    # Calculating percentages for plotting
    percent_length_cols = [f'{col}_perc' for col in length_cols]
    length_percs_agg = length_DF.set_index('sampleID')[percent_length_cols]
    length_percs_agg = merge_dfs(length_percs_agg,exp_factor_DF, 'sampleID', 'sampleID')
    
    #Prep data for plot 13 - % structural category vs %reads greater than 1kb
    length_DF= merge_dfs(length_DF,all_gene_percs_DF, 'sampleID', 'sampleID')
    if args.ALLTABLES:
        length_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_length_summary_w_category_percs.csv'), index=False)
        all_gene_percs_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_structural_category_percs.csv'), index=False)
    
    ##Length for violin plots
    length_DF2 = pd.DataFrame({k: pd.Series(v) for k, v in length_Dct.items()})
    length_DF2 = length_DF2.melt(var_name='sampleID', value_name='length')
    length_DF2 = merge_dfs(length_DF2,exp_factor_DF, 'sampleID', 'sampleID')


    #Make nov_can_percs_DF
    
     
    ##Prep for cv plot(s)
    #cv_DF = cv_DF[cv_DF['count'] >= 3]
    #jxn_exp=3
    jxn_exp = args.JXNEXP
    cv_don_DF = cv_DF[(cv_DF['count'] >= jxn_exp) & (cv_DF['flag_donor'] == 1)]
    cv_acc_DF = cv_DF[(cv_DF['count'] >= jxn_exp) & (cv_DF['flag_acceptor'] == 1)]
    
    if exp_factor == 'temp_factor':
        cv_don_DF = cv_don_DF.copy()
        cv_acc_DF = cv_acc_DF.copy()
        cv_don_DF.loc[:, exp_factor] = cv_don_DF[exp_factor].fillna(0)
        cv_acc_DF.loc[:, exp_factor] = cv_acc_DF[exp_factor].fillna(0)
        
    
    
    cv_don_summary = _cv_category_summary(cv_don_DF, exp_factor, exp_factor_DF)
    cv_acc_summary = _cv_category_summary(cv_acc_DF, exp_factor, exp_factor_DF)
   
    
    cv_acc_percs = cv_acc_summary.copy()
    cv_acc_totals = cv_acc_summary[['ref_match', 'cv_0', 'cv_gt_0']].sum(axis=1)
    # fillna(0) covers samples with zero qualifying junctions (0/0 -> NaN otherwise).
    cv_acc_percs[['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0']] = (
        cv_acc_summary[['ref_match', 'cv_0', 'cv_gt_0']].div(cv_acc_totals, axis=0) * 100).fillna(0)
    cv_acc_percs.drop(columns=['ref_match', 'cv_0', 'cv_gt_0'], inplace=True)


    cv_don_percs = cv_don_summary.copy()
    cv_don_totals = cv_don_summary[['ref_match', 'cv_0', 'cv_gt_0']].sum(axis=1)
    cv_don_percs[['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0']] = (
        cv_don_summary[['ref_match', 'cv_0', 'cv_gt_0']].div(cv_don_totals, axis=0) * 100).fillna(0)
    cv_don_percs.drop(columns=['ref_match', 'cv_0', 'cv_gt_0'], inplace=True)
    
    for sampleID in exp_factor_DF['sampleID']:
        if sampleID not in cv_don_summary['sampleID'].values:
            row = {column: 0 for column in cv_don_summary.columns}
            row['sampleID'] = sampleID  # Set the sampleID for the new row
            row[exp_factor] = exp_factor_DF.loc[exp_factor_DF['sampleID'] == sampleID, exp_factor].iloc[0]
            row = pd.DataFrame([row])
            cv_don_summary = pd.concat([cv_don_summary, row], ignore_index=True)
            reads_logger.info(f"Note: {sampleID} has no donors in annotated genes with > {jxn_exp} reads")
            
        if sampleID not in cv_acc_summary['sampleID'].values:
            row = {column: 0 for column in cv_acc_summary.columns}
            row['sampleID'] = sampleID  # Set the sampleID for the new row
            row[exp_factor] = exp_factor_DF.loc[exp_factor_DF['sampleID'] == sampleID, exp_factor].iloc[0]
            row = pd.DataFrame([row])
            cv_acc_summary = pd.concat([cv_acc_summary, row], ignore_index=True)
            reads_logger.info(f"Note: {sampleID} has no acceptors in annotated genes with > {jxn_exp} reads")
     
    
    if args.ALLTABLES:
       cv_don_summary.to_csv(os.path.join(args.OUT, args.PREFIX + '_cv_don_counts.csv'), index=False) 
       cv_acc_summary.to_csv(os.path.join(args.OUT, args.PREFIX + '_cv_acc_counts.csv'), index=False)
    ##Prep subcategory dataframes
    
    ##FSM
    FSM_DF = FSM_DF.copy()
    
    FSM_perc_DF = FSM_DF.copy()

    # Calculate the sum for the subcategories (excluding 'sampleID' and exp_factor)
    categories = [col for col in FSM_perc_DF.columns if col not in ['sampleID', exp_factor]]
    FSM_perc_DF['total'] = FSM_perc_DF[categories].sum(axis=1)
    
    # Calculate the percentage for each subcategory
    for category in categories:
        FSM_perc_DF[category] = (FSM_perc_DF[category] / FSM_perc_DF['total']) * 100
    FSM_perc_DF.drop('total', axis=1, inplace=True)
    FSM_perc_DF['sampleID'] = FSM_DF['sampleID']
    FSM_perc_DF[exp_factor] = FSM_DF[exp_factor]
    
    ##ISM
    ISM_DF = ISM_DF.copy()
    ISM_perc_DF = ISM_DF.copy()

    # Calculate the sum for the subcategories (excluding 'sampleID' and exp_factor)
    categories = [col for col in ISM_perc_DF.columns if col not in ['sampleID', exp_factor]]
    ISM_perc_DF['total'] = ISM_perc_DF[categories].sum(axis=1)

    # Calculate the percentage for each subcategory
    for category in categories:
        ISM_perc_DF[category] = (ISM_perc_DF[category] / ISM_perc_DF['total']) * 100
    ISM_perc_DF.drop('total', axis=1, inplace=True)
    ISM_perc_DF['sampleID'] = ISM_DF['sampleID']
    ISM_perc_DF[exp_factor] = ISM_DF[exp_factor]

    ##NIC and NNC subcategory percentages (each normalised to 100% within its
    ## structural category, so NIC and NNC read as separate composition plots).
    def _subcat_percent(df):
        out = df.copy()
        cats = [col for col in out.columns if col not in ['sampleID', exp_factor]]
        out['total'] = out[cats].sum(axis=1)
        for category in cats:
            out[category] = (out[category] / out['total']) * 100
        out.drop('total', axis=1, inplace=True)
        out['sampleID'] = df['sampleID']
        out[exp_factor] = df[exp_factor]
        return out

    NIC_DF = NIC_DF.copy()
    NNC_DF = NNC_DF.copy()
    NIC_perc_DF = _subcat_percent(NIC_DF)
    NNC_perc_DF = _subcat_percent(NNC_DF)

    #Make nov_can_perc_DF
    
    nov_can_perc_DF = nov_can_DF.copy()
    categories = [col for col in nov_can_perc_DF.columns if col not in ['sampleID', exp_factor]]
    nov_can_perc_DF['total'] = nov_can_perc_DF[categories].sum(axis=1)
    for category in categories:
        nov_can_perc_DF[category] = (nov_can_perc_DF[category] / nov_can_perc_DF['total']) * 100
    nov_can_perc_DF.drop('total', axis=1, inplace=True)
    nov_can_perc_DF['sampleID'] = nov_can_DF['sampleID']
    nov_can_perc_DF[exp_factor] = nov_can_DF[exp_factor]
    
    ##Prep data for PCA
    pca_err_DF = err_DF.drop(columns=err_DF.filter(regex='^num_reads').columns)
    pca_nov_can_perc_DF = nov_can_perc_DF.rename(columns=lambda x: f'perc_{x}' if x not in ['sampleID', exp_factor] else x)
    pca_cv_don_percs = cv_don_percs.rename(columns=lambda x: f'donors_{x}' if x not in ['sampleID', exp_factor] else x)
    pca_cv_acc_percs = cv_acc_percs.rename(columns=lambda x: f'acceptors_{x}' if x not in ['sampleID', exp_factor] else x)
    pca_ujc_percs_DF = ujc_percs_DF.rename(columns=lambda x: f'perc_ujc_{x}' if x not in ['sampleID', exp_factor] else x)
    pca_ujc_total_DF = ujc_total_DF.rename(columns={'total_counts': 'total_ujcs'})
    pca_length_DF = length_DF.drop(columns=length_DF.filter(regex='^reads').columns)
    pca_length_DF = pca_length_DF.rename(columns=lambda x: f'perc_reads_{x}' if x[0].isupper() else x)
    
    ##Merge all quality metrics for PCA
    quality_DF = pca_err_DF.drop(columns=[exp_factor]).merge(pca_nov_can_perc_DF.drop(columns=[exp_factor]), on='sampleID', how='outer')
    quality_DF = quality_DF.merge(pca_cv_don_percs.drop(columns=[exp_factor]), on='sampleID', how='outer')
    quality_DF = quality_DF.merge(pca_cv_acc_percs.drop(columns=[exp_factor]), on='sampleID', how='outer')
    quality_DF = quality_DF.merge(pca_ujc_percs_DF.drop(columns=[exp_factor]), on='sampleID', how='outer')
    quality_DF = quality_DF.merge(pca_ujc_total_DF, on='sampleID', how='outer')
    quality_DF = quality_DF.merge(pca_length_DF.drop(columns=[exp_factor]), on='sampleID', how='outer')
    quality_DF = quality_DF.merge(exp_factor_DF, on='sampleID', how='outer')
    quality_DF = quality_DF.fillna(0) #!!!
    
    
    pca_features = quality_DF.columns.difference(['sampleID', exp_factor])
    scaler = StandardScaler()
    scaled_feature_DF = scaler.fit_transform(quality_DF[pca_features])
    pca = PCA()
    pca_results = pca.fit_transform(scaled_feature_DF)
    pca_DF = pd.DataFrame(data=pca_results)
    pca_DF['sampleID'] = quality_DF['sampleID']
    pca_DF[exp_factor] = quality_DF[exp_factor]
    #Get PCA loadings
    loadings = pca.components_.T
    loadings_DF = pd.DataFrame(data=loadings, index=pca_features)
    
    #Get variance ratios
    variance_ratio = pca.explained_variance_ratio_
    
    if args.PCATABLES:
        pca_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_pca_values.csv'), index=False)
        loadings_DF2=loadings_DF.reset_index()
        loadings_DF2.to_csv(os.path.join(args.OUT, args.PREFIX + '_pca_loadings.csv'), index=False)
        variance_DF=pd.DataFrame(variance_ratio)
        variance_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_pca_variance.csv'), index=False)

    
    return (all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF, annot_gene_percs_pivot_DF, gene_agg_DF, 
             gene_percs_unstacked, melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF, 
             length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio, 
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_DF, NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_perc_DF, NNC_perc_DF, nov_can_DF, nov_can_perc_DF,
             length_DF2, cv_acc_percs, cv_don_percs)
    
    
    
def render_report_pdf(out_path, all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF, annot_gene_percs_pivot_DF, gene_agg_DF,
             gene_percs_unstacked, melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF,
             length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio,
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_DF, NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_perc_DF, NNC_perc_DF,nov_can_DF, nov_can_perc_DF,
             length_DF2,cv_acc_percs, cv_don_percs, pdf=None, ujc_metrics=None, jxn_offset_metrics=None,
             completeness_metrics=None, scorecard=None, yield_metrics=None, drift_metrics=None,
             tandem_metrics=None, rep_concordance=None, fuzz_concordance=None,
             fuzz_depth_metrics=None, overlap_metrics=None,
             cage_metrics=None, polya_metrics=None, sr_metrics=None, args=None):
    """Render the full SQANTI-reads PDF report.

    One faceted code path serves both modes: with ``--factor`` each page shows
    one panel per factor level; without it, the pipeline's single-level
    ``temp_factor`` placeholder collapses every page to a single panel (and the
    per-facet subtitles are suppressed).
    """

    plt.rcParams.update({'font.size': 13})
    plt.rcParams['pdf.fonttype'] = 42

    exp_factor = args.inFACTOR if args.inFACTOR is not None else 'temp_factor'

    # Stacked-bar pages are rendered by the module-level _render_stacked_bar /
    # _stacked_bars_dict / _stacked_bars_indexed helpers.

    num_factors = all_gene_percs_long_DF[exp_factor].nunique()
    num_samples = all_gene_percs_long_DF['sampleID'].nunique()
    
    # category/subcat palettes and cat_order are module-level constants

    # Suppress warning about too many open figures
    plt.rcParams['figure.max_open_warning'] = 0
    
    #Define sample color palette
    
    unique_sampleIDs = all_gene_percs_long_DF['sampleID'].unique()
    sample_color_palette = {sampleID: sample_seq[i % len(sample_seq)]
                            for i, sampleID in enumerate(unique_sampleIDs)}
    
    with (PdfPages(out_path) if pdf is None else nullcontext(pdf)) as pdf:
        #Cover page
        # Create the title page
        title_fig = plt.figure(figsize=(14, 10))
        title_fig.text(0.5, 0.5, "SQANTI-reads report", ha='center', va='center', fontsize=26)
        pdf.savefig(title_fig)
        plt.close(title_fig)

        # Summary metrics + QC-flag table (same data as the HTML report's top table).
        try:
            from src.utilities.sqanti_reads_report import _compute_summary
            _per_sample, _samples = _compute_summary(
                length_DF, err_DF, all_gene_percs_pivot_DF, gene_agg_DF, args,
                thresholds=args.cfg()['qc_flags'], jxn_offset_metrics=jxn_offset_metrics,
                completeness_metrics=completeness_metrics)
            _render_summary_table_page(pdf, _per_sample, _samples, args.cfg()['qc_flags'])
        except Exception as exc:  # a summary-table hiccup must not sink the report
            reads_logger.warning(f"Could not render the summary-table page: {exc}")

        # Cohort-relative sample-outlier scorecard, right after the summary table
        # (both are top-of-report, sample-level QC overviews). No-ops when disabled
        # (too few samples / no scorable metric).
        try:
            plot_scorecard_page(pdf, scorecard)
        except Exception as exc:
            reads_logger.warning(f"Could not render the sample-outlier scorecard: {exc}")

        # Structural-category composition is shown as faceted stacked bars (below),
        # matching the HTML report. The earlier per-sample stripplot views were removed.

        ##Barplot - structural category % - ALL genes
        _render_stacked_bar(pdf, all_gene_percs_pivot_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=category_color_palette,
                            title="Percentage Reads in Each Structural Category - All Genes",
                            xlabel="Sample ID", ylabel="Percentages",
                            legend_title='Structural Category')

        ##Barplot - structural category % - Annotated genes
        _render_stacked_bar(pdf, annot_gene_percs_pivot_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=category_color_palette,
                            title="Percentage Reads in Each Structural Category - Annotated Genes",
                            xlabel="Sample ID", ylabel="Percentages",
                            legend_title='Structural Category')

             ## Barplot subcategories
        _render_stacked_bar(pdf, FSM_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Number of Reads in Each Sub-Category - FSM ',
                            xlabel='Sample ID', ylabel='Number of reads',
                            legend_title='Structural Category')
        
        _render_stacked_bar(pdf, FSM_perc_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Percentage of FSM Reads in Each Sub-Category ',
                            xlabel='Sample ID', ylabel='Percentage',
                            legend_title='Structural Category')
        
        _render_stacked_bar(pdf, ISM_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Number of Reads in Each Sub-Category - ISM ',
                            xlabel='Sample ID', ylabel='Number of reads',
                            legend_title='Structural Category')
        
        
        _render_stacked_bar(pdf, ISM_perc_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Percentage of ISM Reads in Each Sub-Category ',
                            xlabel='Sample ID', ylabel='Percentage',
                            legend_title='Structural Category')
        
        
        _render_stacked_bar(pdf, NIC_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Number of Reads in Each Sub-Category - NIC ',
                            xlabel='Sample ID', ylabel='Number of reads',
                            legend_title='Structural Category')

        _render_stacked_bar(pdf, NIC_perc_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Percentage of NIC Reads in Each Sub-Category ',
                            xlabel='Sample ID', ylabel='Percentage',
                            legend_title='Structural Category')

        _render_stacked_bar(pdf, NNC_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=nnc_subcat_color_palette,
                            title='Number of Reads in Each Sub-Category - NNC ',
                            xlabel='Sample ID', ylabel='Number of reads',
                            legend_title='Structural Category')

        _render_stacked_bar(pdf, NNC_perc_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=nnc_subcat_color_palette,
                            title='Percentage of NNC Reads in Each Sub-Category ',
                            xlabel='Sample ID', ylabel='Percentage',
                            legend_title='Structural Category')
        
        
        ##Barplot - Genes detected counts
        _render_stacked_bar(pdf, gene_agg_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=[readcount_palette[c] for c in ['100+ reads', '50-100 reads', '11-50 reads', '2-10 reads', '1 read']],
                            categories=['100+ reads', '50-100 reads', '11-50 reads', '2-10 reads', '1 read'],
                            title='Number of Genes Detected',
                            xlabel='Sample ID', ylabel='Number of Genes',
                            legend_title='Read Count')
        
        
        ##Barplot - Genes detected percentages
        _render_stacked_bar(pdf, gene_percs_unstacked, exp_factor=exp_factor,
                            num_factors=num_factors, palette=[readcount_palette[c] for c in ['100+ reads', '50-100 reads', '11-50 reads', '2-10 reads', '1 read']],
                            categories=['100+ reads', '50-100 reads', '11-50 reads', '2-10 reads', '1 read'],
                            title='Number of Genes Detected',
                            xlabel='Sample ID', ylabel='Number of Genes',
                            legend_title='Read Count')
    
            
        
            
        # Boxplots Distribution of % structural category (FSM ISM NIC NNC) 
        melted_annotated_gene_DF = melted_annotated_gene_DF.sort_values(by=[exp_factor, 'sampleID'])
        palette = sample_color_palette
        for category in melted_annotated_gene_DF['category'].unique():
            plt.figure(figsize=(14, 10))  
            for i, exp_factor_val in enumerate(melted_annotated_gene_DF[exp_factor].unique(), start=1):
                # Filter the DataFrame for both category and exp_factor
                df_filtered = melted_annotated_gene_DF[(melted_annotated_gene_DF['category'] == category) &
                                                       (melted_annotated_gene_DF[exp_factor] == exp_factor_val)]
                
                # Create a subplot for this exp_factor
                plt.subplot(1, num_factors, i)
                sns.violinplot(data=df_filtered, x='sampleID', y='percentage',
                               hue='sampleID', palette=palette, legend=False)
                if exp_factor != 'temp_factor':
                    plt.title(exp_factor + f" = {exp_factor_val}")
                plt.xticks(rotation=90)
            
            title = plt.suptitle(f'Gene distribution - {category}', y= 1.02)
            plt.tight_layout()
            matplotlib.rcParams['pdf.fonttype'] = 42
            pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
            plt.close()
        
    
    
        # UJC barplots - percent and counts, structural 
        for stack_by in ['read_category', 'structural_category']:
            ujc_cnts_dct[stack_by] = ujc_cnts_dct[stack_by].sort_values(by=[exp_factor, 'sampleID'])
            ujc_percs_dct[stack_by] = ujc_percs_dct[stack_by].sort_values(by=[exp_factor, 'sampleID'])

            if stack_by == 'read_category':
                readcount_cats = ['100+ reads', '50-100 reads', '11-50 reads', '2-10 reads', '1 read']
                readcount_pal = [readcount_palette[c] for c in readcount_cats]
                _render_stacked_bar(pdf, ujc_cnts_dct[stack_by], exp_factor=exp_factor,
                                    num_factors=num_factors, categories=readcount_cats,
                                    palette=readcount_pal, title="Number of UJCs Detected",
                                    xlabel="Sample ID", ylabel="Number of UJCs",
                                    legend_title='Read count')
                _render_stacked_bar(pdf, ujc_percs_dct[stack_by], exp_factor=exp_factor,
                                    num_factors=num_factors, categories=readcount_cats,
                                    palette=readcount_pal, title="UJCs Detected",
                                    xlabel="Sample ID", ylabel="Percentage of UJCs",
                                    legend_title='Read count')
            elif stack_by == 'structural_category':
                _render_stacked_bar(pdf, ujc_cnts_dct[stack_by], exp_factor=exp_factor,
                                    num_factors=num_factors, palette=category_color_palette,
                                    title="Number of UJCs detected",
                                    xlabel="Sample ID", ylabel="Number of UJCs",
                                    legend_title='Structural Category')
                _render_stacked_bar(pdf, ujc_percs_dct[stack_by], exp_factor=exp_factor,
                                    num_factors=num_factors, palette=category_color_palette,
                                    title="UJCs detected",
                                    xlabel="Sample ID", ylabel="Percentage of UJCs",
                                    legend_title='Structural Category')


        ## Total Mapped reads vs % reads gt 1kb
    
        plt.figure(figsize=(14, 10))
        _scatter_labeled(length_DF, 'perc_reads_gt_1kb', 'total_reads', 'sampleID',
                         'Total Reads vs Percentage of Reads > 1kb',
                         'Percentage of Reads > 1kb', 'Total Reads', factor_col=exp_factor)
        plt.tight_layout()
        matplotlib.rcParams['pdf.fonttype'] = 42
        pdf.savefig()
        plt.close()
        
        ## Bar graph read categories - counts
        _render_stacked_bar(pdf, length_cnts_agg, exp_factor=exp_factor,
                            num_factors=num_factors, palette=[length_palette[c] for c in ['reads_lt_1kb', 'reads_1kb_to_2kb', 'reads_2kb_to_3kb', 'reads_gt_3kb']],
                            categories=['reads_lt_1kb', 'reads_1kb_to_2kb', 'reads_2kb_to_3kb', 'reads_gt_3kb'],
                            title='Lengths of All Mapped Reads',
                            xlabel='Sample ID', ylabel='Number of Reads',
                            legend_title='Read Count')
         
        ## Bar graph read categories - %
        _render_stacked_bar(pdf, length_percs_agg, exp_factor=exp_factor,
                            num_factors=num_factors, palette=[length_palette[c] for c in ['reads_lt_1kb_perc', 'reads_1kb_to_2kb_perc', 'reads_2kb_to_3kb_perc', 'reads_gt_3kb_perc']],
                            categories=['reads_lt_1kb_perc', 'reads_1kb_to_2kb_perc', 'reads_2kb_to_3kb_perc', 'reads_gt_3kb_perc'],
                            title='Lengths of All Mapped Reads',
                            xlabel='Sample ID', ylabel='Percentage',
                            legend_title='Read Count')
        
             
        
        
        g = sns.FacetGrid(length_DF2, col=exp_factor, col_wrap=num_factors, height=8, aspect=0.7, sharex=False, sharey=True)
        g.map_dataframe(sns.violinplot, x='sampleID', y='length', hue='sampleID',
                        palette=sample_color_palette, legend=False)
        g.set_axis_labels("Sample ID", "Length")
        g.set_titles("" if exp_factor == 'temp_factor' else exp_factor + " = {col_name}")
        title = g.fig.suptitle("Read Length Distribution", y=1.02, fontsize=20)
        for ax in g.axes.flatten():
            # Pin tick locations before labelling (avoids the FixedFormatter warning).
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.tight_layout()
        matplotlib.rcParams['pdf.fonttype'] = 42
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.close()

        ## Scatterplot % reads > 1kb vs % structural category
        categories =[cat for cat in ['FSM', 'ISM', 'NIC', 'NNC','GI','GENIC'] if cat in length_DF.columns]
        for category in categories:
            plt.figure(figsize=(14, 10))
            _scatter_labeled(length_DF, category, 'perc_reads_gt_1kb', 'sampleID',
                             'Percentage of Reads > 1kb vs %' + category,
                             '% ' + category, 'Percentage of Reads > 1kb', factor_col=exp_factor)
            plt.tight_layout()
            matplotlib.rcParams['pdf.fonttype'] = 42
            pdf.savefig()
            plt.close()
            
        #Bar graph RTS/intra priming
        _render_grouped_bar(pdf, err_DF, exp_factor=exp_factor, num_factors=num_factors,
                            categories=['num_reads_RTS', 'num_reads_intrapriming', 'num_reads_non-canonical'],
                            palette=list(three_series_palette),
                            title="Number of Artefact Reads",
                            xlabel="Sample ID", ylabel="Number of Reads")

        _render_grouped_bar(pdf, err_DF, exp_factor=exp_factor, num_factors=num_factors,
                            categories=['perc_reads_RTS', 'perc_reads_intrapriming', 'perc_reads_non-canonical'],
                            palette=list(three_series_palette),
                            title="Percentage of Artefact Reads",
                            xlabel="Sample ID", ylabel="Percentage of Reads",
                            aspect=0.7)

        # Per-sample cohort-context pages for the RT-switching and intra-priming
        # rates (the artefact metrics that feed the outlier scorecard): each rate
        # against the cohort median, config warn/fail lines, and scorecard flags.
        try:
            plot_artefact_metric_pages(pdf, err_DF, args.cfg().get('qc_flags', {}),
                                       scorecard=scorecard)
        except Exception as exc:
            reads_logger.warning(f"Could not render artefact metric pages: {exc}")

        # Derived read-quality scalars in cohort context: ISM fragment fraction
        # and novel non-canonical junction burden (both feed the scorecard).
        try:
            plot_quality_metric_pages(pdf, ISM_DF, nov_can_DF,
                                      args.cfg().get('qc_flags', {}), scorecard=scorecard)
        except Exception as exc:
            reads_logger.warning(f"Could not render quality metric pages: {exc}")
        
        
        plt.figure(figsize=(14, 10))
        _scatter_labeled(pca_DF, 0, 1, 'sampleID', 'PCA Plot Based on QC metrics',
                         'Principal Component 1', 'Principal Component 2', factor_col=exp_factor)
        title = plt.gca().title
        plt.subplots_adjust(top=0.85, right=0.8)
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
        # Calculate the cumulative variance and determine the number of components to use
        cumulative_variance = np.cumsum(variance_ratio)
        n_components = np.argmax(cumulative_variance >= args.cfg()['pca_cumulative_variance']) + 1
        
         # Create the plots
        fig, ax = plt.subplots(2, 2, figsize=(20, 20), sharex='col', gridspec_kw={'width_ratios': [10, 3], 'height_ratios': [3, 10]})
        loadings_DF = loadings_DF.iloc[:, :n_components]
        link = linkage(loadings_DF, method='average')
        sorted_idx = leaves_list(link)
        loadings_DF = loadings_DF.iloc[sorted_idx]
        
        # Bar plot for explained variance (Scree Plot)
        x_tick_pos = [i + 0.5 for i in range(n_components)]
        ax[0, 0].bar(x_tick_pos, variance_ratio[:n_components], align='center', label='Individual explained variance')
        ax[0, 0].step(x_tick_pos, cumulative_variance[:n_components], where='mid', label='Cumulative explained variance')
        ax[0, 0].set_xticks(x_tick_pos)
        ax[0, 0].set_xticklabels([])  # Clear x tick labels here
        ax[0, 0].set_ylabel('Variance Explained')
        
        # Set x tick labels for the heatmap
        x_ticks = [f'PC{i+1}' for i in range(n_components)]
        sns.heatmap(loadings_DF, cmap="coolwarm", ax=ax[1, 0], cbar_ax=ax[1, 1], xticklabels=x_ticks)
        ax[1, 0].set_yticks(np.arange(loadings_DF.shape[0]) + 0.5)
        ax[1, 0].set_yticklabels(loadings_DF.index, rotation=0)
        ax[1, 0].set_xlabel('Principal Components')
        
        # Use the ax[0,1] for legend
        ax[0, 1].axis('off')  # Turn off the axis lines and labels
        handles, labels = ax[0, 0].get_legend_handles_labels()
        ax[0, 1].legend(handles, labels, loc='center')  # Place legend at the center of ax[0, 1]
        
        title = fig.suptitle('Variance and Heatmap of PC loadings',y=1.02, fontsize=20)
        plt.tight_layout()
        matplotlib.rcParams['pdf.fonttype'] = 42
        _vectorize_colorbars(fig)
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.close()
        
        
        
        _render_stacked_bar(pdf, nov_can_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=[jxn_palette[c] for c in ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']],
                            categories=['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical'],
                            title='Junctions by Category',
                            xlabel='Sample ID', ylabel='Number of Junctions',
                            height=9)
        
        
        _render_stacked_bar(pdf, nov_can_perc_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=[jxn_palette[c] for c in ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']],
                            categories=['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical'],
                            title='Junctions by Category',
                            xlabel='Sample ID', ylabel='Percentage',
                            height=9)
        
        _render_stacked_bar(pdf, cv_acc_summary, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['ref_match', 'cv_0', 'cv_gt_0'],
                            legend_labels=cv_count_legend_labels,
                            title='Number of Detected Acceptors',
                            xlabel='Sample ID', ylabel='Number of Detected Acceptors',)

        _render_stacked_bar(pdf, cv_don_summary, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['ref_match', 'cv_0', 'cv_gt_0'],
                            legend_labels=cv_count_legend_labels,
                            title='Number of Detected Donors',
                            xlabel='Sample ID', ylabel='Number of Detected Donors',)


        categories = ['perc_ref_match','perc_cv_0','perc_cv_gt_0']
        _render_stacked_bar(pdf, cv_don_percs, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0'],
                            legend_labels=cv_perc_legend_labels,
                            title='Percentage of Detected Donors',
                            xlabel='Sample ID', ylabel='Percentage of Detected Donors',)

        categories = ['perc_ref_match','perc_cv_0','perc_cv_gt_0']
        _render_stacked_bar(pdf, cv_acc_percs, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0'],
                            legend_labels=cv_perc_legend_labels,
                            title='Percentage of Detected Acceptors',
                            xlabel='Sample ID', ylabel='Percentage of Detected Acceptors',)

        # UJC saturation, replicate concordance and UpSet plots
        plot_ujc_metrics_pages(pdf, ujc_metrics, factor_col=exp_factor)

        # A6 — pairwise UJC overlap-index heatmap (between-sample repertoire sharing)
        try:
            plot_ujc_overlap_page(pdf, overlap_metrics)
        except Exception as e:
            reads_logger.warning(f"UJC overlap page skipped: {e}")

        # Splice-site fuzziness: offset spectrum, precision profile, canonical split
        plot_jxn_offset_pages(pdf, jxn_offset_metrics)

        # 5'/3' read-end completeness profiles
        plot_completeness_pages(pdf, completeness_metrics)

        # Between-sample comparison views
        plot_jxn_yield_page(pdf, yield_metrics, scorecard=scorecard)
        plot_composition_drift_page(pdf, drift_metrics)
        plot_tandem_sites_page(pdf, tandem_metrics)
        plot_replicate_concordance_page(pdf, rep_concordance)
        plot_fuzziness_concordance_page(pdf, fuzz_concordance)
        plot_fuzz_depth_page(pdf, fuzz_depth_metrics)

        # B11/B12/B13 orthogonal-support views (no-op without the matching evidence)
        try:
            plot_support_pages(pdf, cage_metrics, polya_metrics, sr_metrics, scorecard=scorecard)
        except Exception as e:
            reads_logger.warning(f"orthogonal-support pages skipped: {e}")


def run_reads_plots(
    ref_gtf: str,
    design_file: str,
    out_dir: str,
    prefix: str,
    factor: Optional[str] = None,
    gene_expression: int = 100,
    jxn_expression: int = 10,
    perc_coverage: int = 20,
    perc_junctions: int = 80,
    factor_level: Optional[str] = None,
    all_tables: bool = False,
    pca_tables: bool = False,
    report: str = 'pdf',
    config: Optional[dict] = None,
    jobs: int = 1
):
    """
    Run the reads plotting pipeline.
    
    This function can be called directly from other modules instead of using subprocess.
    
    Parameters:
    -----------
    ref_gtf : str
        Path to the reference GTF file.
    design_file : str
        Path to the design file (must have sampleID column).
    out_dir : str
        Output directory for saving plots and tables.
    prefix : str
        Output filename prefix.
    factor : str, optional
        Experimental factor column in design file to use for faceting plots.
    gene_expression : int, default=100
        Expression cutoff level for determining underannotated genes.
    jxn_expression : int, default=10
        Coverage threshold for detected reference donors and acceptors.
    perc_coverage : int, default=20
        Percent gene coverage of UJC for determining well-covered transcripts.
    perc_junctions : int, default=80
        Percent of the max junctions in gene for determining near full-length putative novel transcripts.
    factor_level : str, optional
        Factor level to evaluate for underannotation.
    all_tables : bool, default=False
        Export all output tables.
    pca_tables : bool, default=False
        Export table for making PCA plots.
    report : str, default='pdf'
        Report format: 'pdf', 'html', or 'both'.
    """
    # Create args object (threaded explicitly through the pipeline)
    args = ReadsPlotArgs(
        inREF=ref_gtf,
        inDESIGN=design_file,
        OUT=out_dir,
        PREFIX=prefix,
        inFACTOR=factor,
        ANNOTEXP=gene_expression,
        JXNEXP=jxn_expression,
        PERCCOV=perc_coverage,
        PERCMAXJXN=perc_junctions,
        FACTORLVL=factor_level,
        ALLTABLES=all_tables,
        PCATABLES=pca_tables,
        report=report,
        config=config,
        jobs=jobs
    )

    reads_logger.info("Starting SQANTI-reads tables and plots generation...")

    # Run the main pipeline
    main(args)

    reads_logger.info("SQANTI-reads tables and plots generation completed.")


def main(args):
    ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_dfs, nnc_dfs, nov_can_dfs,length_Dct, jxn_offset_dfs, completeness_dfs, site_offset_dfs, fuzz_depth_dfs, jxn_count_by_sample, support_dfs = proc_samples(args, args.inDESIGN, args.inREF)

    gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF,FSM_DF, ISM_DF, NIC_DF, NNC_DF, nov_can_DF, length_Dct, jxn_offset_DF, completeness_DF, site_offset_DF, fuzz_depth_DF, support_DF = prep_tables(args, ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs,
                                                                                                                            fsm_dfs, ism_dfs, nic_dfs, nnc_dfs, nov_can_dfs,length_Dct, jxn_offset_dfs, completeness_dfs, site_offset_dfs, fuzz_depth_dfs, support_dfs)
    dfs_for_plotting = prep_data_4_plots( args, gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_DF, NNC_DF, nov_can_DF, length_Dct )

    # UJC-level metrics (saturation, replicate concordance, UpSet) shared by both
    # renderers. Computed before identify_cand_underannot mutates ujc_count_DF.
    ujc_metrics = compute_ujc_metrics(ujc_count_DF, factor_col=args.inFACTOR)

    # A6 pairwise UJC overlap index (asymmetric |A∩B|/|A|) — between-sample view.
    overlap_metrics = compute_ujc_overlap(ujc_metrics)

    # B11/B12/B13 orthogonal-support views (CAGE 5', polyA 3', short-read). All
    # no-op when the matching evidence wasn't supplied to SQANTI3 (columns all-NaN).
    cage_metrics = compute_cage_metrics(support_DF)
    polya_metrics = compute_polya_metrics(support_DF)
    sr_metrics = compute_sr_metrics(support_DF)

    # Splice-site fuzziness metrics (offset spectrum, precision profile,
    # directionality, per-sample imprecision scalar) shared by both renderers.
    jxn_offset_metrics = compute_jxn_offset_metrics(
        jxn_offset_DF, window=args.cfg().get('jxn_offset_window', 15))

    # 5'/3' completeness metrics (read-end ECDF profiles + per-sample within-window
    # scalars) shared by both renderers.
    completeness_metrics = compute_completeness_metrics(
        completeness_DF, window=args.cfg().get('completeness_window', 50))

    # Between-sample comparison metrics (A4 depth-normalised junction yield):
    # junctions per 1000 reads, so a low value is low complexity rather than
    # shallow sequencing. Read counts come from the length summary.
    read_counts = dict(zip(length_DF["sampleID"].astype(str), length_DF["total_reads"]))
    yield_metrics = compute_jxn_yield_vs_depth(jxn_count_by_sample, read_counts)

    # A3 composition drift: pairwise Jensen–Shannon distance between samples'
    # structural-category composition (all_gene_percs_pivot_DF = dfs_for_plotting[2]).
    drift_metrics = compute_composition_drift(dfs_for_plotting[2])

    # F6 tandem splice sites (NAGNAG): ±3/±4/±6 bp excess in the offset spectrum.
    tandem_metrics = compute_tandem_sites(jxn_offset_metrics)

    # Within-group replicate views (A5). factor_map = {sampleID: level}; None when
    # no --factor, which disables A5/F7 (they only mean something for replicates).
    factor_map = None
    if args.inFACTOR is not None and args.inFACTOR in length_DF.columns:
        factor_map = dict(zip(length_DF["sampleID"].astype(str),
                              length_DF[args.inFACTOR].astype(str)))
    rep_concordance = compute_replicate_concordance(
        dfs_for_plotting[2], length_DF, jxn_offset_metrics, factor_map)
    # F7 replicate concordance of splice-site precision (per reference site).
    fuzz_concordance = compute_fuzziness_concordance(site_offset_DF, factor_map)

    # F4 splice-site imprecision vs junction depth (no-op without per-junction coverage).
    fuzz_depth_metrics = compute_fuzz_depth_metrics(fuzz_depth_DF)

    # Cohort-relative sample-outlier scorecard: robust z of each read-QC metric
    # against the cohort, flagging samples that diverge from their peers. Uses the
    # read-level structural composition (all_gene_percs_pivot_DF = dfs_for_plotting[2])
    # and the junction-category counts already computed above.
    scorecard = compute_sample_scorecard(
        assemble_scorecard_metrics(
            length_DF["sampleID"].astype(str).tolist(),
            length_DF=length_DF, err_DF=err_DF,
            all_gene_percs_pivot_DF=dfs_for_plotting[2], nov_can_DF=nov_can_DF,
            completeness_metrics=completeness_metrics,
            jxn_offset_metrics=jxn_offset_metrics, ISM_DF=ISM_DF,
            yield_metrics=yield_metrics, drift_metrics=drift_metrics,
            support_DF=support_DF),
        args.cfg())

    need_pdf = args.report in ("pdf", "both")
    need_html = args.report in ("html", "both")

    report_pdf = os.path.join(args.OUT, args.PREFIX + '_report.pdf')

    # Both renderers now work on copies of the plotting DataFrames, so HTML/PDF
    # order no longer matters for correctness. HTML is still built first only so
    # its under-annotation CSV tables (written here with plot=False) exist before
    # the HTML report reads them.
    if need_html:
        identify_cand_underannot(args, ujc_count_DF, factor_level=args.FACTORLVL, plot=False)
        from src.utilities.sqanti_reads_report import build_html_report
        report_html = os.path.join(args.OUT, args.PREFIX + '_report.html')
        build_html_report(report_html, dfs_for_plotting, args, ujc_metrics=ujc_metrics,
                          jxn_offset_metrics=jxn_offset_metrics,
                          completeness_metrics=completeness_metrics, scorecard=scorecard,
                          yield_metrics=yield_metrics, drift_metrics=drift_metrics,
                          tandem_metrics=tandem_metrics, rep_concordance=rep_concordance,
                          fuzz_concordance=fuzz_concordance, fuzz_depth_metrics=fuzz_depth_metrics,
                          overlap_metrics=overlap_metrics, cage_metrics=cage_metrics,
                          polya_metrics=polya_metrics, sr_metrics=sr_metrics)

    # Single unified PDF report: QC plots followed by the under-annotation
    # section, all in one PDF. identify_cand_underannot also (re)writes its CSVs;
    # passing the open `pdf` makes it append pages here.
    if need_pdf:
        with PdfPages(report_pdf) as pdf:
            render_report_pdf(report_pdf, *dfs_for_plotting, pdf=pdf, ujc_metrics=ujc_metrics,
                              jxn_offset_metrics=jxn_offset_metrics,
                              completeness_metrics=completeness_metrics, scorecard=scorecard,
                              yield_metrics=yield_metrics, drift_metrics=drift_metrics,
                              tandem_metrics=tandem_metrics, rep_concordance=rep_concordance,
                              fuzz_concordance=fuzz_concordance, fuzz_depth_metrics=fuzz_depth_metrics,
                              overlap_metrics=overlap_metrics, cage_metrics=cage_metrics,
                              polya_metrics=polya_metrics, sr_metrics=sr_metrics, args=args)
            identify_cand_underannot(args, ujc_count_DF, factor_level=args.FACTORLVL, pdf=pdf)
            # Part C — closing "optional orthogonal support" note (only when some
            # support type is absent). Placed last, after the under-annotation section.
            plot_optional_support_note_page(pdf, cage_metrics, polya_metrics, sr_metrics)

    # Close all remaining figures to free memory
    plt.close('all')
