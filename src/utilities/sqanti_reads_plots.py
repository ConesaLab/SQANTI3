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

cat_order = ["FSM", "ISM", "NIC", "NNC", "AS", "FUS", "GENIC", "GI", "INTER"]
cat_order_stacked = ["INTER", "GI", "GENIC", "FUS", "AS", "NNC", "NIC", "ISM", "FSM"]

# Junction classification (known/novel × canonical/non-canonical).
jxn_palette = {
    "known_canonical": "#2c7fb8",
    "known_non_canonical": "#7fcdbb",
    "novel_canonical": "#f03b20",
    "novel_non_canonical": "#feb24c",
}

# Three-series plots (green, orange, yellow) — artefacts (RTS/intra-priming/
# non-canonical) and donor/acceptor CV (ref_match/cv_0/cv_gt_0). Applied to the
# series in column order.
three_series_palette = ["#2CA02C", "#FF7F0E", "#FFC20A"]

# Read-count bins — explicit color per bin so it is stable across datasets and
# identical in the PDF and HTML (values match the HTML default sequence order).
readcount_palette = {
    "100+ reads": "#2E91E5",
    "50-100 reads": "#E15F99",
    "11-50 reads": "#1CA71C",
    "2-10 reads": "#FB0D0D",
    "1 read": "#DA16FF",
}

# Read-length bins (both the count columns and their `_perc` variants).
_length_bins = ["reads_lt_1kb", "reads_1kb_to_2kb", "reads_2kb_to_3kb", "reads_gt_3kb"]
_length_colors = ["#2E91E5", "#E15F99", "#1CA71C", "#FB0D0D"]
length_palette = {b: c for b, c in zip(_length_bins, _length_colors)}
length_palette.update({b + "_perc": c for b, c in zip(_length_bins, _length_colors)})

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


def _stacked_bars_indexed(*args, categories, palette, **kwargs):
    """FacetGrid ``map_dataframe`` callback: stack a fixed ``categories`` list in
    order, colored by the aligned ``palette`` list. A category missing from a
    facet subset is treated as all-zero so stack order/colors stay stable."""
    data = kwargs.pop('data')
    ax = plt.gca()
    bottom = np.zeros(len(data))
    for idx, category in enumerate(categories):
        values = data[category].values if category in data.columns else np.zeros(len(data))
        non_zero_indices = values != 0
        ax.bar(data['sampleID'][non_zero_indices], values[non_zero_indices],
               bottom=bottom[non_zero_indices], color=palette[idx], label=category)
        bottom += values


def _render_stacked_bar(pdf, df, *, exp_factor, num_factors, title, xlabel, ylabel,
                        palette, categories=None, legend_title=None,
                        height=8, aspect=1.3, sort=True):
    """Render one faceted stacked-bar page (one panel per ``exp_factor`` level).

    ``palette`` is either a {column: hex} dict (columns inferred + stacked in
    column order) or an ordered list aligned to ``categories`` (fixed stack).
    With the synthetic single-level ``temp_factor`` this collapses to one panel.
    """
    plot_df = df.sort_values(by=[exp_factor, 'sampleID']) if sort else df
    g = sns.FacetGrid(plot_df, col=exp_factor, col_wrap=num_factors, height=height,
                      aspect=aspect, sharex=False, sharey=True)
    if isinstance(palette, dict):
        g.map_dataframe(_stacked_bars_dict, color_palette=palette, exp_factor=exp_factor)
    else:
        g.map_dataframe(_stacked_bars_indexed, categories=categories, palette=palette)
    for ax, (_name, group) in zip(g.axes.flatten(), plot_df.groupby(exp_factor)):
        # Pin the tick locations before labelling them (avoids matplotlib's
        # "FixedFormatter without FixedLocator" warning).
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

def load_sqanti_file(file,col_Lst, dtype_Dct):
    return pd.read_csv(file, sep="\t", usecols=col_Lst, dtype=dtype_Dct, low_memory=True)


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


def proc_samples(args, design_file, ref):
    # Read design file
    design_DF = pd.read_csv(design_file, sep=",")
    

    
    # Make dictionaries to store each file type
    gene_count_dfs = {}
    ujc_count_dfs = {}
    length_dfs = {}
    err_dfs = {}
    cv_dfs = {}
    fsm_dfs = {}
    ism_dfs = {}
    nic_nnc_dfs = {}
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
                'diff_to_Ref_start_site','diff_to_Ref_end_site','canonical']
    
    jxn_dtypes = {'isoform':'string', 'chrom': 'string', 'strand': 'string', 'junction_number': 'string',
                  'genomic_start_coord':'Int64', 'genomic_end_coord': 'Int64', 'junction_category':'string',
                  'diff_to_Ref_start_site': 'Int64', 'diff_to_Ref_end_site': 'Int64', 'canonical': 'string'}
    
    class_cols = ['isoform','chrom','strand','exons','associated_gene','associated_transcript','structural_category','subcategory',
                    'length', 'RTS_stage','perc_A_downstream_TTS','ref_length','ref_exons','all_canonical', "jxn_string", "jxnHash"]

    class_dtypes = {'isoform': 'string', 'chrom': 'string', 'strand': 'string', 'exons': 'Int64', 'associated_gene': 'string','associated_transcript': 'string', 
                    'structural_category': 'string', 'subcategory': 'string','length': 'Int64', 'RTS_stage': 'boolean', 'perc_A_downstream_TTS': float, 
                    'ref_length': 'Int64','ref_exons': 'Int64', 'all_canonical': 'string', 'jxn_string':'string', "jxnHash":'string'}
    
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
        class_DF = load_sqanti_file(class_file, class_cols, class_dtypes)
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
        NIC_NNC_DF = _summarize_subcategory(class_DF, categories, ['novel_in_catalog', 'novel_not_in_catalog'], sampleID, exp_factor, exp_factor_val)
        ujc_count_DF = _summarize_ujc(class_DF, sampleID, exp_factor, exp_factor_val)
        length_summary_DF = _summarize_length(class_DF, sampleID, exp_factor, exp_factor_val)

        ##Length arrays for violin plot
        length_Dct[sampleID] = np.array(class_DF['length'])

        err_DF = _summarize_errors(class_DF, args.cfg()['intrapriming_perc_A_cutoff'], sampleID, exp_factor, exp_factor_val)

        ##Calculate junction cv for each of the samples
        cv_DF = calc_jxn_cv(jxn_DF, class_DF, ref_DF, dropFlag=True)
        cv_DF['sampleID'] = sampleID
        cv_DF[exp_factor] = exp_factor_val

    
        ##Store all summary dataframes in dictionaries
        gene_count_dfs[sampleID] = gene_count_DF
        ujc_count_dfs[sampleID] = ujc_count_DF
        length_dfs[sampleID] = length_summary_DF
        cv_dfs[sampleID] = cv_DF
        err_dfs[sampleID] = err_DF
        fsm_dfs[sampleID] = FSM_DF
        ism_dfs[sampleID] = ISM_DF
        nic_nnc_dfs[sampleID] = NIC_NNC_DF
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
    ism_dfs, nic_nnc_dfs, nov_can_dfs = _reord(ism_dfs), _reord(nic_nnc_dfs), _reord(nov_can_dfs)
    length_Dct = _reord(length_Dct)

    return(ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct )

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


def prep_tables(args, ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct ):
    
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
    
    #Cat subcategory DFs
    
    FSM_DF = pd.concat(fsm_dfs.values(), sort=False)
    FSM_DF.fillna(0, inplace=True)
    
    ISM_DF = pd.concat(ism_dfs.values(), sort=False)
    ISM_DF.fillna(0, inplace=True)
    
    NIC_NNC_DF = pd.concat(nic_nnc_dfs.values(), sort=False)
    NIC_NNC_DF.fillna(0, inplace=True)
    
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
        ],
        alltables_tables=[
            (err_DF, '_err_counts.csv'),
            (FSM_DF, '_FSM_counts.csv'),
            (ISM_DF, '_ISM_counts.csv'),
            (NIC_NNC_DF, '_NIC_NNC_counts.csv'),
            (nov_can_DF, '_jxn_counts.csv'),
        ],
    )

    return (gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct)




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
        color_mapping = {
            'annotated_with_well_covered_FSM': 'lightblue',
            'annotated_with_low_coverage_FSM': 'purple',
            'underannotated_with_candidate_transcript': 'darkorange',
            'underannotated_no_candidate_transcripts': 'lightsalmon'
        }

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

def prep_data_4_plots(args, gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct):
    
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
        
    
    
    cv_don_summary = cv_don_DF.groupby(['sampleID',exp_factor]).apply(
            lambda x: pd.Series({
                    'ref_match': (x['mean_abs_diff'] == 0).sum(),
                    'cv_0': ((x['cv'] == 0) & (x['mean_abs_diff'] != 0)).sum(),
                    'cv_gt_0': (x['cv'] > 0).sum()
                    }),
                    include_groups=False,
                    ).reset_index()
    cv_don_summary.fillna(0, inplace=True)
    
    cv_acc_summary = cv_acc_DF.groupby(['sampleID',exp_factor]).apply(
            lambda x: pd.Series({
                    'ref_match': (x['mean_abs_diff'] == 0).sum(),
                    'cv_0': ((x['cv'] == 0) & (x['mean_abs_diff'] != 0)).sum(),
                    'cv_gt_0': (x['cv'] > 0).sum()
                    }),
                    include_groups=False,
                    ).reset_index()
    cv_acc_summary.fillna(0, inplace=True)
   
    
    cv_acc_percs = cv_acc_summary.copy()
    cv_acc_totals = cv_acc_summary[['ref_match', 'cv_0', 'cv_gt_0']].sum(axis=1)
    cv_acc_percs[['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0']] = cv_acc_summary[['ref_match', 'cv_0', 'cv_gt_0']].div(cv_acc_totals, axis=0)*100
    cv_acc_percs.drop(columns=['ref_match', 'cv_0', 'cv_gt_0'], inplace=True)
    
 
    cv_don_percs = cv_don_summary.copy()
    cv_don_totals = cv_don_summary[['ref_match', 'cv_0', 'cv_gt_0']].sum(axis=1)
    cv_don_percs[['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0']] = cv_don_summary[['ref_match', 'cv_0', 'cv_gt_0']].div(cv_don_totals, axis=0)*100
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

    ##NIC_NNC
    NIC_NNC_DF = NIC_NNC_DF.copy()
    
    NIC_NNC_perc_DF = NIC_NNC_DF.copy()

    # Calculate the sum for the subcategories (excluding 'sampleID' and exp_factor)
    categories = [col for col in NIC_NNC_perc_DF.columns if col not in ['sampleID', exp_factor]]
    NIC_NNC_perc_DF['total'] = NIC_NNC_perc_DF[categories].sum(axis=1)
    
    # Calculate the percentage for each subcategory
    for category in categories:
        NIC_NNC_perc_DF[category] = (NIC_NNC_perc_DF[category] / NIC_NNC_perc_DF['total']) * 100
    NIC_NNC_perc_DF.drop('total', axis=1, inplace=True)
    NIC_NNC_perc_DF['sampleID'] = NIC_NNC_DF['sampleID']
    NIC_NNC_perc_DF[exp_factor] = NIC_NNC_DF[exp_factor]
    
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
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_NNC_perc_DF, nov_can_DF, nov_can_perc_DF, 
             length_DF2, cv_acc_percs, cv_don_percs)
    
    
    
def render_report_pdf(out_path, all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF, annot_gene_percs_pivot_DF, gene_agg_DF,
             gene_percs_unstacked, melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF,
             length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio,
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_NNC_perc_DF,nov_can_DF, nov_can_perc_DF,
             length_DF2,cv_acc_percs, cv_don_percs, pdf=None, ujc_metrics=None, args=None):
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
        
        
        _render_stacked_bar(pdf, NIC_NNC_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Number of Reads in Each Sub-Category - NIC and NNC ',
                            xlabel='Sample ID', ylabel='Number of reads',
                            legend_title='Structural Category')
        
         
        _render_stacked_bar(pdf, NIC_NNC_perc_DF, exp_factor=exp_factor,
                            num_factors=num_factors, palette=subcat_color_palette,
                            title='Percentage of NIC/NNC Reads in Each Sub-Category ',
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
                            title='Number of Detected Acceptors',
                            xlabel='Sample ID', ylabel='Number of Detected Acceptors',)
        
        _render_stacked_bar(pdf, cv_don_summary, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['ref_match', 'cv_0', 'cv_gt_0'],
                            title='Number of Detected Donors',
                            xlabel='Sample ID', ylabel='Number of Detected Donors',)
        
        
        categories = ['perc_ref_match','perc_cv_0','perc_cv_gt_0']
        _render_stacked_bar(pdf, cv_don_percs, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0'],
                            title='Percentage of Detected Donors',
                            xlabel='Sample ID', ylabel='Percentage of Detected Donors',)
        
        categories = ['perc_ref_match','perc_cv_0','perc_cv_gt_0']
        _render_stacked_bar(pdf, cv_acc_percs, exp_factor=exp_factor,
                            num_factors=num_factors, palette=list(three_series_palette),
                            categories=['perc_ref_match', 'perc_cv_0', 'perc_cv_gt_0'],
                            title='Percentage of Detected Acceptors',
                            xlabel='Sample ID', ylabel='Percentage of Detected Acceptors',)

        # UJC saturation, replicate concordance and UpSet plots
        plot_ujc_metrics_pages(pdf, ujc_metrics, factor_col=exp_factor)


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
    ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct = proc_samples(args, args.inDESIGN, args.inREF)

    gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF,FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct = prep_tables(args, ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs,
                                                                                                                            fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct)
    dfs_for_plotting = prep_data_4_plots( args, gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct )

    # UJC-level metrics (saturation, replicate concordance, UpSet) shared by both
    # renderers. Computed before identify_cand_underannot mutates ujc_count_DF.
    ujc_metrics = compute_ujc_metrics(ujc_count_DF, factor_col=args.inFACTOR)

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
        build_html_report(report_html, dfs_for_plotting, args, ujc_metrics=ujc_metrics)

    # Single unified PDF report: QC plots followed by the under-annotation
    # section, all in one PDF. identify_cand_underannot also (re)writes its CSVs;
    # passing the open `pdf` makes it append pages here.
    if need_pdf:
        with PdfPages(report_pdf) as pdf:
            render_report_pdf(report_pdf, *dfs_for_plotting, pdf=pdf, ujc_metrics=ujc_metrics, args=args)
            identify_cand_underannot(args, ujc_count_DF, factor_level=args.FACTORLVL, pdf=pdf)

    # Close all remaining figures to free memory
    plt.close('all')
