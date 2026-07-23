"""Interactive HTML report for SQANTI-reads.

Builds a single, self-contained ``{prefix}_report.html`` from the same prepared
DataFrames that feed the PDF report (``dfs_for_plotting``), using Plotly. The
page embeds plotly.js inline so it works fully offline (no poppler, no network).

It also writes a machine-readable ``{prefix}_qc_summary.json`` with per-sample
metrics and pass/warn/fail flags for pipeline integration.

Entry point: :func:`build_html_report`.
"""
import json
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import plotly.colors as pcolors
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs

from src.module_logging import reads_logger
from src.utilities.sqanti_reads_plots import (
    category_color_palette,
    subcat_color_palette,
    nnc_subcat_color_palette,
    cat_order,
    jxn_palette,
    three_series_palette,
    readcount_palette,
    length_palette,
    sample_seq,
    underannot_palette,
    compute_upset_intersections,
    missing_support_types,
    _ism_fragment_pct,
    _novel_noncanonical_pct,
)
from src.utilities.sqanti_reads_config import DEFAULT_CONFIG

# Fallback qualitative colors for series that have no entry in an explicit palette.
# Prefix with the shared per-sample sequence so sample colors match the PDF.
_DEFAULT_SEQ = list(sample_seq) + pcolors.qualitative.Light24
# Preferred palette for exactly-three-series plots (green, orange, yellow) — shared.
_THREE_SEQ = three_series_palette

# --- QC flag thresholds --------------------------------------------------------
# Single source of truth is DEFAULT_CONFIG["qc_flags"] in sqanti_reads_config.py;
# used here as the fallback when a run supplies no --config. "warn"/"fail" are
# lower bounds for metrics where higher is worse, upper bounds where lower is worse.
QC_THRESHOLDS = DEFAULT_CONFIG["qc_flags"]

# Legend display names for the donor/acceptor position-consistency categories
# (data columns keep their raw names; only the legend/hover labels are prettified).
_CV_COUNT_LABELS = {"ref_match": "ref_match", "cv_0": "CV=0", "cv_gt_0": "CV>0"}
_CV_PERC_LABELS = {"perc_ref_match": "ref_match", "perc_cv_0": "CV=0", "perc_cv_gt_0": "CV>0"}

_FLAG_COLORS = {"pass": "#4CAF50", "warn": "#FF9800", "fail": "#F44336"}
_PLOT_BG = "#ffffff"


def _flag_for(value, spec):
    """Return 'pass'|'warn'|'fail' for a value given a threshold spec."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "pass"
    if spec["worse"] == "high":
        if value >= spec["fail"]:
            return "fail"
        if value >= spec["warn"]:
            return "warn"
        return "pass"
    else:  # lower is worse
        if value <= spec["fail"]:
            return "fail"
        if value <= spec["warn"]:
            return "warn"
        return "pass"


def _factor_col(args):
    return args.inFACTOR if getattr(args, "inFACTOR", None) else "temp_factor"


def _value_cols(df, exclude):
    return [c for c in df.columns if c not in exclude]


def _base_layout(fig, title, xtitle, ytitle):
    fig.update_layout(
        title=dict(text=title, x=0.01, xanchor="left", font=dict(size=18)),
        xaxis_title=xtitle,
        yaxis_title=ytitle,
        barmode="stack",
        template="plotly_white",
        plot_bgcolor=_PLOT_BG,
        paper_bgcolor=_PLOT_BG,
        margin=dict(l=60, r=30, t=60, b=80),
        height=460,
        legend=dict(orientation="v", yanchor="top", y=1, xanchor="left", x=1.02),
    )
    fig.update_xaxes(tickangle=-45)
    return fig


def _stacked_bar(df, x_col, value_cols, title, xtitle, ytitle, color_map=None,
                 order=None, barmode="stack", facet_col=None, labels=None):
    """Generic stacked/grouped bar: one trace per value column.

    When ``facet_col`` is a real design factor with >1 level, the plot is split
    into one subplot column per factor level (mirroring the PDF's FacetGrid),
    with a single shared legend. ``labels`` optionally maps a column name to the
    display name shown in the legend/hover (data columns are unchanged).
    """
    cols = [c for c in (order or value_cols) if c in value_cols]
    labels = labels or {}

    def _disp(c):
        return labels.get(c, str(c))
    # Stable color per series so a category keeps the same color in every facet
    # and appears exactly once in the legend (auto-colors would drift per subplot).
    # Three-series plots without an explicit palette use green/orange/yellow.
    seq = _THREE_SEQ if (color_map is None and len(cols) == 3) else _DEFAULT_SEQ
    col_color = {}
    for i, c in enumerate(cols):
        col_color[c] = (color_map.get(c) if color_map else None) or seq[i % len(seq)]

    def _add(fig, d, row=None, col=None, show_legend=True):
        x = d[x_col].astype(str).tolist()
        for c in cols:
            fig.add_trace(go.Bar(name=_disp(c), x=x, y=d[c].tolist(),
                                 marker_color=col_color[c], legendgroup=_disp(c),
                                 showlegend=show_legend,
                                 hovertemplate=f"%{{x}}<br>{_disp(c)}: %{{y:.2f}}<extra></extra>"),
                          row=row, col=col)

    faceted = (facet_col and facet_col in df.columns
               and facet_col != "temp_factor" and df[facet_col].nunique() > 1)
    if faceted:
        levels = list(pd.unique(df[facet_col]))
        fig = make_subplots(rows=1, cols=len(levels), shared_yaxes=True,
                            horizontal_spacing=0.08,
                            # short level labels (the report header names the factor)
                            subplot_titles=[str(lv) for lv in levels])
        for i, lv in enumerate(levels, start=1):
            _add(fig, df[df[facet_col] == lv], row=1, col=i, show_legend=(i == 1))
        # No internal chart title on faceted plots: the section heading already
        # names the figure, and a centered title collides with the subplot titles.
        _base_layout(fig, "", xtitle, ytitle)
        # Legend below the panels so the facet columns (and their titles) get the
        # full width and don't overlap at narrow viewports.
        fig.update_layout(barmode=barmode, margin=dict(l=60, r=20, t=50, b=140),
                          legend=dict(orientation="h", yanchor="top", y=-0.28,
                                      xanchor="left", x=0))
        fig.update_annotations(font_size=13)  # shrink subplot (facet) titles
        fig.update_xaxes(tickangle=-45)
        fig.update_yaxes(title_text=ytitle, row=1, col=1)
        return fig

    fig = go.Figure()
    _add(fig, df, show_legend=True)
    _base_layout(fig, title, xtitle, ytitle)
    fig.update_layout(barmode=barmode)
    return fig


# A violin/box KDE is computed client-side from the raw points Plotly embeds, so
# a per-read distribution over millions of reads inflates the self-contained HTML
# to hundreds of MB. A few thousand points per group give a visually identical
# density/box, so each group is randomly subsampled to this cap before plotting.
_VIOLIN_MAX_POINTS = 4000


def _cap_per_group(df, group_cols, cap=_VIOLIN_MAX_POINTS):
    """Randomly subsample each group to at most ``cap`` rows (reproducible), so a
    distribution plot keeps its shape without embedding every raw data point.

    Shuffles once (deterministically) then keeps the first ``cap`` rows of each
    group — equivalent to a per-group random sample, without a groupby.apply."""
    if df is None or df.empty:
        return df
    cols = [c for c in group_cols if c in df.columns]
    if not cols:
        return df.sample(min(cap, len(df)), random_state=0)
    shuffled = df.sample(frac=1, random_state=0)
    keep = shuffled.groupby(cols, observed=True, sort=False).cumcount() < cap
    return shuffled[keep]


def _violin(df, x_col, y_col, title, xtitle, ytitle, facet_col=None):
    """One violin (distribution) per x category; colored by factor when set."""
    fig = go.Figure()
    faceted = (facet_col and facet_col in df.columns
               and facet_col != "temp_factor" and df[facet_col].nunique() > 1)
    # Cap points per violin (per x-group, and per facet level when faceted).
    df = _cap_per_group(df, [x_col] + ([facet_col] if faceted else []))
    if faceted:
        for i, lv in enumerate(pd.unique(df[facet_col])):
            d = df[df[facet_col] == lv]
            fig.add_trace(go.Violin(x=d[x_col].astype(str), y=d[y_col], name=str(lv),
                                    legendgroup=str(lv), scalegroup=str(lv),
                                    box_visible=True, meanline_visible=True,
                                    marker_color=_DEFAULT_SEQ[i % len(_DEFAULT_SEQ)]))
        fig.update_layout(showlegend=True)
    else:
        fig.add_trace(go.Violin(x=df[x_col].astype(str), y=df[y_col],
                                box_visible=True, meanline_visible=True,
                                marker_color="#1f6feb", showlegend=False))
    _base_layout(fig, title, xtitle, ytitle)
    fig.update_layout(barmode="overlay")
    return fig


def _scatter_xy(df, x_col, y_col, label_col, title, xtitle, ytitle, facet_col=None):
    """Scatter with per-point labels; colored by factor when set."""
    fig = go.Figure()
    faceted = (facet_col and facet_col in df.columns
               and facet_col != "temp_factor" and df[facet_col].nunique() > 1)
    groups = pd.unique(df[facet_col]) if faceted else [None]
    for i, lv in enumerate(groups):
        d = df if lv is None else df[df[facet_col] == lv]
        fig.add_trace(go.Scatter(
            x=d[x_col], y=d[y_col], mode="markers+text", text=d[label_col],
            textposition="top center", name=(str(lv) if lv is not None else ""),
            marker=dict(size=13, color=(_DEFAULT_SEQ[i % len(_DEFAULT_SEQ)] if lv is not None else "#1f6feb")),
            hovertemplate=f"%{{text}}<br>{xtitle}: %{{x:.1f}}<br>{ytitle}: %{{y:.0f}}<extra></extra>"))
    _base_layout(fig, title, xtitle, ytitle)
    fig.update_layout(barmode="group", showlegend=faceted)
    return fig


def _heatmap(matrix_df, title, xtitle, ytitle, colorscale="RdBu", zmid=0,
             zmin=None, zmax=None, colorbar_title="value", value_fmt=".2f",
             annotate=False):
    """Heatmap from a DataFrame (rows=y labels, columns=x labels)."""
    fig = go.Figure(go.Heatmap(
        z=matrix_df.values, x=[str(c) for c in matrix_df.columns],
        y=[str(i) for i in matrix_df.index], colorscale=colorscale, zmid=zmid,
        zmin=zmin, zmax=zmax, colorbar=dict(title=colorbar_title),
        hovertemplate="%{y}<br>%{x}: %{z:" + value_fmt + "}<extra></extra>"))
    if annotate:
        import numpy as _np
        vals = matrix_df.values.astype(float)
        # Effective color range, to place each cell on the 0–1 colorscale.
        lo = zmin if zmin is not None else _np.nanmin(vals)
        hi = zmax if zmax is not None else _np.nanmax(vals)
        span = (hi - lo) or 1.0
        # Sequential dark-at-low scales (Viridis etc.) are dark in their lower half,
        # where black text is unreadable — use white there.
        dark_low = str(colorscale).lower() in ("viridis", "cividis", "plasma", "magma", "inferno")
        for yi, idx in enumerate(matrix_df.index):
            for xi, col in enumerate(matrix_df.columns):
                v = vals[yi, xi]
                norm = (v - lo) / span
                text_color = "#ffffff" if (dark_low and norm < 0.45) else "#222"
                fig.add_annotation(x=str(col), y=str(idx),
                                   text=format(v, value_fmt),
                                   showarrow=False, font=dict(size=11, color=text_color))
    _base_layout(fig, title, xtitle, ytitle)
    fig.update_layout(barmode="group",
                      height=max(460, 26 * len(matrix_df.index) + 160))
    return fig


def _upset_figure(upset_data, samples):
    """Plotly UpSet mirroring the PDF: intersection bars (top-right, stacked by
    structural category), per-sample set-size bars (bottom-left, also stacked by
    SC), and a membership-dot matrix (bottom-right)."""
    inter = upset_data["intersections"]
    set_sizes = upset_data["set_sizes"]
    sc_order = upset_data["sc_order"]
    n = len(inter)
    ns = len(samples)
    x = list(range(n))
    fig = make_subplots(rows=2, cols=2, row_heights=[0.7, 0.3],
                        column_widths=[0.22, 0.78],
                        shared_xaxes=True, shared_yaxes=True,
                        horizontal_spacing=0.015, vertical_spacing=0.05,
                        specs=[[None, {}], [{}, {}]])
    for c in sc_order:
        color = category_color_palette.get(c, "#969696")
        # intersection bars (top-right)
        fig.add_trace(go.Bar(x=x, y=[d["sc_counts"].get(c, 0) for d in inter],
                             name=str(c), marker_color=color, legendgroup=str(c),
                             hovertemplate=f"{c}: %{{y}}<extra></extra>"), row=1, col=2)
        # set-size bars (bottom-left), horizontal, stacked by SC
        fig.add_trace(go.Bar(y=list(range(ns)),
                             x=[set_sizes.get(s, {}).get(c, 0) for s in samples],
                             orientation="h", marker_color=color, legendgroup=str(c),
                             showlegend=False,
                             hovertemplate=f"{c}: %{{x}}<extra></extra>"), row=2, col=1)
    # membership dot matrix (bottom-right)
    idx = {s: i for i, s in enumerate(samples)}
    grey_x, grey_y, dark_x, dark_y = [], [], [], []
    for j, d in enumerate(inter):
        members = set(d["combo"])
        for s in samples:
            (dark_x if s in members else grey_x).append(j)
            (dark_y if s in members else grey_y).append(idx[s])
        present = sorted(idx[s] for s in members)
        if len(present) > 1:
            fig.add_trace(go.Scatter(x=[j, j], y=[present[0], present[-1]],
                                     mode="lines", line=dict(color="#333", width=3),
                                     showlegend=False, hoverinfo="skip"), row=2, col=2)
    fig.add_trace(go.Scatter(x=grey_x, y=grey_y, mode="markers",
                             marker=dict(color="#dddddd", size=14), showlegend=False,
                             hoverinfo="skip"), row=2, col=2)
    fig.add_trace(go.Scatter(x=dark_x, y=dark_y, mode="markers",
                             marker=dict(color="#333333", size=14), showlegend=False,
                             hoverinfo="skip"), row=2, col=2)
    fig.update_yaxes(title_text="# UJCs", row=1, col=2)
    # samples on the shared row-2 y-axis; labels on the left (set-size) panel
    fig.update_yaxes(tickmode="array", tickvals=list(range(ns)), ticktext=samples,
                     autorange="reversed", row=2, col=1)
    fig.update_xaxes(title_text="Set size", autorange="reversed", row=2, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=2)
    fig.update_xaxes(showticklabels=False, row=2, col=2)
    fig.update_layout(barmode="stack", template="plotly_white", height=560,
                      plot_bgcolor=_PLOT_BG, paper_bgcolor=_PLOT_BG,
                      margin=dict(l=90, r=30, t=40, b=40),
                      legend=dict(title="Structural category", orientation="v",
                                  x=1.02, y=1))
    return fig


def _section(title, fig, interpretation, div_id):
    """Render one report section (heading + interpretation + Plotly div)."""
    plot_div = pio.to_html(fig, full_html=False, include_plotlyjs=False,
                           div_id=div_id, config={"displaylogo": False,
                                                   "responsive": True})
    return f"""
    <section class="card">
      <h2>{title}</h2>
      <p class="interp">{interpretation}</p>
      {plot_div}
    </section>
    """


def _compute_summary(length_DF, err_DF, all_gene_percs_pivot_DF, gene_agg_DF, args,
                     thresholds=QC_THRESHOLDS, jxn_offset_metrics=None,
                     completeness_metrics=None):
    """Per-sample metrics + flags. Returns (summary_dict, samples_in_order)."""
    exp_factor = _factor_col(args)
    samples = length_DF["sampleID"].astype(str).tolist()

    length_idx = length_DF.set_index("sampleID")
    err_idx = err_DF.set_index("sampleID")
    fsm = all_gene_percs_pivot_DF.set_index("sampleID")["FSM"] \
        if "FSM" in all_gene_percs_pivot_DF.columns else None
    gene_cols = _value_cols(gene_agg_DF, {"sampleID", exp_factor})
    genes_detected = gene_agg_DF.set_index("sampleID")[gene_cols].sum(axis=1)

    # Per-sample splice-site imprecision scalar (share of near-reference site
    # observations off by >0 bp), keyed by sampleID; None-safe when unavailable.
    imprecise = {}
    if jxn_offset_metrics and not jxn_offset_metrics.get("per_sample", None) is None:
        ps = jxn_offset_metrics["per_sample"]
        if not ps.empty:
            imprecise = ps.set_index("sampleID")["perc_imprecise"].astype(float).to_dict()

    # Per-sample 5'/3' completeness scalars (% of reads within the window),
    # keyed by sampleID; None-safe when completeness could not be computed.
    compl_5p, compl_3p = {}, {}
    if completeness_metrics and not completeness_metrics.get("per_sample", None) is None:
        cp = completeness_metrics["per_sample"]
        if not cp.empty:
            cpi = cp.set_index("sampleID")
            compl_5p = cpi["perc_5p_within_window"].astype(float).to_dict()
            compl_3p = cpi["perc_3p_within_window"].astype(float).to_dict()

    per_sample = {}
    for s in samples:
        metrics = {
            "total_reads": int(length_idx.at[s, "total_reads"]),
            "genes_detected": int(genes_detected.get(s, 0)),
            "median_length": float(length_idx.at[s, "median_length"]),
            "average_length": float(length_idx.at[s, "average_length"]),
            "perc_reads_gt_1kb": float(length_idx.at[s, "perc_reads_gt_1kb"]),
            "perc_FSM": float(fsm.get(s)) if fsm is not None else None,
            "perc_reads_intrapriming": float(err_idx.at[s, "perc_reads_intrapriming"]),
            "perc_reads_RTS": float(err_idx.at[s, "perc_reads_RTS"]),
            "perc_reads_non-canonical": float(err_idx.at[s, "perc_reads_non-canonical"]),
            "perc_sites_imprecise": float(imprecise[s]) if s in imprecise else None,
            "perc_5p_within_window": float(compl_5p[s]) if s in compl_5p else None,
            "perc_3p_within_window": float(compl_3p[s]) if s in compl_3p else None,
        }
        flags = {m: _flag_for(metrics.get(m), spec) for m, spec in thresholds.items()}
        metrics["flags"] = flags
        metrics["overall_flag"] = ("fail" if "fail" in flags.values()
                                    else "warn" if "warn" in flags.values() else "pass")
        per_sample[s] = metrics
    return per_sample, samples


def _summary_html(per_sample, samples, thresholds=QC_THRESHOLDS):
    """Summary metric table + QC flag panel."""
    metric_labels = {spec["label"]: m for m, spec in thresholds.items()}
    # Metrics table
    head_cols = ["Sample", "Reads", "Genes", "Median len", "% &gt;1kb", "% FSM",
                 "% intra-prim", "% RTS", "% non-canon", "% 5' compl", "% 3' compl",
                 "% fuzzy sites", "Overall"]
    rows = []
    for s in samples:
        m = per_sample[s]
        badge = f'<span class="badge {m["overall_flag"]}">{m["overall_flag"].upper()}</span>'
        # (value, flaggable-metric-name or None) per column.
        cells = [
            (s, None),
            (f'{m["total_reads"]:,}', None),
            (m["genes_detected"], None),
            (f'{m["median_length"]:.0f}', "median_length"),
            (f'{m["perc_reads_gt_1kb"]:.1f}', None),
            ((f'{m["perc_FSM"]:.1f}' if m["perc_FSM"] is not None else "—"), None),
            (f'{m["perc_reads_intrapriming"]:.1f}', "perc_reads_intrapriming"),
            (f'{m["perc_reads_RTS"]:.1f}', "perc_reads_RTS"),
            (f'{m["perc_reads_non-canonical"]:.1f}', "perc_reads_non-canonical"),
            ((f'{m["perc_5p_within_window"]:.1f}' if m.get("perc_5p_within_window") is not None else "—"),
             "perc_5p_within_window"),
            ((f'{m["perc_3p_within_window"]:.1f}' if m.get("perc_3p_within_window") is not None else "—"),
             "perc_3p_within_window"),
            ((f'{m["perc_sites_imprecise"]:.1f}' if m["perc_sites_imprecise"] is not None else "—"),
             "perc_sites_imprecise"),
        ]
        tds = []
        for val, metric in cells:
            fl = m["flags"].get(metric) if metric else None
            if fl in ("warn", "fail"):
                # Bold + flag-tinted the cells that trigger a warn/fail.
                tds.append(f'<td class="flag-{fl}"><strong>{val}</strong></td>')
            else:
                tds.append(f"<td>{val}</td>")
        rows.append("<tr>" + "".join(tds) + f"<td>{badge}</td></tr>")
    table = ("<table class='summary'><thead><tr>"
             + "".join(f"<th>{c}</th>" for c in head_cols)
             + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>")

    # Flag legend
    legend = ("<div class='flag-legend'>"
              + "".join(f"<span class='badge {k}'>{k.upper()}</span>" for k in
                        ["pass", "warn", "fail"])
              + " &nbsp; thresholds: "
              + "; ".join(f"{spec['label']} warn≥{spec['warn']}/fail≥{spec['fail']}"
                          if spec["worse"] == "high" else
                          f"{spec['label']} warn≤{spec['warn']}/fail≤{spec['fail']}"
                          for spec in thresholds.values())
              + "</div>")
    return f"<section class='card'><h2>Summary &amp; QC flags</h2>{table}{legend}</section>"


def _gene_classification_section(out_dir, prefix):
    """Optional under-annotation gene-category bar, read from the CSV on disk."""
    path = os.path.join(out_dir, prefix + "_gene_classfication.csv")
    if not os.path.isfile(path):
        return ""
    df = pd.read_csv(path)
    if df.empty or "gene_category" not in df.columns:
        return ("<section class='card'><h2>Under-annotation analysis</h2>"
                "<p class='interp'>No genes met the under-annotation criteria.</p></section>")
    counts = df["gene_category"].value_counts()
    fig = go.Figure(go.Bar(
        x=[c.replace("_", " ") for c in counts.index], y=counts.values,
        marker_color=[underannot_palette.get(c, "#969696") for c in counts.index],
        hovertemplate="%{x}<br>%{y} genes<extra></extra>"))
    _base_layout(fig, "Gene annotation categories", "Category", "Number of genes")
    fig.update_layout(showlegend=False, barmode="group")
    interp = (
        "SQANTI-reads looks for genes that may be <b>under-annotated</b> — where the reads "
        "suggest a transcript the reference is missing. Only genes with enough read support "
        "(the <code>-ge</code> gene-expression cut-off) are considered, and each such gene is "
        "placed in one of four categories by comparing its reads against the annotation: "
        "<b>annotated_with_well_covered_FSM</b> — a well-supported full-splice-match to a known "
        "transcript was seen, so the gene looks correctly annotated; "
        "<b>annotated_with_low_coverage_FSM</b> — a known transcript matched but only weakly; "
        "<b>underannotated_with_candidate_transcript</b> — no well-covered known transcript, but "
        "a well-supported <i>novel</i> unique junction chain is present (a strong candidate for "
        "a real transcript missing from the annotation, worth manual inspection); "
        "<b>underannotated_no_candidate_transcripts</b> — supported reads but no clear "
        "well-covered candidate. The bars count how many genes fall in each category. It is a "
        "discovery aid, not a pass/fail: the 'candidate' class is where to look for genuinely "
        "novel transcripts, but individual calls should be confirmed (e.g. in a genome browser).")
    return _section("Under-annotation analysis", fig, interp, "fig-underannot")


def _optional_support_footer(cage_metrics, polya_metrics, sr_metrics):
    """Part C — a muted closing note naming the orthogonal-support types not
    supplied and the (optional) SQANTI3 flags that would add them. Returns "" when
    every support type is present. Deliberately non-judgemental: the evidence is
    optional and the report is complete without it."""
    missing = missing_support_types(cage_metrics, polya_metrics, sr_metrics)
    if not missing:
        return ""
    items = "".join(
        f"<li>{label} — add with <code>{flags}</code></li>" for label, flags in missing)
    return f"""
    <section class="card" style="opacity:0.75;border-style:dashed;">
      <h2 style="font-size:1.05em;">Optional orthogonal support</h2>
      <p class="interp">This report is complete as it stands. The independent
      end/junction evidence below was not supplied to SQANTI3 and is entirely
      <b>optional</b> — providing it would add extra cross-checks of transcript
      5'/3' ends and splice junctions in a future run:</p>
      <ul class="interp">{items}</ul>
      <p class="interp">None of these are <b>required</b> to interpret the QC
      results above.</p>
    </section>
    """


PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script>{plotlyjs}</script>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
          margin: 0; background: #f4f6f8; color: #1a2027; }}
  header {{ background: #15918A; color: #fff; padding: 20px 28px; }}
  header h1 {{ margin: 0 0 4px; font-size: 22px; }}
  header .meta {{ opacity: .85; font-size: 13px; }}
  main {{ max-width: 1100px; margin: 0 auto; padding: 20px; }}
  .card {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 10px;
           padding: 18px 20px; margin: 0 0 20px; box-shadow: 0 1px 2px rgba(0,0,0,.04); }}
  .card h2 {{ margin: 0 0 6px; font-size: 17px; color: #1f3a5f; }}
  .interp {{ margin: 0 0 12px; color: #52606d; font-size: 13.5px; line-height: 1.5; }}
  table.summary {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  table.summary th, table.summary td {{ border-bottom: 1px solid #edf0f2;
           padding: 7px 9px; text-align: right; }}
  table.summary th:first-child, table.summary td:first-child {{ text-align: left; }}
  table.summary thead th {{ background: #f7f9fb; color: #3a4653; }}
  table.summary td.flag-warn {{ background: #FFF3E0; color: #B45309; }}
  table.summary td.flag-fail {{ background: #FFEBEE; color: #B71C1C; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 10px;
            color: #fff; font-size: 11px; font-weight: 600; }}
  .badge.pass {{ background: {pass_c}; }} .badge.warn {{ background: {warn_c}; }}
  .badge.fail {{ background: {fail_c}; }}
  .flag-legend {{ margin-top: 12px; font-size: 12px; color: #6b7683; }}
  .toc {{ font-size: 13px; }} .toc a {{ color: #1f6feb; text-decoration: none; margin-right: 14px; }}
  footer {{ text-align: center; color: #8a97a5; font-size: 12px; padding: 20px; }}
</style></head><body>
<header>
  <h1>SQANTI-reads QC report</h1>
  <div class="meta">{n_samples} samples &middot; prefix <code>{prefix}</code>{factor_note}</div>
</header>
<main>
{body}
</main>
<footer>Generated by SQANTI-reads &middot; interactive Plotly report</footer>
</body></html>"""


def _jxn_offset_figures(m, samples):
    """Return (spectrum_fig, profile_fig, byclass_fig) Plotly figures for the
    splice-site fuzziness section, or (None, None, None) if no data."""
    if not m or not m.get("samples"):
        return None, None, None
    palette = {s: _DEFAULT_SEQ[i % len(_DEFAULT_SEQ)] for i, s in enumerate(samples)}
    window = m["window"]

    # 1. Signed-offset spectrum, donor + acceptor subplots (exact excluded).
    spec = m["spectrum"]
    sfig = make_subplots(rows=1, cols=2, subplot_titles=("Donor", "Acceptor"),
                         shared_yaxes=True)
    for col, site in ((1, "donor"), (2, "acceptor")):
        ss = spec[(spec["site"] == site) & (spec["offset"] != 0)]
        for s in samples:
            d = ss[ss["sampleID"] == s].sort_values("offset")
            sfig.add_trace(go.Scatter(x=d["offset"], y=d["count"], mode="lines",
                                      name=str(s), legendgroup=str(s),
                                      showlegend=(col == 1),
                                      line=dict(color=palette[s]),
                                      hovertemplate="offset %{x} bp: %{y}<extra></extra>"),
                           row=1, col=col)
        sfig.update_xaxes(title_text=f"Offset from reference {site} (bp)",
                          range=[-window, window], row=1, col=col)
    sfig.update_yaxes(title_text="Site observations (exact excluded)", row=1, col=1)
    _base_layout(sfig, "Splice-site offset spectrum", "", "")

    # 2. Precision profile.
    prof = m["profile"]
    pfig = go.Figure()
    for s in samples:
        d = prof[prof["sampleID"] == s].sort_values("k")
        pfig.add_trace(go.Scatter(x=d["k"], y=d["perc_within"], mode="lines+markers",
                                  name=str(s), line=dict(color=palette[s]),
                                  hovertemplate="within ±%{x} bp: %{y:.1f}%<extra></extra>"))
    _base_layout(pfig, "Splice-site precision profile",
                 "Window half-width k (bp)", "% within ±k of reference")

    # 3. Canonical vs non-canonical imprecision.
    bc = m["by_class"]
    cfig = None
    if bc is not None and not bc.empty:
        classes = [c for c in ["canonical", "non_canonical"] if c in set(bc["canonical"])]
        classes += [c for c in sorted(set(bc["canonical"])) if c not in classes]
        class_colors = {"canonical": "#4CAF50", "non_canonical": "#C51B7D"}
        cfig = go.Figure()
        for can in classes:
            vals = []
            for s in samples:
                v = bc[(bc["sampleID"] == s) & (bc["canonical"] == can)]["perc_imprecise"]
                vals.append(float(v.iloc[0]) if len(v) else 0.0)
            cfig.add_trace(go.Bar(x=samples, y=vals, name=can.replace("_", " "),
                                  marker_color=class_colors.get(can, "#969696"),
                                  hovertemplate="%{y:.1f}% imprecise<extra></extra>"))
        cfig.update_layout(barmode="group")
        _base_layout(cfig, "Splice-site fuzziness by canonical class",
                     "Sample", "% site observations imprecise")
    return sfig, pfig, cfig


def _completeness_figure(m, samples):
    """Return the 5'/3' read-end completeness profile Plotly figure (|distance|
    ECDFs, one subplot per end), or None if no data."""
    if not m or not m.get("samples"):
        return None
    palette = {s: _DEFAULT_SEQ[i % len(_DEFAULT_SEQ)] for i, s in enumerate(samples)}
    window = m["window"]
    prof = m["profile"]
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("5' completeness", "3' completeness"),
                        shared_yaxes=True)
    for col, end in ((1, "5prime"), (2, "3prime")):
        ed = prof[prof["end"] == end]
        for s in samples:
            d = ed[ed["sampleID"] == s].sort_values("k")
            fig.add_trace(go.Scatter(x=d["k"], y=d["perc_within"], mode="lines+markers",
                                     name=str(s), legendgroup=str(s),
                                     showlegend=(col == 1), line=dict(color=palette[s]),
                                     hovertemplate="within %{x} bp: %{y:.1f}%<extra></extra>"),
                          row=1, col=col)
        lbl = "5'" if end == "5prime" else "3'"
        fig.update_xaxes(title_text=f"|distance to gene {lbl} end| (bp)", row=1, col=col)
    fig.update_yaxes(title_text="Cumulative % of reads", row=1, col=1)
    _base_layout(fig, "Read-end completeness profile", "", "")
    return fig


def _scorecard_figure(sc):
    """Return the cohort-relative sample-outlier scorecard heatmap (signed robust
    z; red = worse than cohort), or None when disabled."""
    if not sc or not sc.get("enabled"):
        return None
    samples = sc["samples"]
    metrics = sc["metrics"]
    Z = [[sc["z"][s][mt] for mt in metrics] for s in samples]
    labels = [mt.replace("perc_", "%").replace("_", " ") for mt in metrics]
    text = [[f"{sc['z'][s][mt]:.1f}" for mt in metrics] for s in samples]
    fig = go.Figure(data=go.Heatmap(
        z=Z, x=labels, y=samples, text=text, texttemplate="%{text}",
        zmid=0, colorscale="RdBu_r", zmin=-3.5, zmax=3.5,
        colorbar=dict(title="robust z<br>(+ = worse)"),
        hovertemplate="%{y} · %{x}: z=%{z:.2f}<extra></extra>"))
    _base_layout(fig, "Sample-outlier scorecard (cohort-relative robust z)",
                 "Metric", "Sample")
    return fig


def _metric_cohort_figure(err_DF, metric, title, ytitle, scorecard=None,
                          threshold=None, color="#F58A53"):
    """Per-sample bar of one QC rate in cohort context (median + config warn/fail
    lines + scorecard flag markers), or None if the metric is unavailable.

    The scorecard heatmap answers "which samples diverge"; this per-metric figure
    shows the raw rate behind the verdict and where the cohort sits."""
    if err_DF is None or metric not in getattr(err_DF, "columns", []):
        return None
    d = err_DF[["sampleID", metric]].dropna()
    if d.empty:
        return None
    samples = d["sampleID"].astype(str).tolist()
    vals = d[metric].astype(float).tolist()

    # Per-bar color: flag color when the scorecard flagged this sample on this
    # metric, else the base color.
    flag_color = {"warn": "#FF9800", "fail": "#F44336"}
    bar_colors, line_widths = [], []
    scored = bool(scorecard and scorecard.get("enabled")
                  and metric in scorecard.get("metrics", []))
    for s in samples:
        fl = scorecard["cell_flags"].get(s, {}).get(metric, "pass") if scored else "pass"
        bar_colors.append(flag_color.get(fl, color))
        line_widths.append(2.5 if fl in flag_color else 0)

    fig = go.Figure(go.Bar(
        x=samples, y=vals, marker_color=bar_colors,
        marker_line=dict(color="#333333", width=line_widths),
        hovertemplate="%{x}: %{y:.2f}%<extra></extra>"))

    med = float(np.median(vals))
    fig.add_hline(y=med, line_dash="dash", line_color="#777777",
                  annotation_text=f"cohort median ({med:.1f})",
                  annotation_position="top left")
    if threshold:
        for lvl, dash in (("warn", "dot"), ("fail", "dashdot")):
            if lvl in threshold:
                fig.add_hline(y=float(threshold[lvl]), line_dash=dash,
                              line_color=flag_color.get(lvl, "#F44336"),
                              annotation_text=f"{lvl} ({threshold[lvl]:.0f})",
                              annotation_position="top right")
    _base_layout(fig, title, "Sample", ytitle)
    fig.update_layout(barmode="group")
    return fig


def _composition_drift_figure(drift):
    """A3 pairwise composition-drift heatmap (Jensen–Shannon distance), or None."""
    if (not drift or not drift.get("samples")
            or drift.get("matrix") is None or drift["matrix"].empty):
        return None
    return _heatmap(drift["matrix"], "Composition drift (pairwise)",
                    "Sample", "Sample", colorscale="Viridis", zmid=None, zmin=0,
                    colorbar_title="JS distance", annotate=True)


def _tandem_sites_figure(tandem):
    """F6 near-zero offset histogram (summed over samples) with ±3/±4/±6 tandem
    bars highlighted, or None."""
    if (not tandem or not tandem.get("samples")
            or tandem.get("by_offset") is None or tandem["by_offset"].empty):
        return None
    agg = (tandem["by_offset"].groupby("offset", as_index=False)["count"].sum()
           .sort_values("offset"))
    tand = {3, -3, 4, -4, 6, -6}
    colors = ["#EE446F" if o in tand else "#BBBBBB" for o in agg["offset"]]
    fig = go.Figure(go.Bar(x=agg["offset"].tolist(), y=agg["count"].tolist(),
                           marker_color=colors,
                           hovertemplate="offset %{x} bp<br>%{y} obs<extra></extra>"))
    _base_layout(fig, "Tandem splice sites (NAGNAG)",
                 "signed splice-site offset (bp)", "imprecise observations")
    fig.update_layout(barmode="group", showlegend=False)
    return fig


def _ujc_overlap_figure(overlap):
    """A6 pairwise UJC overlap-index heatmap, or None when <2 samples carry UJCs."""
    if not overlap or overlap.get("matrix") is None:
        return None
    return _heatmap(overlap["matrix"],
                    "Pairwise UJC overlap index (row ∩ col / row)",
                    "fraction also detected in this sample", "fraction of this sample's UJCs",
                    colorscale="Viridis", zmid=None, zmin=0, zmax=1,
                    colorbar_title="overlap", value_fmt=".2f", annotate=True)


def _replicate_concordance_figure(concordance):
    """A5 multi-axis replicate-concordance grouped bar, or None when disabled."""
    m = concordance
    if not m or not m.get("enabled") or m.get("per_sample") is None or m["per_sample"].empty:
        return None
    ps = m["per_sample"]
    axes = m["axes"]
    x = [f"{s} ({g})" for s, g in zip(ps["sampleID"], ps["group"])]
    colors = {"category_composition": "#15918A", "length_profile": "#F58A53",
              "imprecision": "#74CDF0"}
    # A8 companion — the group's robust within-group spread (MAD) for each axis,
    # surfaced on hover next to the agreement value.
    mad_of = {"category_composition": "mad_category_composition",
              "length_profile": "mad_length_profile", "imprecision": "mad_imprecision"}
    mad_units = {"category_composition": "JSD", "length_profile": "bp", "imprecision": "pp"}
    fig = go.Figure()
    for axis in axes:
        madcol = mad_of.get(axis)
        if madcol and madcol in ps.columns:
            cust = [("n/a" if pd.isna(v) else f"{v:.3g} {mad_units[axis]}")
                    for v in ps[madcol]]
        else:
            cust = ["n/a"] * len(ps)
        fig.add_trace(go.Bar(
            name=axis.replace("_", " "), x=x, y=ps[axis].fillna(0).tolist(),
            marker_color=colors.get(axis, "#999999"), customdata=cust,
            hovertemplate="%{x}<br>agreement: %{y:.3f}<br>group spread (MAD): "
                          "%{customdata}<extra>" + axis.replace("_", " ") + "</extra>"))
    _base_layout(fig, "Replicate concordance (multi-axis)",
                 "Sample (group)", "within-group agreement (1 = matches replicates)")
    fig.update_layout(barmode="group", yaxis_range=[0, 1.05])
    return fig


def _fuzziness_concordance_figure(concordance):
    """F7 per-sample splice-site precision agreement bar, or None when disabled."""
    m = concordance
    if not m or not m.get("enabled") or m.get("per_sample") is None or m["per_sample"].empty:
        return None
    ps = m["per_sample"]
    x = [f"{s} ({g})" for s, g in zip(ps["sampleID"], ps["group"])]
    fig = go.Figure(go.Bar(x=x, y=ps["site_precision_agreement"].fillna(0).tolist(),
                           marker_color="#9F7BB8",
                           hovertemplate="%{x}<br>agreement %{y:.2f}<extra></extra>"))
    _base_layout(fig, "Replicate concordance of splice-site precision",
                 "Sample (group)", "site-precision agreement")
    fig.update_layout(barmode="group", yaxis_range=[0, 1.05], showlegend=False)
    return fig


def _fuzz_depth_figure(metrics):
    """F4 imprecision-vs-depth line per sample, or None (fixture has no coverage)."""
    m = metrics
    if not m or not m.get("samples") or m.get("profile") is None or m["profile"].empty:
        return None
    prof = m["profile"].copy()
    prof["depth_bin"] = prof["depth_bin"].astype(str)
    fig = go.Figure()
    for i, (s, g) in enumerate(prof.groupby("sampleID")):
        fig.add_trace(go.Scatter(x=g["depth_bin"].tolist(), y=g["perc_imprecise"].tolist(),
                                 mode="lines+markers", name=str(s),
                                 line=dict(color=_DEFAULT_SEQ[i % len(_DEFAULT_SEQ)])))
    _base_layout(fig, "Splice-site imprecision vs junction depth",
                 "junction read depth (total_coverage_unique)",
                 "% imprecise splice-site observations")
    fig.update_layout(barmode="group")
    return fig


def build_html_report(out_path, dfs_for_plotting, args, ujc_metrics=None,
                      jxn_offset_metrics=None, completeness_metrics=None,
                      scorecard=None, yield_metrics=None, drift_metrics=None,
                      tandem_metrics=None, rep_concordance=None, fuzz_concordance=None,
                      fuzz_depth_metrics=None, overlap_metrics=None,
                      cage_metrics=None, polya_metrics=None, sr_metrics=None):
    """Build the interactive HTML report and the qc_summary.json sidecar.

    Parameters
    ----------
    out_path : str
        Destination ``*_report.html`` path.
    dfs_for_plotting : tuple
        The tuple returned by ``prep_data_4_plots`` (same as feeds the PDF).
    args : ReadsPlotArgs
        Plotting config (thresholds, factor, prefix, OUT).
    """
    (all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF,
     annot_gene_percs_pivot_DF, gene_agg_DF, gene_percs_unstacked,
     melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF,
     length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio,
     cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_DF, NNC_DF, FSM_perc_DF,
     ISM_perc_DF, NIC_perc_DF, NNC_perc_DF, nov_can_DF, nov_can_perc_DF, length_DF2,
     cv_acc_percs, cv_don_percs) = dfs_for_plotting

    exp_factor = _factor_col(args)
    drop = {"sampleID", exp_factor}
    facet = getattr(args, "inFACTOR", None)  # None -> single (non-faceted) plots

    # QC flag thresholds come from the (optionally user-supplied) config.
    qc_flags = args.cfg().get("qc_flags", QC_THRESHOLDS) if hasattr(args, "cfg") else QC_THRESHOLDS

    per_sample, samples = _compute_summary(length_DF, err_DF, all_gene_percs_pivot_DF,
                                           gene_agg_DF, args, thresholds=qc_flags,
                                           jxn_offset_metrics=jxn_offset_metrics,
                                           completeness_metrics=completeness_metrics)

    sections = [_summary_html(per_sample, samples, thresholds=qc_flags)]

    # 0. Cohort-relative sample-outlier scorecard, right below the summary table.
    _scfig = _scorecard_figure(scorecard)
    if _scfig is not None:
        sections.append(_section(
            "Sample-outlier scorecard", _scfig,
            "Each read-QC metric is turned into a robust z-score against the cohort "
            "(median/MAD across samples), oriented so a positive value (red) means "
            "worse than peers. A sample is flagged only when it diverges on several "
            "independent metrics at once — the score is relative, so it stays quiet "
            "when all samples agree and needs no dataset-specific thresholds.",
            "fig-scorecard"))
    elif scorecard and scorecard.get("reason"):
        sections.append(
            f"<div class='section'><h2>Sample-outlier scorecard</h2>"
            f"<p class='muted'>Not shown: {scorecard['reason']}.</p></div>")

    subcat_order = list(subcat_color_palette.keys())
    readcount_order = ["100+ reads", "50-100 reads", "11-50 reads", "2-10 reads", "1 read"]
    length_order = ["reads_lt_1kb_perc", "reads_1kb_to_2kb_perc",
                    "reads_2kb_to_3kb_perc", "reads_gt_3kb_perc"]

    # 1. Structural category % (all genes)
    sections.append(_section(
        "Reads per structural category — all genes",
        _stacked_bar(all_gene_percs_pivot_DF, "sampleID",
                     _value_cols(all_gene_percs_pivot_DF, drop),
                     "Percent reads per structural category (all genes)", "Sample",
                     "Percentage", color_map=category_color_palette, order=cat_order,
                     facet_col=facet),
        "Distribution of reads across SQANTI3 structural categories. A high FSM fraction "
        "indicates reads matching known transcripts; high ISM can signal degradation/5'–3' "
        "truncation; NIC/NNC capture novelty.", "fig-cat-all"))

    # 2. Structural category % (annotated genes)
    sections.append(_section(
        "Reads per structural category — annotated genes",
        _stacked_bar(annot_gene_percs_pivot_DF, "sampleID",
                     _value_cols(annot_gene_percs_pivot_DF, drop),
                     "Percent reads per structural category (annotated genes)", "Sample",
                     "Percentage", color_map=category_color_palette, order=cat_order,
                     facet_col=facet),
        "As above but restricted to annotated genes, so samples are compared on the same "
        "gene set.", "fig-cat-annot"))

    # 2b. Per-gene distribution of structural-category % (violin per sample)
    if melted_annotated_gene_DF is not None and not melted_annotated_gene_DF.empty:
        mdf = melted_annotated_gene_DF.copy()
        mdf["category"] = (mdf["category"].astype(str)
                           .str.replace("percent_", "", regex=False)
                           .str.replace("_", " ", regex=False))
        # Cap per (sample, category) so the box spread stays intact without
        # embedding every gene's value (this figure was the second-largest).
        mdf = _cap_per_group(mdf, ["sampleID", "category"])
        vio = go.Figure()
        for i, s in enumerate(samples):
            d = mdf[mdf["sampleID"] == s]
            vio.add_trace(go.Violin(x=d["category"], y=d["percentage"], name=str(s),
                                    legendgroup=str(s), scalegroup=str(s),
                                    box_visible=True, meanline_visible=True, points=False,
                                    marker_color=_DEFAULT_SEQ[i % len(_DEFAULT_SEQ)]))
        _base_layout(vio, "Per-gene structural-category composition", "Structural category",
                     "% of a gene's reads")
        vio.update_layout(violinmode="group")
        sections.append(_section(
            "Per-gene structural-category distribution", vio,
            "For each annotated gene, the % of its reads in each category; the box shows the "
            "spread across genes per sample. Wide/​shifted boxes flag heterogeneity between "
            "samples.", "fig-cat-box"))

    # 3. Genes detected (%) by read-count
    sections.append(_section(
        "Genes detected, by read support",
        _stacked_bar(gene_percs_unstacked, "sampleID",
                     _value_cols(gene_percs_unstacked, drop),
                     "Percent of detected genes by read-count bin", "Sample", "Percentage",
                     color_map=readcount_palette, order=readcount_order, facet_col=facet),
        "How much of each sample's detected genes are well-supported. A large '1 read' / "
        "'2-10 reads' fraction means many genes are seen only shallowly.", "fig-genes"))

    # 3b. Genes detected — absolute counts
    sections.append(_section(
        "Genes detected, by read support (counts)",
        _stacked_bar(gene_agg_DF, "sampleID", _value_cols(gene_agg_DF, drop),
                     "Number of detected genes by read-count bin", "Sample", "Number of genes",
                     color_map=readcount_palette, order=readcount_order, facet_col=facet),
        "Absolute number of genes detected in each read-support bin, per sample.", "fig-genes-n"))

    # 4. UJCs detected (%) by structural category
    if "structural_category" in ujc_percs_dct:
        ujc_sc = ujc_percs_dct["structural_category"]
        sections.append(_section(
            "Unique junction chains, by structural category",
            _stacked_bar(ujc_sc, "sampleID", _value_cols(ujc_sc, drop),
                         "Percent of UJCs per structural category", "Sample", "Percentage",
                         color_map=category_color_palette, order=cat_order, facet_col=facet),
            "Distribution of distinct splicing patterns (UJCs) across categories — complements "
            "the read-based view by weighting each junction chain once.", "fig-ujc"))

    # 4b. UJC counts by structural category
    if "structural_category" in ujc_cnts_dct:
        ujc_scn = ujc_cnts_dct["structural_category"]
        sections.append(_section(
            "Unique junction chains, by structural category (counts)",
            _stacked_bar(ujc_scn, "sampleID", _value_cols(ujc_scn, drop),
                         "Number of UJCs per structural category", "Sample", "Number of UJCs",
                         color_map=category_color_palette, order=cat_order, facet_col=facet),
            "Absolute number of distinct junction chains per category.", "fig-ujc-n"))

    # 4c. UJC by read support (%)
    if "read_category" in ujc_percs_dct:
        ujc_rc = ujc_percs_dct["read_category"]
        sections.append(_section(
            "UJCs, by read support",
            _stacked_bar(ujc_rc, "sampleID", _value_cols(ujc_rc, drop),
                         "Percent of UJCs by read-count bin", "Sample", "Percentage",
                         color_map=readcount_palette, order=readcount_order, facet_col=facet),
            "How well-supported the distinct junction chains are. Many single-read UJCs suggest "
            "noise or under-sequencing.", "fig-ujc-read"))

    # 5. Read length distribution (%)
    sections.append(_section(
        "Read length distribution",
        _stacked_bar(length_percs_agg, "sampleID", _value_cols(length_percs_agg, drop),
                     "Percent of reads by length bin", "Sample", "Percentage",
                     color_map=length_palette, order=length_order, facet_col=facet),
        "Length composition per sample. A large &lt;1kb fraction can indicate degradation or "
        "short-fragment libraries.", "fig-length"))

    # 5b. Raw read-length distribution (violin)
    if length_DF2 is not None and not length_DF2.empty and "length" in length_DF2.columns:
        sections.append(_section(
            "Read length distribution (raw)",
            _violin(length_DF2, "sampleID", "length", "Read length per sample", "Sample",
                    "Read length (bp)", facet_col=facet),
            "Per-read length distribution (box + density), with read length on the "
            "vertical axis. Compare medians and spread; a long tail reaching toward the "
            "bottom (low values) indicates an excess of short/degraded reads. For a compact "
            "interactive file, up to a few thousand reads per sample are drawn at random — the "
            "density is representative of the full set.",
            "fig-length-violin"))

    # 5c. Total reads vs % reads > 1kb (scatter)
    if {"perc_reads_gt_1kb", "total_reads"}.issubset(length_DF.columns):
        sections.append(_section(
            "Sequencing depth vs read length",
            _scatter_xy(length_DF, "perc_reads_gt_1kb", "total_reads", "sampleID",
                        "Total reads vs % reads &gt; 1kb", "% reads &gt; 1kb", "Total reads",
                        facet_col=facet),
            "Each sample by depth (total reads) against the fraction of long reads. Outliers on "
            "either axis are worth a closer look.", "fig-reads-len"))

    # 6. Subcategory breakdowns (% of reads), one per major category
    for label, sub_df, div in [("FSM", FSM_perc_DF, "fig-sub-fsm"),
                               ("ISM", ISM_perc_DF, "fig-sub-ism"),
                               ("NIC", NIC_perc_DF, "fig-sub-nic"),
                               ("NNC", NNC_perc_DF, "fig-sub-nnc")]:
        if sub_df is None or sub_df.empty:
            continue
        # NNC shares sub-category names with NIC but is colored red (matching its
        # structural-category color), so it uses a dedicated reds palette.
        sub_palette = nnc_subcat_color_palette if label == "NNC" else subcat_color_palette
        sub_order = list(sub_palette.keys())
        sections.append(_section(
            f"{label} reads by sub-category",
            _stacked_bar(sub_df, "sampleID", _value_cols(sub_df, drop),
                         f"Percent of {label} reads per sub-category", "Sample", "Percentage",
                         color_map=sub_palette, order=sub_order, facet_col=facet),
            f"Sub-category composition of {label} reads (e.g. reference match, alternative "
            f"5'/3' ends, fragments, intron retention). Shifts here pinpoint the kind of "
            f"partial/novel structure driving the category.", div))

    # 6a2. ISM fragment fraction in cohort context: one scalar per sample (share
    # of ISM reads that are 3'/5'/internal fragments) against the cohort median +
    # scorecard flags. A SHAPE signal, read relative to peers -- not a quality
    # verdict (fragment-shaped ISMs can be alternative isoforms, not degradation).
    _frag = _ism_fragment_pct(ISM_DF)
    if _frag:
        _fdf = pd.DataFrame({"sampleID": list(_frag.keys()),
                             "perc_ISM_fragments": list(_frag.values())})
        _ffig = _metric_cohort_figure(
            _fdf, "perc_ISM_fragments", "ISM fragment fraction (shape, not quality)",
            "% of ISM reads that are 3'/5'/internal fragments", scorecard=scorecard,
            threshold=qc_flags.get("perc_ISM_fragments"), color="#c4531d")
        if _ffig is not None:
            sections.append(_section(
                "ISM fragment fraction", _ffig,
                "Of each sample's incomplete-splice-match reads, the share that are "
                "3'/5'/internal fragments (rather than mono-exon or intron-retention "
                "ISMs). This describes the <em>shape</em> of the ISM population, not its "
                "quality: fragment-shaped ISMs arise from truncation or degradation but "
                "equally from genuine alternative isoforms (alternative TSS/TTS, shorter "
                "isoforms). Read it cohort-relatively — a sample whose ISM shape diverges "
                "from its peers is worth a look; a high value shared by the whole cohort "
                "is likely just the biology of the sample type, which is why this metric "
                "carries no absolute threshold.",
                "fig-ism-frag"))

    # 6b. Junction classification (known/novel × canonical/non-canonical)
    jxn_order = ["known_canonical", "known_non_canonical", "novel_canonical", "novel_non_canonical"]
    if nov_can_perc_DF is not None and not nov_can_perc_DF.empty:
        sections.append(_section(
            "Junction classification (%)",
            _stacked_bar(nov_can_perc_DF, "sampleID", _value_cols(nov_can_perc_DF, drop),
                         "Percent of junctions by class", "Sample", "Percentage",
                         color_map=jxn_palette, order=jxn_order, facet_col=facet),
            "Splice junctions split into known/novel and canonical/non-canonical. A rise in "
            "novel non-canonical junctions can indicate mapping noise or genuine novelty.",
            "fig-jxn"))
    if nov_can_DF is not None and not nov_can_DF.empty:
        sections.append(_section(
            "Junction classification (counts)",
            _stacked_bar(nov_can_DF, "sampleID", _value_cols(nov_can_DF, drop),
                         "Number of junctions by class", "Sample", "Number of junctions",
                         color_map=jxn_palette, order=jxn_order, facet_col=facet),
            "Absolute junction counts per class.", "fig-jxn-n"))

    # 6c. Novel non-canonical junction burden in cohort context: the single
    # junction class most likely to be alignment/calling artefacts, pulled out of
    # the classification above into one comparative scalar per sample.
    _nnc = _novel_noncanonical_pct(nov_can_DF)
    if _nnc:
        _ndf = pd.DataFrame({"sampleID": list(_nnc.keys()),
                             "perc_novel_noncanonical_jxn": list(_nnc.values())})
        _nfig = _metric_cohort_figure(
            _ndf, "perc_novel_noncanonical_jxn", "Novel non-canonical junction burden",
            "% of junctions that are novel & non-canonical", scorecard=scorecard,
            threshold=qc_flags.get("perc_novel_noncanonical_jxn"), color="#FC8D59")
        if _nfig is not None:
            sections.append(_section(
                "Novel non-canonical junction burden", _nfig,
                "Percent of all junctions that are both novel and non-canonical. This "
                "class is enriched for alignment and junction-calling artefacts, but it "
                "is not exclusively artefactual — non-canonical splicing does occur "
                "biologically (e.g. minor-spliceosome introns). The dashed grey line is the "
                "cohort median; the dotted and dash-dot lines are the configured <b>warn</b> "
                "and <b>fail</b> thresholds — fixed absolute cut-offs from the QC config "
                "(independent of the cohort), so a bar above them exceeds that preset level, "
                "while bars are outlined when the outlier scorecard flags the sample. Read it "
                "cohort-relatively: a sample well above its peers is worth inspecting for a "
                "mapping or library difference, while a level shared across the cohort may "
                "simply reflect the sample type or a less complete annotation.", "fig-jxn-novnc"))

    # 7. Artefacts (grouped bar of %)
    art_cols = ["perc_reads_RTS", "perc_reads_intrapriming", "perc_reads_non-canonical"]
    sections.append(_section(
        "Potential artefacts",
        _stacked_bar(err_DF, "sampleID", [c for c in art_cols if c in err_DF.columns],
                     "Reads with artefact evidence (%)", "Sample", "Percentage",
                     barmode="group", facet_col=facet),
        "Percent of reads flagged for RT-switching, intra-priming (high genomic-A downstream) "
        "or non-canonical junctions. Consistently elevated values across samples may indicate "
        "a systematic library issue.", "fig-artefacts"))

    # 7b. Per-sample cohort-context pages for the two artefact rates that feed the
    # outlier scorecard: each rate against the cohort median, config warn/fail
    # thresholds, and scorecard flag markers (bars outlined/colored when flagged).
    _art_pages = [
        ("perc_reads_RTS", "RT-switching reads per sample",
         "% reads with RT-switching evidence", "#F58A53", "fig-rts",
         "Per-sample rate of reads flagged by the RT-switching heuristic (a sequence "
         "signature of template switching at direct repeats around a junction). "
         "The dashed grey line is the cohort median; the dotted and dash-dot lines are the "
         "configured <b>warn</b> and <b>fail</b> thresholds — fixed absolute cut-offs from "
         "the QC config (not cohort-relative), so a bar rising above them exceeds that "
         "preset level. Bars are also outlined/marked when the outlier scorecard flags this "
         "sample. The flag marks candidate artefacts, not confirmed ones, so a sample well "
         "above its peers is worth checking for a library or mapping difference."),
        ("perc_reads_intrapriming", "Intra-priming reads per sample",
         "% reads with intra-priming evidence", "#EE446F", "fig-intraprim",
         "Per-sample intra-priming rate (reads with a genomic poly-A run downstream of "
         "the TTS, likely primed off genomic A-stretches rather than the polyA tail). "
         "The dashed grey line is the cohort median; the dotted and dash-dot lines are the "
         "configured warn / fail thresholds (fixed QC-config cut-offs, not cohort-relative). "
         "Bars are marked when the scorecard flags this sample. A sample elevated relative to "
         "its peers is worth checking for a sample-specific priming or template difference."),
    ]
    for metric, title, ytitle, color, div_id, interp in _art_pages:
        _mfig = _metric_cohort_figure(err_DF, metric, title, ytitle,
                                      scorecard=scorecard,
                                      threshold=qc_flags.get(metric), color=color)
        if _mfig is not None:
            sections.append(_section(title, _mfig, interp, div_id))

    # 8. Junction donor/acceptor CV summary (%)
    if cv_don_percs is not None and not cv_don_percs.empty:
        sections.append(_section(
            "Splice-donor detection variation",
            _stacked_bar(cv_don_percs, "sampleID", _value_cols(cv_don_percs, drop),
                         "Detected donors by position-variation category (%)", "Sample",
                         "Percentage", facet_col=facet, labels=_CV_PERC_LABELS),
            "For every reference splice donor, we look at where the reads that use it actually "
            "place the site relative to the annotated position, and classify the donor by the "
            "spread of those read offsets (using their coefficient of variation, CV): "
            "<b>ref_match</b> — every read lands exactly on the annotated donor (offset 0); "
            "<b>CV=0</b> — the reads land at a single, consistent position that is <i>shifted</i> "
            "from the annotation (a reproducible offset — no spread, but not on the reference); "
            "<b>CV&gt;0</b> — the reads land at variable positions (a real spread of offsets). "
            "So ref_match and CV=0 are <i>not</i> the same: both are consistent, but ref_match "
            "is on-target while CV=0 is consistently off-target; CV&gt;0 is genuinely imprecise. "
            "The bars show the % of a sample's donors in each category — a larger CV&gt;0 share "
            "means noisier donor boundaries.", "fig-cvdon"))

    # 8b. Splice-acceptor CV (%)
    if cv_acc_percs is not None and not cv_acc_percs.empty:
        sections.append(_section(
            "Splice-acceptor detection variation",
            _stacked_bar(cv_acc_percs, "sampleID", _value_cols(cv_acc_percs, drop),
                         "Detected acceptors by position-variation category (%)", "Sample",
                         "Percentage", facet_col=facet, labels=_CV_PERC_LABELS),
            "The acceptor counterpart of the donor plot above: every reference splice acceptor "
            "is classified by how consistently the reads that use it land relative to the "
            "annotated position — <b>ref_match</b> (all reads exactly on the acceptor), "
            "<b>CV=0</b> (all reads at one consistent but shifted position), <b>CV&gt;0</b> "
            "(reads at variable positions). A larger CV&gt;0 share means noisier acceptor "
            "boundaries; ref_match vs CV=0 distinguishes on-target from a reproducible offset.",
            "fig-cvacc"))

    # 8c. Detected donor / acceptor counts
    det_order = ["ref_match", "cv_0", "cv_gt_0"]
    for label, cv_df, div in [("donors", cv_don_summary, "fig-det-don"),
                              ("acceptors", cv_acc_summary, "fig-det-acc")]:
        if cv_df is None or cv_df.empty:
            continue
        sections.append(_section(
            f"Detected splice {label} (counts)",
            _stacked_bar(cv_df, "sampleID", [c for c in det_order if c in cv_df.columns],
                         f"Number of detected {label}", "Sample", f"Number of {label}",
                         facet_col=facet, labels=_CV_COUNT_LABELS),
            f"The absolute number of reference splice {label} detected in each sample, split by "
            f"the same position-consistency categories as the variation plot above: "
            f"<b>ref_match</b> (reads exactly on the annotated {label[:-1]}), <b>CV=0</b> (reads "
            f"at one consistent but shifted position) and <b>CV&gt;0</b> (reads at variable "
            f"positions). This is the count view — it shows how many {label} were recovered "
            f"overall and how that total splits across the categories, whereas the percentage "
            f"plot above normalises each sample to 100%.", div))

    # 9. PCA of QC metrics (colored by factor when set)
    if pca_DF is not None and pca_DF.shape[0] >= 2 and 0 in pca_DF.columns and 1 in pca_DF.columns:
        vr = list(variance_ratio) if variance_ratio is not None else [0, 0]
        fig = go.Figure()
        if facet and exp_factor in pca_DF.columns and pca_DF[exp_factor].nunique() > 1:
            for lv in pd.unique(pca_DF[exp_factor]):
                d = pca_DF[pca_DF[exp_factor] == lv]
                fig.add_trace(go.Scatter(
                    x=d[0], y=d[1], mode="markers+text", name=f"{exp_factor}={lv}",
                    text=d["sampleID"], textposition="top center", marker=dict(size=13),
                    hovertemplate="%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>"))
            show_leg = True
        else:
            fig.add_trace(go.Scatter(
                x=pca_DF[0], y=pca_DF[1], mode="markers+text", text=pca_DF["sampleID"],
                textposition="top center", marker=dict(size=13, color="#1f6feb"),
                hovertemplate="%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>"))
            show_leg = False
        _base_layout(fig, "PCA of per-sample QC metrics",
                     f"PC1 ({vr[0]*100:.1f}%)" if len(vr) > 0 else "PC1",
                     f"PC2 ({vr[1]*100:.1f}%)" if len(vr) > 1 else "PC2")
        fig.update_layout(barmode="group", showlegend=show_leg)
        sections.append(_section(
            "Sample similarity (PCA)", fig,
            "Each sample is described by its full vector of QC metrics; the metrics are "
            "standardised and reduced by <b>principal component analysis</b> to the two axes "
            "(PC1, PC2) that capture the most variation across samples, and each sample is "
            "plotted at its score on those axes (the % on each axis is the share of total "
            "variance it explains). Samples close together have similar overall QC profiles — "
            "replicates should cluster; a sample sitting apart differs systematically across "
            "many metrics at once.", "fig-pca"))

        # 9b. PCA loadings (which metrics drive PC1/PC2)
        if loadings_DF is not None and loadings_DF.shape[0] > 0 and {0, 1}.issubset(loadings_DF.columns):
            load = loadings_DF[[0, 1]].copy()
            load.columns = ["PC1", "PC2"]
            # order metrics by strongest absolute PC1 loading for readability
            load = load.reindex(load["PC1"].abs().sort_values(ascending=True).index)
            sections.append(_section(
                "PCA loadings",
                _heatmap(load, "Metric loadings on PC1 / PC2", "Principal component", "QC metric"),
                "Principal components are weighted combinations of the QC metrics; this heatmap "
                "shows each metric's weight (loading) on PC1 and PC2. The metrics with the "
                "largest absolute loadings are the ones driving how samples separate along each "
                "axis — read it alongside the PCA plot to see which metrics an outlier sample is "
                "extreme on.", "fig-pca-load"))

    # 9c. UJC-level metrics: saturation, replicate concordance, UpSet
    if ujc_metrics:
        m_samples = ujc_metrics["samples"]

        sat = ujc_metrics["saturation"]
        fig = go.Figure()
        for i, s in enumerate(m_samples):
            d = sat[sat["sampleID"] == s]
            fig.add_trace(go.Scatter(x=d["depth"], y=d["unique_ujcs"],
                                     mode="lines+markers", name=str(s),
                                     marker=dict(size=4),
                                     line=dict(color=sample_seq[i % len(sample_seq)])))
        _base_layout(fig, "UJC saturation (rarefaction)", "Reads sampled",
                     "Expected unique UJCs")
        fig.update_layout(barmode="group")
        sections.append(_section(
            "UJC saturation", fig,
            "Rarefaction curve: the expected number of distinct junction chains recovered as "
            "reads are subsampled to increasing depths. The expectation is computed "
            "analytically from the per-UJC read counts (the expected count of UJCs seen at "
            "least once at each depth — deterministic, not a random resample). A curve that "
            "plateaus means the library is saturated (more depth finds few new isoforms); a "
            "still-rising curve means deeper sequencing would keep discovering junction chains.",
            "fig-saturation"))

        # Saturation per structural category (same total-depth x-axis)
        sat_c = ujc_metrics.get("saturation_by_category")
        if sat_c is not None and not sat_c.empty:
            present = set(sat_c["structural_category"])
            cats = [c for c in cat_order if c in present]
            cats += [c for c in sorted(present) if c not in cats]
            ncols = min(3, len(cats))
            nrows = (len(cats) + ncols - 1) // ncols
            cfig = make_subplots(rows=nrows, cols=ncols, subplot_titles=cats,
                                 shared_xaxes=False, vertical_spacing=0.10,
                                 horizontal_spacing=0.06)
            for k, cat in enumerate(cats):
                r, cc = k // ncols + 1, k % ncols + 1
                cd = sat_c[sat_c["structural_category"] == cat]
                for i, s in enumerate(m_samples):
                    d = cd[cd["sampleID"] == s]
                    cfig.add_trace(go.Scatter(
                        x=d["depth"], y=d["unique_ujcs"], mode="lines", name=str(s),
                        legendgroup=str(s), showlegend=(k == 0),
                        line=dict(color=sample_seq[i % len(sample_seq)])), row=r, col=cc)
            cfig.update_layout(template="plotly_white", plot_bgcolor=_PLOT_BG,
                               paper_bgcolor=_PLOT_BG, height=260 * nrows + 80,
                               margin=dict(l=60, r=30, t=50, b=50),
                               legend=dict(title="Sample", x=1.02, y=1))
            cfig.update_xaxes(title_text="Reads sampled")
            cfig.update_yaxes(title_text="Unique UJCs")
            sections.append(_section(
                "UJC saturation by structural category", cfig,
                "The overall saturation split by structural category (same total-depth "
                "x-axis, so the panels decompose the curve above). Known categories (FSM) "
                "typically plateau early; novel ones (NIC/NNC) often keep rising, i.e. more "
                "depth is still discovering novel junction chains.", "fig-saturation-cat"))

        conc = ujc_metrics["concordance"]
        sections.append(_section(
            "Replicate concordance",
            _heatmap(conc, "Replicate concordance (per-UJC read counts)", "Sample",
                     "Sample", colorscale="Viridis", zmid=None, zmin=0, zmax=1,
                     colorbar_title="Pearson r", annotate=True),
            "Each unique junction chain's read count is tabulated per sample (a UJC × sample "
            "matrix), and every pair of sample columns is compared with the <b>Pearson "
            "correlation</b> of those counts — how linearly the two samples' per-UJC abundances "
            "track each other (1 = identical profile). Replicates of the same condition should "
            "correlate highly; a low value flags an outlier or a swapped/mismatched sample.",
            "fig-concordance"))

        up_data = compute_upset_intersections(ujc_metrics["upset"], m_samples)
        if up_data:
            sections.append(_section(
                "Shared UJCs across samples (UpSet)",
                _upset_figure(up_data, m_samples),
                "Which junction chains are shared across samples. Each column is a "
                "combination of samples (filled dots = members); the bar shows how many "
                "UJCs fall in that combination, stacked by structural category. Large "
                "sample-specific bars indicate low overlap between samples.", "fig-upset"))

        _ovfig = _ujc_overlap_figure(overlap_metrics)
        if _ovfig is not None:
            sections.append(_section(
                "Pairwise UJC overlap index", _ovfig,
                "For every ordered pair of samples we count how many unique junction chains "
                "(UJCs) they share and divide by the <b>row</b> sample's own UJC count, giving "
                "<b>|A ∩ B| / |A|</b> — the fraction of sample A's junction-chain repertoire that "
                "is also detected in sample B (diagonal = 1). Read along a row to see how much of "
                "that sample's repertoire each other sample recovers. The matrix is <i>asymmetric</i>: "
                "a small, deeply-shared sample can have a high row (most of its UJCs appear "
                "elsewhere) yet a low column (it recovers only a fraction of a larger sample's "
                "repertoire), which is exactly how it distinguishes shallow or low-complexity "
                "samples from rich ones. This complements the UpSet view (exact sharing pattern) "
                "and the read-count concordance (abundance agreement).", "fig-ujc-overlap"))

    # 9d. Splice-site fuzziness: offset spectrum, precision profile, canonical split
    if jxn_offset_metrics and jxn_offset_metrics.get("samples"):
        f_samples = jxn_offset_metrics["samples"]
        sfig, pfig, cfig = _jxn_offset_figures(jxn_offset_metrics, f_samples)
        if sfig is not None:
            sections.append(_section(
                "Splice-site offset spectrum", sfig,
                "For every detected donor and acceptor we measure its <b>signed distance to "
                "the nearest annotated reference site</b> — 0 = exactly on the reference, "
                "negative = a few bp upstream, positive = downstream — and count how many "
                "observations fall at each offset (donors and acceptors shown separately). "
                "The curve is that histogram of offsets. Exact matches (offset 0) are excluded "
                "so the plot shows only the imprecise tail, and only sites within ±window bp of "
                "a reference are included (farther ones are genuinely novel, not imprecise). "
                "A tight peak hugging 0 means the boundaries the reads place are very close to "
                "the annotation; a broad or one-sided skirt flags systematic imprecision (a "
                "consistent shift); and regular sub-peaks at ±3/±4/±6 bp are the signature of "
                "tandem-site (e.g. NAGNAG) usage — reads legitimately choosing a nearby "
                "alternative site.", "fig-offset-spectrum"))
        if pfig is not None:
            sections.append(_section(
                "Splice-site precision profile", pfig,
                "The same offsets as the spectrum above, read as a cumulative curve: for each "
                "distance k (in bp) it plots the <b>% of a sample's donor/acceptor observations "
                "that land within ±k bp of a reference site</b>. It starts at k=0 (the "
                "exact-match rate) and climbs to 100% by ±window bp. A steeper curve that "
                "reaches 100% sooner means splice boundaries sit closer to the annotation "
                "(more precise); a sample whose curve lags below the others has a larger "
                "fraction of its sites placed a few bp off the annotated position. It is the "
                "cumulative (ECDF) companion to the spectrum — easier for comparing samples at "
                "a glance.", "fig-offset-profile"))
        if cfig is not None:
            sections.append(_section(
                "Splice-site fuzziness by canonical class", cfig,
                "Share of site observations off by >0 bp, split by canonical vs "
                "non-canonical splice sites. Non-canonical junctions are typically far "
                "fuzzier; a high canonical fuzziness is unusual and worth inspecting.",
                "fig-offset-byclass"))

    # 9e. Read-end completeness profiles (5'/3' distance-to-gene-end ECDFs)
    if completeness_metrics and completeness_metrics.get("samples"):
        cmp_fig = _completeness_figure(completeness_metrics, completeness_metrics["samples"])
        if cmp_fig is not None:
            sections.append(_section(
                "Read-end completeness profile", cmp_fig,
                "Cumulative % of reads whose 5'/3' end lands within a given distance of "
                "the annotated gene end. Steeper curves mean more reads reach the annotated "
                "boundary. A sample whose 5' curve lags its peers has systematically shorter "
                "5' ends — consistent with truncation or RNA degradation, but also with "
                "alternative transcription start sites or incomplete reference annotation, "
                "so treat it as a comparative signal rather than a degradation verdict. "
                "3' ends are usually tight (polyA-anchored), so 5' completeness is the more "
                "informative axis for comparing samples.", "fig-completeness"))

    # 9f. Between-sample comparison views (A4 depth-normalised junction yield)
    if yield_metrics and yield_metrics.get("samples"):
        _yfig = _metric_cohort_figure(
            yield_metrics["per_sample"], "jxn_per_1k_reads",
            "Junction yield per 1000 reads (depth-normalised)",
            "distinct junctions / 1k reads", scorecard=scorecard, color="#74CDF0")
        if _yfig is not None:
            sections.append(_section(
                "Junction yield vs depth", _yfig,
                "Distinct junctions recovered per 1000 reads. Because it is normalised to "
                "depth, a sample low here has genuinely lower junction complexity "
                "(consistent with degradation or a simpler transcriptome), not merely fewer "
                "reads. Cohort-relative, no absolute threshold.", "fig-jxn-yield"))

    _dfig = _composition_drift_figure(drift_metrics)
    if _dfig is not None:
        sections.append(_section(
            "Composition drift (pairwise)", _dfig,
            "Each sample's structural-category percentages (FSM, ISM, NIC, NNC, …) are read "
            "as a probability distribution, and every pair of samples is compared with the "
            "<b>Jensen–Shannon distance</b> — a symmetric 0–1 measure of how different two "
            "distributions are (0 = identical category mix, 1 = no categories in common). The "
            "heatmap shows those pairwise distances. A sample far from all the others (a bright "
            "row/column — a high average distance to the rest of the cohort) has an atypical "
            "composition; worth checking whether that is a genuine biological difference "
            "(e.g. a distinct condition) or a technical one. The distance itself is not a "
            "pass/fail.", "fig-comp-drift"))

    _tfig = _tandem_sites_figure(tandem_metrics)
    if _tfig is not None:
        sections.append(_section(
            "Tandem splice sites (NAGNAG)", _tfig,
            "This zooms into the near-zero part of the splice-site offset distribution (the same "
            "signed offsets as the spectrum above) and highlights the bars at <b>±3 bp</b> "
            "(and ±4/±6). A NAGNAG tandem site is a pair of splice sites exactly 3 nt apart, so "
            "when reads legitimately use the alternative site their offsets pile up at ±3 bp — "
            "a sharp excess there (rather than a smooth decay away from 0) is the signature of "
            "tandem donor/acceptor usage, a real biological phenomenon distinct from random "
            "boundary imprecision. The note under the plot gives each sample's <b>tandem "
            "fraction</b> = the share of its imprecise (offset ≠ 0) observations that sit at "
            "±3 bp. This is shown descriptively — a higher tandem fraction is a property of the "
            "sample's biology, not a quality defect, so it carries no threshold.", "fig-tandem"))

    _rfig = _replicate_concordance_figure(rep_concordance)
    if _rfig is not None:
        sections.append(_section(
            "Replicate concordance (multi-axis)", _rfig,
            "For samples sharing a design-factor level (replicates), each sample is scored "
            "0–1 on how well it matches its group-mates on three axes, where <b>1 = identical "
            "to its replicates</b>. <b>Composition</b>: 1 − the mean Jensen–Shannon distance "
            "of its structural-category mix to each group-mate. <b>Read-length profile</b>: "
            "1 − the mean relative difference of its (q25, median, q75) read lengths to each "
            "group-mate. <b>Splice-site precision</b>: 1 − the mean relative difference of its "
            "imprecise-site fraction to each group-mate. A replicate low on one axis diverges "
            "from its peers on that specific property — more sensitive than UJC overlap alone. "
            "Shown only when the design has a factor with ≥2 replicates in a level. "
            "<b>Hover</b> a bar to also see the group's <b>robust spread (MAD)</b> on that axis — "
            "the median absolute deviation across the group of a per-sample scalar (composition: "
            "Jensen–Shannon distance to the group-mean mix, in JSD; length: median read length, "
            "in bp; imprecision: % imprecise sites, in pp), where 0 means the replicates are "
            "identical on that axis. The MAD is a whole-group companion to the per-sample "
            "agreement and stays informative even for exactly two replicates. "
            "<i>Note:</i> in a group of exactly two replicates each sample's only group-mate is "
            "the other, and all three distances are symmetric, so the two bars are identical by "
            "construction (the score reads as the pair's mutual agreement); values diverge "
            "per-sample only with ≥3 replicates in a group.",
            "fig-rep-concord"))

    _ffcfig = _fuzziness_concordance_figure(fuzz_concordance)
    if _ffcfig is not None:
        sections.append(_section(
            "Replicate concordance of splice-site precision", _ffcfig,
            "Whether replicates place the <i>same</i> reference splice sites at the <i>same</i> "
            "sub-bp offsets. For each reference donor/acceptor we take each sample's median "
            "signed offset (how many bp its reads land from the annotated site), keep the sites "
            "seen in ≥2 replicates of a group, and score each sample as <b>1 − the mean absolute "
            "offset difference to its group-mates, capped at ±3 bp</b> (so 1 = places every "
            "shared site identically, 0 = ≥3 bp apart on average). High agreement means the "
            "imprecision pattern is reproducible (a property of the sites/protocol); a replicate "
            "that disagrees has sample-specific boundary placement. Shown only with a factor and "
            "≥2 replicates in a level.", "fig-fuzz-concord"))

    _fdfig = _fuzz_depth_figure(fuzz_depth_metrics)
    if _fdfig is not None:
        sections.append(_section(
            "Splice-site imprecision vs junction depth", _fdfig,
            "Share of imprecise splice-site observations as a function of the junction's "
            "read depth. Imprecision concentrated at low-depth junctions is consistent "
            "with limited support (few reads to pin the boundary); imprecision that "
            "persists at high depth points to a systematic placement offset. Requires "
            "per-junction coverage (<code>total_coverage_unique</code>); the section is "
            "omitted when that column is not populated.", "fig-fuzz-depth"))

    # 9e. Orthogonal support (B11/B12/B13) — CAGE 5', polyA 3', short-read. Each
    # section renders only when that evidence was supplied to SQANTI3; the fixture
    # has none, so all three are omitted there.
    if cage_metrics and cage_metrics.get("samples"):
        cfig = _metric_cohort_figure(
            cage_metrics["per_sample"], "perc_within_cage",
            "CAGE 5'-end support: reads with a TSS inside a CAGE peak",
            "% of reads within a CAGE peak", color="#74CDF0")
        if cfig is not None:
            sections.append(_section(
                "CAGE 5'-end support (B11)", cfig,
                "For each sample we take the fraction of reads whose 5' end (TSS) falls inside a "
                "CAGE peak — an orthogonal, cap-based measurement of genuine transcription start "
                "sites — and show it against the cohort median. This is independent evidence that "
                "a sample's 5' ends are real rather than truncated or degraded; read "
                "cohort-relatively, a sample well below its peers has less CAGE-supported 5' ends. "
                "Requires CAGE data supplied to SQANTI3 (the <code>within_CAGE_peak</code> column); "
                "omitted otherwise.", "fig-cage"))
    if polya_metrics and polya_metrics.get("samples"):
        pfig = _metric_cohort_figure(
            polya_metrics["per_sample"], "perc_within_polya",
            "PolyA 3'-end support: reads with a TTS at a polyA site",
            "% of reads at a polyA site", color="#EE446F")
        if pfig is not None:
            sections.append(_section(
                "PolyA 3'-end support (B12)", pfig,
                "For each sample we take the fraction of reads whose 3' end (TTS) sits at a known "
                "polyA site — orthogonal, cleavage-based evidence of a genuine transcript end — "
                "shown against the cohort median. Read cohort-relatively: a sample below its peers "
                "has fewer independently-supported 3' ends (consistent with internal priming or 3' "
                "degradation, worth inspecting alongside the intra-priming rate). Requires polyA "
                "data supplied to SQANTI3 (the <code>within_polyA_site</code> column); omitted "
                "otherwise.", "fig-polya"))
    if sr_metrics and sr_metrics.get("samples"):
        sfig = _metric_cohort_figure(
            sr_metrics["per_sample"], "perc_jxn_SR_supported",
            "Short-read junction support: junctions with short-read coverage",
            "% of junctions with short-read coverage", color="#15918A")
        if sfig is not None:
            sections.append(_section(
                "Short-read junction support (B13)", sfig,
                "For each sample we take the fraction of junctions that carry Illumina short-read "
                "coverage (<code>total_coverage_unique &gt; 0</code>) — orthogonal, high-accuracy "
                "confirmation that a splice junction is real — shown against the cohort median. "
                "Read cohort-relatively: a sample below its peers has more junctions lacking "
                "short-read backing (worth checking against its novel-junction and RT-switching "
                "rates). Requires short-read coverage supplied to SQANTI3; omitted otherwise.",
                "fig-sr-jxn"))

    # 10. Under-annotation section (from CSV on disk)
    sections.append(_gene_classification_section(args.OUT, args.PREFIX))

    # 11. Part C — closing optional-support note (only when some type is absent).
    sections.append(_optional_support_footer(cage_metrics, polya_metrics, sr_metrics))

    factor_note = (f" &middot; faceted by <code>{args.inFACTOR}</code>"
                   if getattr(args, "inFACTOR", None) else "")
    html = PAGE_TEMPLATE.format(
        title=f"{args.PREFIX} — SQANTI-reads report",
        plotlyjs=get_plotlyjs(),
        n_samples=len(samples), prefix=args.PREFIX, factor_note=factor_note,
        pass_c=_FLAG_COLORS["pass"], warn_c=_FLAG_COLORS["warn"], fail_c=_FLAG_COLORS["fail"],
        body="\n".join(sections),
    )
    with open(out_path, "w") as fh:
        fh.write(html)
    reads_logger.info(f"Interactive HTML report saved as {out_path}")

    # qc_summary.json sidecar
    json_path = os.path.join(args.OUT, args.PREFIX + "_qc_summary.json")
    summary = {
        "prefix": args.PREFIX,
        "n_samples": len(samples),
        "factor": getattr(args, "inFACTOR", None),
        "thresholds": qc_flags,
        "samples": per_sample,
    }
    # Cohort-relative sample-outlier verdict (machine-readable), when computed.
    if scorecard is not None:
        if scorecard.get("enabled"):
            summary["sample_scorecard"] = {
                "enabled": True,
                "metrics": scorecard["metrics"],
                "z": scorecard["z"],
                "cell_flags": scorecard["cell_flags"],
                "n_flagged": scorecard["n_flagged"],
                "overall": scorecard["overall"],
            }
        else:
            summary["sample_scorecard"] = {"enabled": False,
                                           "reason": scorecard.get("reason")}
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    reads_logger.info(f"QC summary written to {json_path}")
    return out_path
