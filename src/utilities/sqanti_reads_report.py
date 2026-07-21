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
    cat_order,
    jxn_palette,
    three_series_palette,
    readcount_palette,
    length_palette,
    sample_seq,
    compute_upset_intersections,
)

# Fallback qualitative colors for series that have no entry in an explicit palette.
# Prefix with the shared per-sample sequence so sample colors match the PDF.
_DEFAULT_SEQ = list(sample_seq) + pcolors.qualitative.Light24
# Preferred palette for exactly-three-series plots (green, orange, yellow) — shared.
_THREE_SEQ = three_series_palette

# --- QC flag thresholds --------------------------------------------------------
# Per-sample thresholds for the QC flag panel. "warn"/"fail" are lower bounds for
# metrics where higher is worse, except read length where lower is worse.
# NOTE: these defaults will move into the Phase-4 YAML config; edit here for now.
QC_THRESHOLDS = {
    "perc_reads_intrapriming": {"warn": 20.0, "fail": 40.0, "worse": "high",
                                 "label": "Intra-priming reads (%)"},
    "perc_reads_RTS": {"warn": 15.0, "fail": 30.0, "worse": "high",
                        "label": "RT-switching reads (%)"},
    "perc_reads_non-canonical": {"warn": 10.0, "fail": 20.0, "worse": "high",
                                  "label": "Non-canonical reads (%)"},
    "median_length": {"warn": 500.0, "fail": 250.0, "worse": "low",
                       "label": "Median read length (bp)"},
}

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
                 order=None, barmode="stack", facet_col=None):
    """Generic stacked/grouped bar: one trace per value column.

    When ``facet_col`` is a real design factor with >1 level, the plot is split
    into one subplot column per factor level (mirroring the PDF's FacetGrid),
    with a single shared legend.
    """
    cols = [c for c in (order or value_cols) if c in value_cols]
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
            fig.add_trace(go.Bar(name=str(c), x=x, y=d[c].tolist(),
                                 marker_color=col_color[c], legendgroup=str(c),
                                 showlegend=show_legend,
                                 hovertemplate=f"%{{x}}<br>{c}: %{{y:.2f}}<extra></extra>"),
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


def _violin(df, x_col, y_col, title, xtitle, ytitle, facet_col=None):
    """One violin (distribution) per x category; colored by factor when set."""
    fig = go.Figure()
    faceted = (facet_col and facet_col in df.columns
               and facet_col != "temp_factor" and df[facet_col].nunique() > 1)
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
        for yi, idx in enumerate(matrix_df.index):
            for xi, col in enumerate(matrix_df.columns):
                fig.add_annotation(x=str(col), y=str(idx),
                                   text=format(matrix_df.values[yi, xi], value_fmt),
                                   showarrow=False, font=dict(size=11, color="#222"))
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
                     thresholds=QC_THRESHOLDS):
    """Per-sample metrics + flags. Returns (summary_dict, samples_in_order)."""
    exp_factor = _factor_col(args)
    samples = length_DF["sampleID"].astype(str).tolist()

    length_idx = length_DF.set_index("sampleID")
    err_idx = err_DF.set_index("sampleID")
    fsm = all_gene_percs_pivot_DF.set_index("sampleID")["FSM"] \
        if "FSM" in all_gene_percs_pivot_DF.columns else None
    gene_cols = _value_cols(gene_agg_DF, {"sampleID", exp_factor})
    genes_detected = gene_agg_DF.set_index("sampleID")[gene_cols].sum(axis=1)

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
                 "% intra-prim", "% RTS", "% non-canon", "Overall"]
    rows = []
    for s in samples:
        m = per_sample[s]
        badge = f'<span class="badge {m["overall_flag"]}">{m["overall_flag"].upper()}</span>'
        rows.append(
            "<tr>" + "".join(f"<td>{v}</td>" for v in [
                s, f'{m["total_reads"]:,}', m["genes_detected"],
                f'{m["median_length"]:.0f}', f'{m["perc_reads_gt_1kb"]:.1f}',
                (f'{m["perc_FSM"]:.1f}' if m["perc_FSM"] is not None else "—"),
                f'{m["perc_reads_intrapriming"]:.1f}',
                f'{m["perc_reads_RTS"]:.1f}',
                f'{m["perc_reads_non-canonical"]:.1f}',
            ]) + f"<td>{badge}</td></tr>"
        )
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
    colors = {
        "annotated_with_well_covered_FSM": "#6BAED6",
        "annotated_with_low_coverage_FSM": "#9E9AC8",
        "underannotated_with_candidate_transcript": "#FD8D3C",
        "underannotated_no_candidate_transcripts": "#FDBE85",
    }
    fig = go.Figure(go.Bar(
        x=[c.replace("_", " ") for c in counts.index], y=counts.values,
        marker_color=[colors.get(c, "#969696") for c in counts.index],
        hovertemplate="%{x}<br>%{y} genes<extra></extra>"))
    _base_layout(fig, "Gene annotation categories", "Category", "Number of genes")
    fig.update_layout(showlegend=False, barmode="group")
    interp = ("Genes classified by whether a well-covered annotated (FSM) transcript was seen. "
              "<b>underannotated_with_candidate_transcript</b> highlights genes where a "
              "well-covered novel UJC suggests a missing annotation.")
    return _section("Under-annotation analysis", fig, interp, "fig-underannot")


PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script>{plotlyjs}</script>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
          margin: 0; background: #f4f6f8; color: #1a2027; }}
  header {{ background: #1f3a5f; color: #fff; padding: 20px 28px; }}
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


def build_html_report(out_path, dfs_for_plotting, args, ujc_metrics=None):
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
     cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_NNC_DF, FSM_perc_DF,
     ISM_perc_DF, NIC_NNC_perc_DF, nov_can_DF, nov_can_perc_DF, length_DF2,
     cv_acc_percs, cv_don_percs) = dfs_for_plotting

    exp_factor = _factor_col(args)
    drop = {"sampleID", exp_factor}
    facet = getattr(args, "inFACTOR", None)  # None -> single (non-faceted) plots

    # QC flag thresholds come from the (optionally user-supplied) config.
    qc_flags = args.cfg().get("qc_flags", QC_THRESHOLDS) if hasattr(args, "cfg") else QC_THRESHOLDS

    per_sample, samples = _compute_summary(length_DF, err_DF, all_gene_percs_pivot_DF,
                                           gene_agg_DF, args, thresholds=qc_flags)

    sections = [_summary_html(per_sample, samples, thresholds=qc_flags)]

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
            "Full per-read length distribution (box + density). Compare medians and spread; a "
            "long left tail indicates short/degraded reads.", "fig-length-violin"))

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
                               ("NIC/NNC", NIC_NNC_perc_DF, "fig-sub-nicnnc")]:
        if sub_df is None or sub_df.empty:
            continue
        sections.append(_section(
            f"{label} reads by sub-category",
            _stacked_bar(sub_df, "sampleID", _value_cols(sub_df, drop),
                         f"Percent of {label} reads per sub-category", "Sample", "Percentage",
                         color_map=subcat_color_palette, order=subcat_order, facet_col=facet),
            f"Sub-category composition of {label} reads (e.g. reference match, alternative "
            f"5'/3' ends, fragments, intron retention). Shifts here pinpoint the kind of "
            f"partial/novel structure driving the category.", div))

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

    # 8. Junction donor/acceptor CV summary (%)
    if cv_don_percs is not None and not cv_don_percs.empty:
        sections.append(_section(
            "Splice-donor detection variation",
            _stacked_bar(cv_don_percs, "sampleID", _value_cols(cv_don_percs, drop),
                         "Detected donors by position-variation category (%)", "Sample",
                         "Percentage", facet_col=facet),
            "Reference splice donors split by how variable their detected position is "
            "(ref_match / cv=0 / cv&gt;0). More cv&gt;0 means noisier junction ends.",
            "fig-cvdon"))

    # 8b. Splice-acceptor CV (%)
    if cv_acc_percs is not None and not cv_acc_percs.empty:
        sections.append(_section(
            "Splice-acceptor detection variation",
            _stacked_bar(cv_acc_percs, "sampleID", _value_cols(cv_acc_percs, drop),
                         "Detected acceptors by position-variation category (%)", "Sample",
                         "Percentage", facet_col=facet),
            "As for donors, but for reference splice acceptors.", "fig-cvacc"))

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
                         facet_col=facet),
            f"Absolute number of reference {label} detected, by position-variation category.",
            div))

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
            "Each point is a sample in the space of its QC metrics. Replicates should cluster; "
            "an outlier flags a sample that differs systematically.", "fig-pca"))

        # 9b. PCA loadings (which metrics drive PC1/PC2)
        if loadings_DF is not None and loadings_DF.shape[0] > 0 and {0, 1}.issubset(loadings_DF.columns):
            load = loadings_DF[[0, 1]].copy()
            load.columns = ["PC1", "PC2"]
            # order metrics by strongest absolute PC1 loading for readability
            load = load.reindex(load["PC1"].abs().sort_values(ascending=True).index)
            sections.append(_section(
                "PCA loadings",
                _heatmap(load, "Metric loadings on PC1 / PC2", "Principal component", "QC metric"),
                "How much each QC metric contributes to PC1/PC2. Metrics with the largest "
                "absolute loadings explain why samples separate along each axis.", "fig-pca-load"))

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
            "Rarefaction: expected distinct junction chains as reads are subsampled. "
            "A curve that plateaus means the library is saturated (more depth finds few "
            "new isoforms); a still-rising curve means deeper sequencing would help.",
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
            "Pearson correlation of per-UJC read counts between samples. Replicates of the "
            "same condition should correlate highly; a low value flags an outlier or a "
            "swapped/mismatched sample.", "fig-concordance"))

        up_data = compute_upset_intersections(ujc_metrics["upset"], m_samples)
        if up_data:
            sections.append(_section(
                "Shared UJCs across samples (UpSet)",
                _upset_figure(up_data, m_samples),
                "Which junction chains are shared across samples. Each column is a "
                "combination of samples (filled dots = members); the bar shows how many "
                "UJCs fall in that combination, stacked by structural category. Large "
                "sample-specific bars indicate low overlap between samples.", "fig-upset"))

    # 10. Under-annotation section (from CSV on disk)
    sections.append(_gene_classification_section(args.OUT, args.PREFIX))

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
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    reads_logger.info(f"QC summary written to {json_path}")
    return out_path
