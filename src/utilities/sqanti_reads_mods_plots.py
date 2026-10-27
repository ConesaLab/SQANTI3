import numpy as np
import sys
import os
import pandas as pd 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from src.module_logging import reads_logger
from src.utilities.conesa_colors import color_conesa, conesa_cmap

#CONSTANTS
CAT_ABBR = {"full-splice_match": "FSM", "incomplete-splice_match": "ISM",
            "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC",
            "antisense": "AS", "fusion": "FUS", "genic": "GENIC",
            "genic_intron": "GI", "intergenic": "INTER"}
CAT_ORDER = ["FSM", "ISM", "NIC", "NNC", "AS", "FUS", "GENIC", "GI", "INTER"]
CAT_PALETTE = {"FSM": "#6BAED6", "ISM": "#FC8D59", "NIC": "#78C679", "NNC": "#EE6A50",
               "GENIC": "#969696", "AS": "#66C2A4", "FUS": "#ffc125",
               "INTER": "#e9967a", "GI": "#41B6C4"}

READ_COLS = ["isoform", "chrom", "read_start", "read_end", "strand", "length", "exons",
             "structural_category", "subcategory", "associated_gene",
             "associated_transcript", "jxn_string", "jxnHash"]

METAGENE_BINS = 100
def load_reads(classif_path, sampleID, args):
    # We load the classification table with the information we need
    n_mod_cols = [f"{mod}_n_mod" for mod in args.MODTYPES]
    # Clear error if the file comes from a version without these columns (e.g. before read_start/read_end)
    header = pd.read_csv(classif_path, sep="\t", nrows=0).columns
    missing = [col for col in READ_COLS + n_mod_cols if col not in header]
    if missing:
        reads_logger.error(f"[ERROR] {classif_path} lacks the columns {', '.join(missing)}")
        sys.exit(-1)

    df = pd.read_csv(classif_path, sep="\t", usecols=READ_COLS + n_mod_cols,
                     dtype=str, keep_default_na=False)
    # For reads with coordinates only
    df = df[df["read_start"] != ""].copy()
    # We convert to numeric type
    num_cols = ["read_start", "read_end", "length", "exons"] + n_mod_cols
    df[num_cols] = df[num_cols].astype(int)
    # We abbreviate structural categories
    df["category"] = df["structural_category"].map(CAT_ABBR)
    # Number of mods per kb
    for mod in args.MODTYPES:
        df[f"{mod}_per_kb"] = df[f"{mod}_n_mod"] / df["length"] * 1000
    # We add a column with the sampleID
    df.insert(0, "sampleID", sampleID)
    return df

def load_sites(sites_path, reads, sampleID):
    '''One row per modified site of a UJC of a sample, with the attributes of its UJC'''
    sites = pd.read_csv(sites_path, sep="\t")
    # Adding useful information from the classification file
    ujc_info = reads.drop_duplicates("jxnHash")[["jxnHash", "category", "exons",
                                                 "associated_gene", "jxn_string"]]
    sites = sites.merge(ujc_info, on="jxnHash", how="left")
    sites.insert(0, "sampleID", sampleID)
    return sites

def load_design_samples(design_file, args):
    design = pd.read_csv(design_file, sep=",")
    all_reads, all_sites = [] , []
    for _, row in design.iterrows():
        reads_logger.info(f"Loading modifications of {row['sampleID']}")
        reads = load_reads(row["mod_classification_file"], row["sampleID"], args)
        all_reads.append(reads)
        all_sites.append(load_sites(row["ujc_mod_sites_file"], reads, row["sampleID"]))
    return pd.concat(all_reads,ignore_index=True), pd.concat(all_sites, ignore_index=True)

def filter_sites(sites, min_cov=20, min_mod=1):
    '''Sites covered by at least min_cov reads of their UJC and modified in at least min_mod of them'''
    return sites[(sites["n_cov"] >= min_cov) & (sites["n_mod"] >= min_mod)]

def filter_label(n_sites, min_cov, min_mod):
    # ≥ is ≥, written as an escape so the title does not depend on the file encoding
    return f"n_cov ≥ {min_cov}, n_mod ≥ {min_mod}, {n_sites} sites"

# To group by number of exons
EXON_BIN_EDGES = [0, 1, 2, 5, 10, np.inf]
EXON_BIN_LABELS = ["1", "2", "3-5", "6-10", ">10"]

def add_exon_bin(df):
    '''Groups the number of exons in the bins used by the plots'''
    df = df.copy()
    df["exon_bin"] = pd.cut(df["exons"], EXON_BIN_EDGES, labels=EXON_BIN_LABELS).astype(str)
    return df

def _group_palette(df, group_col, sample_palette=None):
    # With one sample the bars are coloured by group: SQANTI colours for the structural categories,
    # Conesa colours for the exon bins. With several samples each sample gets its own colour
    if df["sampleID"].nunique() == 1:
        if group_col == "category":
            return group_col, CAT_PALETTE
        return group_col, dict(zip(EXON_BIN_LABELS, color_conesa(len(EXON_BIN_LABELS))))
    return "sampleID", sample_palette


def _panel_grid(n_panels, ncols=4, width=4, height=3.5, **kwargs):
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows),
                             squeeze=False, **kwargs)
    for i in range(n_panels, nrows * ncols):
        axes.flat[i].set_visible(False)
        # With a shared x axis the panel above an empty one would be left without tick labels
        if i - ncols >= 0:
            axes.flat[i - ncols].tick_params(labelbottom=True)
    return fig, list(axes.flat)


def plot_reads_by_group(reads, modtypes, group_col="category", order=None, metric="per_kb", palette=None):
    '''
    Per read modification metric by group (structural category or exon bin), a panel per type.
    metric="per_kb": mean modifications per kb of read; metric="pct": % of reads with >=1 modification
    palette: colour per sample, used when there are several samples
    '''
    order = order or [g for g in CAT_ORDER if g in reads[group_col].unique()]
    hue, palette = _group_palette(reads, group_col, palette)
    fig, axes = _panel_grid(len(modtypes), sharex=True)
    for ax, mod in zip(axes, modtypes):
        data = reads[["sampleID", group_col]].copy()
        data["value"] = (reads[f"{mod}_per_kb"] if metric == "per_kb"
                         else (reads[f"{mod}_n_mod"] > 0) * 100)
        sns.barplot(data=data, x=group_col, y="value", hue=hue, order=order, palette=palette,
                    errorbar="se", ax=ax, legend=(hue == "sampleID" and ax is axes[0]))
        ax.set_title(mod)
        ax.set_xlabel("")
        ax.set_ylabel("Mods per kb (mean ± se)" if metric == "per_kb" else "% reads with ≥1 mod")
        ax.tick_params(axis="x", rotation=90)
    title = "Modification density" if metric == "per_kb" else "Reads with modifications"
    fig.suptitle(f"{title} by {group_col}")
    fig.tight_layout()
    return fig


def plot_occupancy_by_group(sites, modtypes, group_col="category", order=None, min_cov=10, min_mod=1,
                            palette=None):
    '''Occupancy of UJC sites by group (structural category or exon bin), a panel per type.
    palette: colour per sample, used when there are several samples'''
    covered = filter_sites(sites, min_cov, min_mod)
    order = order or [g for g in CAT_ORDER if g in covered[group_col].unique()]
    mods = [mod for mod in modtypes if mod in covered["mod_type"].unique()]
    hue, palette = _group_palette(covered, group_col, palette)
    fig, axes = _panel_grid(len(mods), sharex=True, sharey=True)
    for ax, mod in zip(axes, mods):
        sns.boxplot(data=covered[covered["mod_type"] == mod], x=group_col, y="occupancy",
                    hue=hue, order=order, palette=palette, showfliers=False, ax=ax,
                    legend=(hue == "sampleID" and ax is axes[0]))
        ax.set_title(mod)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=90)
    fig.suptitle(f"Site occupancy by {group_col} ({filter_label(len(covered), min_cov, min_mod)})")
    fig.tight_layout()
    return fig

def parse_introns(jxn_string):
    '''(start, end) of the introns of a UJC, 1-based and inclusive'''
    parts = jxn_string.split("_")
    if "monoexon" in parts:
        return []
    # We look for the position of the strand
    strand_idx = next(i for i, part in enumerate(parts) if part in ("+", "-"))
    # Coords will be a list of all the intron starting positions
    coords = list(map(int, parts[strand_idx + 1:]))
    # Now we return a list of tuples (start,end) of the introns of the UJC
    return list(zip(coords[::2], coords[1::2]))

def exon_model(jxn_string, ujc_start, ujc_end):
    '''Genomic exons of a UJC: We go from introns to exon'''
    # We generate a list of tuples (start,end) of the introns
    introns = parse_introns(jxn_string)
    # The exon start is the intron_end+1 for every intron 
    ex_starts = np.array([ujc_start] + [end + 1 for _, end in introns])
    # The exon end is the intron start -1 for every intron
    ex_ends = np.array([start - 1 for start, _ in introns] + [ujc_end])
    return ex_starts, ex_ends

def to_ujc_coords(pos, ex_starts, ex_ends, strand):
    '''
    Takes a position in a read and gives 0-based position from the 5' end of 
    the UJC model and exon index in 5'->3' order. Non exonic positions get -1
    '''
    # lengths of all the exons
    lengths = ex_ends - ex_starts + 1
    # Where do exons 1, 2 ,3,...  start?
    # Exon 1 starts at 0 and the rest at the cumulated sum of exons
    offsets = np.concatenate(([0], np.cumsum(lengths)[:-1]))
    # We look for the last exon that starts before the position
    idx = np.searchsorted(ex_starts, pos, side="right") - 1
    # Cannot be higher than the length of exons - 1
    idx_c = idx.clip(0, len(ex_starts) - 1)
    # We check that pos is really in the exon
    exonic = (idx >= 0) & (pos <= ex_ends[idx_c])
    # bps of previous exons + distance inside the exon
    tx = offsets[idx_c] + pos - ex_starts[idx_c]
    total = int(lengths.sum())
    # If strand is - 
    if strand == "-":
        tx = total - 1 - tx
        idx_c = len(ex_starts) - 1 - idx_c
    return np.where(exonic, tx, -1), idx_c, total

def junction_distance(tx, exon_idx, lengths_5to3):
    '''
    Auxiliar function for add_ujc_position
    The input is the position of each site in the UJC (taken from
    to_ujc_coords), the exon in which every site falls and the length of each
    exon.
    It computes the distance to the nearest exon-exon junction
    '''
    # 
    n_ex = len(lengths_5to3)
    if n_ex == 1:
        
        return np.full(len(tx), np.nan)
    # Offset is the start of each exon
    offsets = np.concatenate(([0], np.cumsum(lengths_5to3)[:-1]))
    # Distance to previous junction
    up = np.where(exon_idx > 0, tx - offsets[exon_idx] + 1, np.inf)
    # Distance to next junction
    down = np.where(exon_idx < n_ex - 1, offsets[exon_idx] + lengths_5to3[exon_idx] - tx, np.inf)
    # If up is nearer it keeps those distances, if its down it is negative
    dist = np.where(up <= down, up, -down).astype(float)
    return np.where(tx >= 0, dist, np.nan)


def junction_background(lengths_5to3):
    '''
    Auxiliar function for add_ujc_position
    Counts number of nts at each side of the junctions of a UJC. Recieves
    length of exons from 5' to 3'
    Returns 2 lists that help us correct for different exon distances
    '''
    n_ex = len(lengths_5to3) 
    if n_ex == 1:
        return [], []
    # Down_side staets with last exon and up with first
    down_side, up_side = [lengths_5to3[-1]], [lengths_5to3[0]]
    for length in lengths_5to3[1:-1]:
        # Internal exons are split at their middle, ties go to the upstream junction
        down_side.append((length + 1) // 2)
        up_side.append(length - (length + 1) // 2)
    return down_side, up_side


def add_ujc_positions(sites, reads, n_bins=METAGENE_BINS, max_dist=100):
    '''
    Returns 3 tables:
    - Sites with added columns (tx_pos: pos fron 5' in UJC, ujc_length, rel_pos: relative position
    from 0 (5' end) to 100 (3' end), dist_3p, dist_junction and junction_side)
    -Coverage (for each sample and stretch of 1/n_bins of rel_pos, what fraction of reads covers it)
    - Background: For each sample and distance from -100 to 100, how many possible exonic positions are there
    '''
    sites = sites.reset_index(drop=True).copy()
    tx_all = np.full(len(sites), -1)
    length_all = np.full(len(sites), -1)
    dist_all = np.full(len(sites), np.nan)
    coverage, background = [], []
    # For each sample 
    for sampleID, s_reads in reads.groupby("sampleID"):
        
        s_rows = np.flatnonzero((sites["sampleID"] == sampleID).to_numpy())
        site_groups = sites.iloc[s_rows].groupby("jxnHash").indices
        site_pos = sites["mod_pos"].to_numpy()[s_rows]
        r_start, r_end = s_reads["read_start"].to_numpy(), s_reads["read_end"].to_numpy()
        r_strand, r_jxn = s_reads["strand"].to_numpy(), s_reads["jxn_string"].to_numpy()
        cov_diff = np.zeros(n_bins + 1)
        down_lengths, up_lengths = [], []
        # For each exon 
        for jxn, r_idx in s_reads.groupby("jxnHash").indices.items():
            
            strand = r_strand[r_idx[0]]
            ex_starts, ex_ends = exon_model(r_jxn[r_idx[0]], r_start[r_idx].min(), r_end[r_idx].max())
            lengths = ex_ends - ex_starts + 1
            lengths_5to3 = lengths if strand == "+" else lengths[::-1]

            # Read coverage along the UJC, as +1/-1 at the bins where each read starts/ends
            t1, _, total = to_ujc_coords(r_start[r_idx], ex_starts, ex_ends, strand)
            t2, _, _ = to_ujc_coords(r_end[r_idx], ex_starts, ex_ends, strand)
            scale = n_bins / max(total, 1)
            b5 = (np.minimum(t1, t2) * scale).astype(int).clip(0, n_bins - 1)
            b3 = (np.maximum(t1, t2) * scale).astype(int).clip(0, n_bins - 1)
            np.add.at(cov_diff, b5, 1)
            np.add.at(cov_diff, b3 + 1, -1)

            if jxn not in site_groups:
                continue
            down, up = junction_background(lengths_5to3)
            down_lengths += down
            up_lengths += up
            rows = site_groups[jxn]
            tx, exon_idx, total = to_ujc_coords(site_pos[rows], ex_starts, ex_ends, strand)
            tx_all[s_rows[rows]] = tx
            length_all[s_rows[rows]] = total
            dist_all[s_rows[rows]] = junction_distance(tx, exon_idx, lengths_5to3)

        coverage.append(pd.DataFrame({"sampleID": sampleID,
                                      "rel_pos": (np.arange(n_bins) + 0.5) * 100 / n_bins,
                                      "read_fraction": np.cumsum(cov_diff)[:n_bins] / len(s_reads)}))
        # Positions at distance d: stretches of length >= d, counted with a reverse cumulative sum
        dists = np.arange(1, max_dist + 1)
        n_down = np.bincount(np.minimum(np.asarray(down_lengths, dtype=int), max_dist),
                             minlength=max_dist + 1)[::-1].cumsum()[::-1][1:]
        n_up = np.bincount(np.minimum(np.asarray(up_lengths, dtype=int), max_dist),
                           minlength=max_dist + 1)[::-1].cumsum()[::-1][1:]
        background.append(pd.DataFrame({"sampleID": sampleID,
                                        "dist_junction": np.concatenate((-dists[::-1], dists)),
                                        "n_positions": np.concatenate((n_up[::-1], n_down))}))

    sites["tx_pos"] = tx_all
    sites["ujc_length"] = length_all
    sites["rel_pos"] = np.where(tx_all >= 0, tx_all * 100 / np.maximum(length_all - 1, 1), np.nan)
    sites["dist_3p"] = np.where(tx_all >= 0, length_all - 1 - tx_all, -1)
    sites["dist_junction"] = dist_all
    sites["junction_side"] = np.select([dist_all < 0, dist_all > 0], ["donor", "acceptor"], "none")
    return sites, pd.concat(coverage, ignore_index=True), pd.concat(background, ignore_index=True)

def plot_junction_profile(sites, background, modtypes, max_dist=100, min_cov=1, min_mod=1, smooth=5,
                          metric="density", palette=None):
    '''
    A panel per type and a line per sample (palette: colour per sample).
    -if metric="density": sites per exonic position available at each distance, relative to the mean of
    the window (+- max_dist). It counts detected sites, so it also rises where coverage is higher.
    -elif metric="occupancy": mean occupancy of the sites at each distance. 
    '''
    modtypes = modtypes or sorted(sites["mod_type"].unique())
    # We filter sites to a minimum coverage and number of mods
    near = filter_sites(sites, min_cov, min_mod)
    # We filter out all sites that our further away than the selected max_dist
    near = near[near["dist_junction"].abs() <= max_dist]
    mods = [mod for mod in modtypes if mod in near["mod_type"].unique()]
    
    fig, axes = _panel_grid(len(mods))
    # One panel for each mod
    for ax, mod in zip(axes, mods):
        # For each sample
        for sampleID, bg in background.groupby("sampleID"):
            # We group sites from the sample that have the mod by distance to junction
            by_dist = near[(near["mod_type"] == mod) & (near["sampleID"] == sampleID)].groupby("dist_junction")
            # Distances to junction
            dists = bg["dist_junction"]
            if metric == "density":
                # Number of counts for each distance is the size of the group
                counts = by_dist.size().reindex(dists, fill_value=0).to_numpy()
                # Density is the number of counts divided by possible position at that distance
                density = counts / np.maximum(bg["n_positions"].to_numpy(), 1)
                values = pd.Series(density / density.mean()).rolling(smooth, center=True, min_periods=1).mean()
            else:
                # Smoothed as a weighted mean: occupancy summed over the window / sites in the window
                occ_sum = by_dist["occupancy"].sum().reindex(dists, fill_value=0).reset_index(drop=True)
                n_sites = by_dist.size().reindex(dists, fill_value=0).reset_index(drop=True)
                values = (occ_sum.rolling(smooth, center=True, min_periods=1).sum()
                          / n_sites.rolling(smooth, center=True, min_periods=1).sum())
                # Reference line at the mean occupancy of the type in the whole window
                ax.axhline(occ_sum.sum() / max(n_sites.sum(), 1), color="grey", lw=0.8, ls=":")
            ax.plot(dists, values, label=sampleID, color=palette[sampleID] if palette else None)
        ax.axvline(0, color="grey", lw=0.8, ls="--")
        if metric == "density":
            ax.axhline(1, color="grey", lw=0.8, ls=":")
        ax.set_xlim(-max_dist, max_dist)
        ax.set_title(mod)
        ax.set_ylabel("Relative site density" if metric == "density" else "Mean occupancy")
    axes[0].legend(title="sampleID", fontsize=8)
    # One x label for the whole figure, a long label per panel overlaps the panels below
    fig.supxlabel("nt to the nearest exon-exon junction (- donor side, + acceptor side)")
    what = "Modified sites" if metric == "density" else "Mean site occupancy"
    fig.suptitle(f"{what} around exon-exon junctions ({filter_label(len(near), min_cov, min_mod)})")
    fig.tight_layout()
    return fig


def plot_metagene(sites, coverage, modtypes, n_bins=METAGENE_BINS, min_cov=10, min_mod=1, multiexon_only=False,
                  palette=None):
    '''
    Sites along the 5'->3' UJC model (0-100). Rows: share of the sites of each type per bin,
    mean occupancy per bin (sites with enough coverage) and fraction of reads covering each bin.
    A column per sample (palette: colour per modification type)
    '''
    data = sites[sites["rel_pos"].notna()]
    if multiexon_only:
        data = data[data["exons"] > 1]
    hue_order = [mod for mod in modtypes if mod in data["mod_type"].unique()]
    edges = np.linspace(0, 100, n_bins + 1)
    data = data.assign(bin_center=pd.cut(data["rel_pos"], edges, include_lowest=True,
                                         labels=(edges[:-1] + edges[1:]) / 2).astype(float))
    samples = data["sampleID"].unique()
    fig, axes = plt.subplots(3, len(samples), figsize=(8 * len(samples), 11), sharex=True, squeeze=False)
    for col, sampleID in enumerate(samples):
        s = data[data["sampleID"] == sampleID]
        share = (s.groupby(["mod_type", "bin_center"]).size()
                  .groupby(level="mod_type").transform(lambda x: x / x.sum()).rename("share").reset_index())
        sns.lineplot(data=share, x="bin_center", y="share", hue="mod_type", hue_order=hue_order,
                     palette=palette, ax=axes[0, col])
        occ = (filter_sites(s, min_cov, min_mod).groupby(["mod_type", "bin_center"])["occupancy"]
               .mean().rename("mean_occupancy").reset_index())
        sns.lineplot(data=occ, x="bin_center", y="mean_occupancy", hue="mod_type", hue_order=hue_order,
                     palette=palette, ax=axes[1, col], legend=False)
        cov = coverage[coverage["sampleID"] == sampleID]
        axes[2, col].plot(cov["rel_pos"], cov["read_fraction"], color="black")
        axes[0, col].set_title(sampleID)
        axes[0, col].set_ylabel("Share of the sites of the type")
        axes[1, col].set_ylabel(f"Mean occupancy (n_cov ≥ {min_cov}, n_mod ≥ {min_mod})")
        axes[2, col].set_ylabel("Fraction of reads covering")
        axes[2, col].set_xlabel("Position along the UJC (% from 5' to 3')")
    fig.suptitle("Modifications along the UJC" + (" (multiexonic UJCs)" if multiexon_only else ""))
    fig.tight_layout()
    return fig

def plot_dist_3p(sites, modtypes, max_nt=3000, bin_nt=100, min_cov=10, min_mod=1, palette=None):
    '''Sites by distance to the 3' end of their UJC: number of sites and mean occupancy per bin.
    More stable than rel_pos in direct RNA, where the 3' ends of the reads vary much less than the 5' ends.
    palette: colour per modification type'''
    data = sites[(sites["dist_3p"] >= 0) & (sites["dist_3p"] < max_nt)]
    data = data.assign(bin_nt=(data["dist_3p"] // bin_nt) * bin_nt + bin_nt / 2)
    hue_order = [mod for mod in modtypes if mod in data["mod_type"].unique()]
    samples = data["sampleID"].unique()
    fig, axes = plt.subplots(2, len(samples), figsize=(8 * len(samples), 8), sharex=True, squeeze=False)
    for col, sampleID in enumerate(samples):
        s = data[data["sampleID"] == sampleID]
        counts = s.groupby(["mod_type", "bin_nt"]).size().rename("n_sites").reset_index()
        sns.lineplot(data=counts, x="bin_nt", y="n_sites", hue="mod_type", hue_order=hue_order,
                     palette=palette, ax=axes[0, col])
        occ = (filter_sites(s, min_cov, min_mod).groupby(["mod_type", "bin_nt"])["occupancy"]
               .mean().rename("mean_occupancy").reset_index())
        sns.lineplot(data=occ, x="bin_nt", y="mean_occupancy", hue="mod_type", hue_order=hue_order,
                     palette=palette, ax=axes[1, col], legend=False)
        axes[0, col].set_yscale("log")
        axes[0, col].set_title(sampleID)
        axes[0, col].set_ylabel("Number of sites")
        axes[1, col].set_ylabel(f"Mean occupancy (n_cov ≥ {min_cov}, n_mod ≥ {min_mod})")
        axes[1, col].set_xlabel("nt from the 3' end of the UJC")
        axes[1, col].invert_xaxis()
    fig.suptitle("Modifications by distance to the 3' end")
    fig.tight_layout()
    return fig

def plot_sample_correlation(sites, modtypes, min_cov=10):
    '''Spearman correlation of the occupancy of the sites shared by all samples with n_cov >= min_cov,
    a heatmap per type. Returns None with less than two samples'''
    if sites["sampleID"].nunique() < 2:
        reads_logger.info("Only one sample, the correlation between samples needs at least two")
        return None
    covered = sites[sites["n_cov"] >= min_cov]
    occ = covered.pivot_table(index=["jxnHash", "mod_type", "mod_pos"], columns="sampleID",
                              values="occupancy").dropna()
    mods = [mod for mod in modtypes if mod in occ.index.get_level_values("mod_type")]
    fig, axes = _panel_grid(len(mods), width=4.5, height=4)
    for ax, mod in zip(axes, mods):
        corr = occ.xs(mod, level="mod_type").corr(method="spearman")
        n_shared = len(occ.xs(mod, level="mod_type"))
        sns.heatmap(corr, vmin=0, vmax=1, annot=True, fmt=".2f", cmap=conesa_cmap("warm"), ax=ax)
        ax.set_title(f"{mod} ({n_shared} shared sites)")
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle(f"Occupancy correlation between samples (n_cov ≥ {min_cov} in all)")
    fig.tight_layout()
    return fig

def run_mods_plots(design_file, out_dir, prefix, args, min_cov=20, min_mod=2):
    '''Modification report: loads the reads and UJC sites of every sample of the processed design,
    adds their positions along the UJCs and writes <prefix>_mods_plots.pdf in out_dir'''
    reads, sites = load_design_samples(design_file, args)
    reads_logger.info("Placing the modified sites along their UJCs...")
    sites_pos, coverage, junction_bg = add_ujc_positions(sites, reads)
    reads = add_exon_bin(reads)
    sites_pos = add_exon_bin(sites_pos)
    modtypes = args.MODTYPES
    # Fixed colours
    samples = sorted(reads["sampleID"].unique())
    mod_palette = dict(zip(modtypes, color_conesa(len(modtypes))))
    sample_palette = dict(zip(samples, color_conesa(len(samples))))
    figures = [
        # B. Structural category
        ("category % reads", lambda: plot_reads_by_group(reads, modtypes, "category", metric="pct",
                                                         palette=sample_palette)),
        ("category per kb", lambda: plot_reads_by_group(reads, modtypes, "category", metric="per_kb",
                                                        palette=sample_palette)),
        ("category occupancy", lambda: plot_occupancy_by_group(sites_pos, modtypes, "category", min_cov=min_cov,
                                                               min_mod=min_mod, palette=sample_palette)),
        # C. Number of exons
        ("exons % reads", lambda: plot_reads_by_group(reads, modtypes, "exon_bin", order=EXON_BIN_LABELS, metric="pct",
                                                      palette=sample_palette)),
        ("exons per kb", lambda: plot_reads_by_group(reads, modtypes, "exon_bin", order=EXON_BIN_LABELS,
                                                     metric="per_kb", palette=sample_palette)),
        ("exons occupancy", lambda: plot_occupancy_by_group(sites_pos, modtypes, "exon_bin", order=EXON_BIN_LABELS,
                                                            min_cov=min_cov, min_mod=min_mod, palette=sample_palette)),
        # D-E. Junctions, metagene and 3' end
        ("junction density", lambda: plot_junction_profile(sites_pos, junction_bg, modtypes, palette=sample_palette)),
        ("junction occupancy", lambda: plot_junction_profile(sites_pos, junction_bg, modtypes, min_cov=min_cov,
                                                             min_mod=min_mod, metric="occupancy", smooth=9,
                                                             palette=sample_palette)),
        ("metagene", lambda: plot_metagene(sites_pos, coverage, modtypes, min_cov=min_cov, min_mod=min_mod,
                                           palette=mod_palette)),
        ("distance to 3' end", lambda: plot_dist_3p(sites_pos, modtypes, min_cov=min_cov, min_mod=min_mod,
                                                    palette=mod_palette)),
        # F. Between samples (skipped with one sample)
        ("correlation between samples", lambda: plot_sample_correlation(sites_pos, min_cov=min_cov, modtypes=modtypes))
    ]
    pdf_path = os.path.join(out_dir, f"{prefix}_mods_plots.pdf")
    with PdfPages(pdf_path) as pdf:
        for name, make in figures:
            # A failing figure is reported and skipped, so it does not lose the whole report
            try:
                fig = make()
                if fig is None:
                    continue
                pdf.savefig(fig)
                plt.close(fig)
            except Exception as err:
                reads_logger.warning(f"[WARNING] Figure '{name}' could not be drawn: {err}")
    reads_logger.info(f"Modification plots written to {pdf_path}")
    return pdf_path
