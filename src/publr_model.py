"""
publr_model.py

Probability that a NIC/NNC intron chain is a real isoform, from the
pubLR model shipped inside the pantranscriptome catalog. pubLR is a
LightGBM classifier trained with positive-unlabeled learning on the
intron chains of public long-read samples: chains seen in at least 3
studies and 2 platforms are the positives. Its variables describe the
chain only (structure, splice motifs, PWM strength of the sites, novelty
and distance of the sites to the annotated ones, and direct repeats at
the junctions), so every multi-exon NIC/NNC can be scored without
reads. They are computed here exactly as in training.

Version: 1.0
Author: Adrian Berenguer-Agustí (I2SysBio-CSIC)
Date: 2026-11-07
"""

import gzip
import re

import lightgbm as lgb
import numpy as np
import pandas as pd

# Windows of the PWMs (as in MaxEntScan): donor 3 exonic + 6 intronic
# nt, acceptor 20 intronic + 3 exonic nt.
DONOR_EXON = 3
DONOR_INTRON = 6
ACCEPTOR_INTRON = 20
ACCEPTOR_EXON = 3
# Canonical and semi-canonical donor-acceptor dinucleotides, and their
# codes as computed by motif_codes().
CANONICAL_MOTIFS = ("GTAG", "GCAG", "ATAC")
CANONICAL_CODES = np.array([sum("ACGT".index(base) * 4 ** (3 - i)
                                for i, base in enumerate(motif))
                            for motif in CANONICAL_MOTIFS])
# Longest direct repeat measured on each side of a junction.
MAX_REPEAT = 30
# Novel sites within this distance (nt) of an annotated site count as
# shifted.
SHIFT_NT = 10
# Cap of the distance (nt) of a novel site to an annotated one.
MAX_DISTANCE = 100000
# Code of each base in the genome arrays (anything else is N = 4).
BASE_CODE = np.full(256, 4, dtype=np.uint8)
for _code, _base in enumerate(b"ACGT"):
    BASE_CODE[_base] = _code
    BASE_CODE[ord(chr(_base).lower())] = _code
# Complement of each base code.
COMPLEMENT_CODE = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
# Pattern of the intron chain key (see intron_chain_key()).
CHAIN_PATTERN = re.compile(r"^(.+?)([+-])(\d+-\d+(?:_\d+-\d+)*)$")
# SQANTI3 subcategories of NIC and NNC, coded as 0/1 variables.
SUBCATEGORIES = ("combination_of_known_junctions",
                 "combination_of_known_splicesites", "intron_retention",
                 "at_least_one_novel_splicesite",
                 "mono-exon_by_intron_retention")


def genome_codes(genome_dict: dict | None, fasta_path: str,
                 chroms: set[str]) -> dict[str, np.ndarray]:
    """Base codes of some chromosomes.

    Args:
        genome_dict: Chromosome name -> SeqRecord, as read by QC, or
            None to read the FASTA.
        fasta_path: Genome FASTA (plain or gzipped), read only when
            genome_dict is None.
        chroms: Chromosomes to keep.

    Returns:
        Base codes (0-3 for ACGT, 4 for N) of each chromosome found;
        index 0 is position 1.
    """
    genome = {}
    if genome_dict is not None:
        for name in chroms & set(genome_dict):
            genome[name] = BASE_CODE[np.frombuffer(
                bytes(genome_dict[name].seq), np.uint8)]
    else:
        opener = gzip.open if fasta_path.endswith(".gz") else open
        name = None
        chunks = []
        handle = opener(fasta_path, "rb")
        for line in handle:
            if line.startswith(b">"):
                if name in chroms:
                    genome[name] = BASE_CODE[np.frombuffer(
                        b"".join(chunks), np.uint8)]
                name = line[1:].split()[0].decode()
                chunks = []
            elif name in chroms:
                chunks.append(line.rstrip())
        if name in chroms:
            genome[name] = BASE_CODE[np.frombuffer(b"".join(chunks),
                                                   np.uint8)]
        handle.close()
    result = genome

    return result


def annotated_introns(gtf_path: str) -> pd.DataFrame:
    """Read the introns of an annotation GTF.

    Args:
        gtf_path: Annotation GTF (plain or gzipped).

    Returns:
        Unique introns with chrom, strand, start and end (1-based,
        inclusive).
    """
    exons = pd.read_csv(gtf_path, sep="\t", header=None, comment="#",
                        usecols=[0, 2, 3, 4, 6, 8],
                        names=["chrom", "feature", "start", "end", "strand",
                               "attr"])
    exons = exons[exons["feature"] == "exon"].copy()
    exons["transcript"] = exons["attr"].str.extract(
        r'transcript_id "([^"]+)"', expand=False)
    exons = exons.sort_values(["transcript", "start"])
    exons["next_start"] = exons.groupby("transcript")["start"].shift(-1)
    introns = exons.dropna(subset=["next_start"])
    result = pd.DataFrame({
        "chrom": introns["chrom"].to_numpy(),
        "strand": introns["strand"].to_numpy(),
        "start": introns["end"].to_numpy() + 1,
        "end": introns["next_start"].to_numpy().astype(np.int64) - 1,
    }).drop_duplicates().reset_index(drop=True)

    return result


def split_chains(chains: pd.Series) -> pd.DataFrame:
    """Split intron chain keys into one row per junction.

    Args:
        chains: Intron chain keys.

    Returns:
        One row per junction with columns ujc, chrom, strand, rank
        (0-based, genomic order), start and end (1-based, inclusive
        intron coordinates).

    Raises:
        ValueError: If a key does not match CHAIN_PATTERN.
    """
    parts = chains.str.extract(CHAIN_PATTERN)
    if parts.isna().any().any():
        raise ValueError("Malformed intron chain key: "
                         f"{chains[parts[0].isna()].iloc[0]}")
    parts.columns = ["chrom", "strand", "chain"]
    parts["ujc"] = chains.to_numpy()
    junctions = parts.assign(intron=parts["chain"].str.split("_")).explode(
        "intron")
    junctions["rank"] = junctions.groupby("ujc").cumcount()
    bounds = junctions["intron"].str.split("-", expand=True).astype(np.int64)
    # The key holds the last exonic base and the first base of the next
    # exon; the intron lies strictly between them.
    junctions["start"] = bounds[0] + 1
    junctions["end"] = bounds[1] - 1
    result = junctions[["ujc", "chrom", "strand", "rank", "start",
                        "end"]].reset_index(drop=True)

    return result


def window_codes(sequence: np.ndarray, anchors: np.ndarray, before: int,
                 after: int, minus: np.ndarray) -> np.ndarray:
    """Extract fixed windows around sites, in the sense of the gene.

    Args:
        sequence: Base codes of one chromosome (index 0 = position 1).
        anchors: 1-based positions of the first base after the window's
            "before" part, on the plus strand.
        before: Bases taken before the anchor.
        after: Bases taken from the anchor on.
        minus: Whether each site is on the minus strand; those windows
            are reverse-complemented.

    Returns:
        Array (sites x window) of base codes; positions off the
        chromosome are N.
    """
    offsets = np.arange(-before, after)
    index = anchors[:, None] - 1 + offsets[None, :]
    inside = (index >= 0) & (index < len(sequence))
    codes = np.where(inside, sequence[np.clip(index, 0, len(sequence) - 1)],
                     4).astype(np.uint8)
    codes[minus] = COMPLEMENT_CODE[codes[minus][:, ::-1]]
    result = codes

    return result


def site_windows(sequence: np.ndarray, junctions: pd.DataFrame
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Donor and acceptor windows of junctions of one chromosome.

    Args:
        sequence: Base codes of the chromosome.
        junctions: Junctions with strand, start and end.

    Returns:
        Donor windows (DONOR_EXON + DONOR_INTRON) and acceptor windows
        (ACCEPTOR_INTRON + ACCEPTOR_EXON), in the sense of the gene.
    """
    minus = (junctions["strand"] == "-").to_numpy()
    start = junctions["start"].to_numpy()
    end = junctions["end"].to_numpy()
    # On the plus strand the donor is at the intron start; on the minus
    # strand, at its end, and the window is then reverse-complemented,
    # so the anchor moves to keep the same exonic/intronic split.
    donor_plus = window_codes(sequence, start, DONOR_EXON, DONOR_INTRON,
                              np.zeros(len(start), bool))
    donor_minus = window_codes(sequence, end - DONOR_INTRON + 1, 0,
                               DONOR_INTRON + DONOR_EXON, minus)
    acceptor_plus = window_codes(sequence, end - ACCEPTOR_INTRON + 1,
                                 0, ACCEPTOR_INTRON + ACCEPTOR_EXON,
                                 np.zeros(len(start), bool))
    acceptor_minus = window_codes(sequence, start - ACCEPTOR_EXON,
                                  0, ACCEPTOR_INTRON + ACCEPTOR_EXON, minus)
    donors = np.where(minus[:, None], donor_minus, donor_plus)
    acceptors = np.where(minus[:, None], acceptor_minus, acceptor_plus)
    result = (donors, acceptors)

    return result


def repeat_length(sequence: np.ndarray, junctions: pd.DataFrame
                  ) -> np.ndarray:
    """Length of the direct repeat that makes each junction ambiguous.

    A junction can slide left while the last exonic base equals the
    last intronic base, and right while the first intronic base equals
    the first base of the next exon. The sum of both slides is the
    length of the direct repeat at the junction, which RT switching
    leaves behind.

    Args:
        sequence: Base codes of the chromosome.
        junctions: Junctions with start and end.

    Returns:
        Repeat length of each junction, at most 2 * MAX_REPEAT.
    """
    start = junctions["start"].to_numpy() - 1
    end = junctions["end"].to_numpy() - 1
    size = len(sequence)
    left = np.zeros(len(start), dtype=np.int64)
    right = np.zeros(len(start), dtype=np.int64)
    open_left = np.ones(len(start), dtype=bool)
    open_right = np.ones(len(start), dtype=bool)
    for shift in range(MAX_REPEAT):
        a = np.clip(start - 1 - shift, 0, size - 1)
        b = np.clip(end - shift, 0, size - 1)
        open_left &= (sequence[a] == sequence[b]) & (sequence[a] < 4)
        left += open_left
        c = np.clip(start + shift, 0, size - 1)
        d = np.clip(end + 1 + shift, 0, size - 1)
        open_right &= (sequence[c] == sequence[d]) & (sequence[c] < 4)
        right += open_right
    result = left + right

    return result


def motif_codes(donors: np.ndarray, acceptors: np.ndarray) -> np.ndarray:
    """Code the donor-acceptor dinucleotides of junctions.

    Args:
        donors: Donor windows, in the sense of the gene.
        acceptors: Acceptor windows, in the sense of the gene.

    Returns:
        Code of the four motif bases in base 4, or -1 if one is N.
    """
    bases = np.column_stack([donors[:, DONOR_EXON],
                             donors[:, DONOR_EXON + 1],
                             acceptors[:, ACCEPTOR_INTRON - 2],
                             acceptors[:, ACCEPTOR_INTRON - 1]]).astype(int)
    codes = bases @ np.array([64, 16, 4, 1])
    result = np.where((bases < 4).all(axis=1), codes, -1)

    return result


def pwm_score(windows: np.ndarray, pwm: np.ndarray) -> np.ndarray:
    """Score windows with a PWM (N scores 0).

    Args:
        windows: Base codes (sites x window).
        pwm: Log2-odds matrix (window x 4).

    Returns:
        Score of each window.
    """
    padded = np.vstack([pwm.T, np.zeros(pwm.shape[0])])
    result = padded[windows, np.arange(pwm.shape[0])[None, :]].sum(axis=1)

    return result


def nearest_distance(sites: np.ndarray, annotated: np.ndarray
                     ) -> np.ndarray:
    """Distance from each site to the nearest annotated site.

    Args:
        sites: Positions to look up.
        annotated: Sorted annotated positions (may be empty).

    Returns:
        Distance to the nearest annotated position, MAX_DISTANCE if
        there is none.
    """
    distance = np.full(len(sites), MAX_DISTANCE, dtype=np.int64)
    if len(annotated):
        index = np.searchsorted(annotated, sites)
        lower = annotated[np.clip(index - 1, 0, len(annotated) - 1)]
        upper = annotated[np.clip(index, 0, len(annotated) - 1)]
        distance = np.minimum(np.abs(sites - lower), np.abs(upper - sites))
        distance = np.minimum(distance, MAX_DISTANCE)
    result = distance

    return result


def junction_table(junctions: pd.DataFrame, genome: dict[str, np.ndarray],
                   introns: pd.DataFrame,
                   pwms: tuple[np.ndarray, np.ndarray]) -> pd.DataFrame:
    """Compute the per-junction variables.

    Args:
        junctions: Output of split_chains().
        genome: Base codes of the chromosomes.
        introns: Introns of the reference annotation.
        pwms: Donor and acceptor PWMs.

    Returns:
        The junctions with motif, PWM scores, repeat length, whether the
        junction is annotated and the distance of each site to the
        nearest annotated one.
    """
    table = junctions.merge(introns, on=["chrom", "strand", "start", "end"],
                            how="left", indicator=True)
    table["known_junction"] = table.pop("_merge") == "both"
    parts = []
    for (chrom, strand), group in table.groupby(["chrom", "strand"]):
        group = group.copy()
        sequence = genome.get(chrom)
        if sequence is None:
            sequence = np.full(int(group["end"].max()) + 1000, 4, np.uint8)
        donors, acceptors = site_windows(sequence, group)
        group["noncanonical"] = ~np.isin(motif_codes(donors, acceptors),
                                         CANONICAL_CODES)
        group["donor_score"] = pwm_score(donors, pwms[0])
        group["acceptor_score"] = pwm_score(acceptors, pwms[1])
        group["repeat"] = repeat_length(sequence, group)
        same = introns[(introns["chrom"] == chrom)
                       & (introns["strand"] == strand)]
        starts = np.unique(same["start"].to_numpy())
        ends = np.unique(same["end"].to_numpy())
        start_distance = nearest_distance(group["start"].to_numpy(), starts)
        end_distance = nearest_distance(group["end"].to_numpy(), ends)
        donor_is_start = strand == "+"
        group["donor_distance"] = (start_distance if donor_is_start
                                   else end_distance)
        group["acceptor_distance"] = (end_distance if donor_is_start
                                      else start_distance)
        parts.append(group)
    result = pd.concat(parts, ignore_index=True)

    return result


def chain_features(chains: pd.Series, genome: dict[str, np.ndarray],
                   introns: pd.DataFrame,
                   pwms: tuple[np.ndarray, np.ndarray]) -> pd.DataFrame:
    """Compute the variables of a set of intron chains.

    Args:
        chains: Unique intron chain keys.
        genome: Base codes of the chromosomes.
        introns: Introns of the reference annotation.
        pwms: Donor and acceptor PWMs.

    Returns:
        One row per chain (index ujc) with the structure, motif,
        strength, annotation and repeat variables.
    """
    table = junction_table(split_chains(chains), genome, introns, pwms)
    table["intron_log"] = np.log10(table["end"] - table["start"] + 1)
    table = table.sort_values(["ujc", "rank"])
    # Internal exon: from the end of one intron to the start of the next.
    table["next_start"] = table.groupby("ujc")["start"].shift(-1)
    table["internal_exon"] = table["next_start"] - table["end"] - 1
    novel_donor = table["donor_distance"] > 0
    novel_acceptor = table["acceptor_distance"] > 0
    table["novel_junction"] = ~table["known_junction"]
    table["novel_donor"] = novel_donor
    table["novel_acceptor"] = novel_acceptor
    table["shifted_sites"] = (
        (novel_donor & (table["donor_distance"] <= SHIFT_NT)).astype(int)
        + (novel_acceptor & (table["acceptor_distance"] <= SHIFT_NT)))
    # Distance of the closest novel site of the junction to an annotated
    # site of its kind (NaN if both sites are annotated).
    nearest = np.minimum(
        np.where(novel_donor, table["donor_distance"], np.inf),
        np.where(novel_acceptor, table["acceptor_distance"], np.inf))
    table["novel_site_log_distance"] = np.log10(
        1 + np.where(np.isinf(nearest), np.nan, nearest))
    table["novel_repeat"] = np.where(table["novel_junction"],
                                     table["repeat"], np.nan)
    groups = table.groupby("ujc")
    result = pd.DataFrame({
        "n_junctions": groups.size(),
        "intron_log_min": groups["intron_log"].min(),
        "intron_log_max": groups["intron_log"].max(),
        "intron_log_mean": groups["intron_log"].mean(),
        "internal_exon_min": groups["internal_exon"].min(),
        "n_noncanonical": groups["noncanonical"].sum(),
        "donor_score_min": groups["donor_score"].min(),
        "donor_score_mean": groups["donor_score"].mean(),
        "acceptor_score_min": groups["acceptor_score"].min(),
        "acceptor_score_mean": groups["acceptor_score"].mean(),
        "n_novel_junctions": groups["novel_junction"].sum(),
        "n_novel_donors": groups["novel_donor"].sum(),
        "n_novel_acceptors": groups["novel_acceptor"].sum(),
        "n_shifted_sites": groups["shifted_sites"].sum(),
        "novel_site_log_distance_min": groups[
            "novel_site_log_distance"].min(),
        "repeat_max": groups["repeat"].max(),
        "novel_repeat_max": groups["novel_repeat"].max(),
    })
    result["frac_novel_junctions"] = (result["n_novel_junctions"]
                                      / result["n_junctions"])

    return result


def model_probability(isoforms: pd.DataFrame, genome: dict[str, np.ndarray],
                      introns: pd.DataFrame, model: dict) -> np.ndarray:
    """pubLR probability of some NIC/NNC isoforms.

    Args:
        isoforms: One row per isoform with ujc, structural_category and
            subcategory.
        genome: Base codes of the chromosomes.
        introns: Introns of the reference annotation.
        model: booster (model text), pwm_donor, pwm_acceptor, features
            and c, as stored in the catalog.

    Returns:
        Probability of each isoform, in the order of the rows: the
        LightGBM score divided by the PU constant c, capped at 1.
    """
    probability = np.zeros(len(isoforms))
    if len(isoforms):
        chains = chain_features(isoforms["ujc"].drop_duplicates(), genome,
                                introns, (model["pwm_donor"],
                                          model["pwm_acceptor"]))
        features = chains.loc[isoforms["ujc"]].set_axis(isoforms.index)
        features["is_NNC"] = (isoforms["structural_category"]
                              == "novel_not_in_catalog").astype(int)
        for name in SUBCATEGORIES:
            features[f"sub_{name}"] = (isoforms["subcategory"]
                                       == name).astype(int)
        booster = lgb.Booster(model_str=model["booster"])
        raw = booster.predict(features[model["features"]])
        probability = np.minimum(1.0, raw / model["c"])
    result = probability

    return result
