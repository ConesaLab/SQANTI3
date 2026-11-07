"""
test_pantranscriptome.py

Unit tests of the pantranscriptome annotation of SQANTI3 QC: catalog
packing, intron chain keys, catalog building and lookup, the pubLR
variables and the columns added to NIC/NNC.

Version: 1.1
Author: Adrian Berenguer-Agustí (I2SysBio-CSIC)
Date: 2026-11-07
"""

import json
import os
import sys

import lightgbm as lgb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.classification_steps import intron_chain_key  # noqa: E402
from src.qc_classes import (  # noqa: E402
    Pantranscriptome, genePredRecord, myQueryTranscripts
)
from src.publr_model import (  # noqa: E402
    SUBCATEGORIES, annotated_introns, chain_features, genome_codes
)
from src.qc_computations import pantranscriptome_annotation  # noqa: E402
from src.utilities.pantranscriptome_build import build_catalog  # noqa: E402
from src.utils import pack_array, unpack_array  # noqa: E402

# Chains of the test catalog: samples, studies and platforms.
CHAINS = {
    "chr1+150-200_300-400": (10, 4, 2),   # recurrent
    "chr1-150-200": (1, 1, 1),            # seen in one sample
    "chr2+10-20_30-40_50-60": (3, 2, 1),  # seen, not recurrent
    "chr3+5-9": (2, 3, 1),                # 3 studies, one platform
}
# Variables of the test pubLR model, in the order of the real one.
FEATURES = (["n_junctions", "intron_log_min", "intron_log_max",
             "intron_log_mean", "internal_exon_min", "n_noncanonical",
             "donor_score_min", "donor_score_mean", "acceptor_score_min",
             "acceptor_score_mean", "n_novel_junctions", "n_novel_donors",
             "n_novel_acceptors", "n_shifted_sites",
             "novel_site_log_distance_min", "repeat_max",
             "novel_repeat_max", "frac_novel_junctions", "is_NNC"]
            + [f"sub_{name}" for name in SUBCATEGORIES])
# Length of the chromosomes of the test genome.
CHROM_LENGTH = 1000


@pytest.fixture
def model_dir(tmp_path) -> str:
    """Train a small pubLR-like model on random data and save it."""
    rng = np.random.default_rng(1)
    data = lgb.Dataset(rng.normal(size=(200, len(FEATURES))),
                       label=rng.integers(0, 2, 200),
                       feature_name=FEATURES)
    booster = lgb.train({"objective": "binary", "verbose": -1,
                         "min_data_in_leaf": 5}, data, num_boost_round=5)
    result = str(tmp_path / "model")
    os.makedirs(result)
    booster.save_model(os.path.join(result, "model.txt"))
    np.save(os.path.join(result, "pwm_donor.npy"), rng.normal(size=(9, 4)))
    np.save(os.path.join(result, "pwm_acceptor.npy"),
            rng.normal(size=(23, 4)))
    with open(os.path.join(result, "meta.json"), "w") as handle:
        json.dump({"features": FEATURES, "c": 0.5, "threshold": 0.5},
                  handle)

    return result


@pytest.fixture
def reference(tmp_path) -> tuple[str, str]:
    """Random genome FASTA and a GTF with the intron 151-299 of chr1."""
    rng = np.random.default_rng(2)
    fasta = tmp_path / "genome.fa"
    records = []
    for chrom in ("chr1", "chr2", "chr3", "chr9"):
        sequence = "".join(rng.choice(list("ACGT"), CHROM_LENGTH))
        records.append(f">{chrom}\n{sequence}\n")
    fasta.write_text("".join(records))
    gtf = tmp_path / "reference.gtf"
    attr = 'gene_id "G1"; transcript_id "T1";'
    gtf.write_text(f"chr1\ttest\texon\t100\t150\t.\t+\t.\t{attr}\n"
                   f"chr1\ttest\texon\t300\t400\t.\t+\t.\t{attr}\n")
    result = (str(fasta), str(gtf))

    return result


@pytest.fixture
def catalog_path(tmp_path, model_dir: str) -> str:
    """Build a catalog of CHAINS and return its path."""
    keys = list(CHAINS)
    counts = np.array([CHAINS[key] for key in keys])
    table = pa.table({"ujc": keys, "n_samples": counts[:, 0],
                      "n_studies": counts[:, 1],
                      "n_platforms": counts[:, 2]})
    parquet_path = str(tmp_path / "evidence.parquet")
    pq.write_table(table, parquet_path)
    result = str(tmp_path / "catalog.npz")
    build_catalog(parquet_path, model_dir, result, "GRCh38")

    return result


@pytest.mark.parametrize("dtype", ["uint8", "uint16", "uint64"])
def test_pack_round_trip(dtype: str) -> None:
    values = np.array([0, 1, 7, 200, 255], dtype=dtype)
    assert np.array_equal(unpack_array(pack_array(values), dtype), values)
    empty = np.array([], dtype=dtype)
    assert len(unpack_array(pack_array(empty), dtype)) == 0


def test_intron_chain_key() -> None:
    # genePred starts are 0-based: exons 100-150, 201-300 and 401-500.
    record = genePredRecord(id="t1", chrom="chr1", strand="+", txStart=99,
                            txEnd=500, cdsStart=99, cdsEnd=99, exonCount=3,
                            exonStarts=[99, 200, 400],
                            exonEnds=[150, 300, 500])
    assert intron_chain_key(record) == "chr1+150-201_300-401"


def test_catalog_lookup(catalog_path: str) -> None:
    catalog = Pantranscriptome(catalog_path)
    assert catalog.assembly == "GRCh38"
    keys = [*CHAINS, "chr9+1-2"]
    expected = np.array([*CHAINS.values(), (0, 0, 0)])
    assert np.array_equal(catalog.lookup(keys), expected)
    assert catalog.model["features"] == FEATURES
    assert catalog.model["c"] == 0.5
    lgb.Booster(model_str=catalog.model["booster"])


def test_chain_features(reference: tuple[str, str]) -> None:
    genome = genome_codes(None, reference[0], {"chr1"})
    assert len(genome["chr1"]) == CHROM_LENGTH
    introns = annotated_introns(reference[1])
    assert introns[["start", "end"]].values.tolist() == [[151, 299]]
    # One annotated intron, and the same chain with the acceptor moved
    # 5 nt and a second, novel intron.
    chains = pd.Series(["chr1+150-300", "chr1+150-305_400-500"])
    pwms = (np.zeros((9, 4)), np.zeros((23, 4)))
    features = chain_features(chains, genome, introns, pwms)
    known = features.loc["chr1+150-300"]
    assert known["n_junctions"] == 1 and known["n_novel_junctions"] == 0
    assert np.isnan(known["novel_site_log_distance_min"])
    shifted = features.loc["chr1+150-305_400-500"]
    assert shifted["n_novel_junctions"] == 2
    assert shifted["n_novel_donors"] == 1
    assert shifted["n_shifted_sites"] == 1
    # Internal exon: bases 305-400.
    assert shifted["internal_exon_min"] == 96


def test_not_a_catalog(tmp_path) -> None:
    path = tmp_path / "not_a_catalog.npz"
    path.write_text("ujc\tn_samples\n")
    with pytest.raises(ValueError):
        Pantranscriptome(str(path))


def test_annotation(catalog_path: str, reference: tuple[str, str]
                    ) -> None:
    isoforms = {
        "nic": myQueryTranscripts(isoform="nic",
                                  structural_category="novel_in_catalog",
                                  ujc="chr1+150-200_300-400"),
        "nnc": myQueryTranscripts(isoform="nnc",
                                  structural_category="novel_not_in_catalog",
                                  ujc="chr3+5-9"),
        "new": myQueryTranscripts(isoform="new",
                                  structural_category="novel_not_in_catalog",
                                  ujc="chr9+100-200"),
        "fsm": myQueryTranscripts(isoform="fsm",
                                  structural_category="full-splice_match",
                                  genes=["G1"], transcripts=["T1"]),
    }
    for name in ("nic", "nnc", "new"):
        isoforms[name].chrom = isoforms[name].ujc.split("+")[0]
    pantranscriptome_annotation(isoforms, catalog_path, reference[1],
                                reference[0])
    labels = {name: iso.pantranscriptome for name, iso in isoforms.items()}
    assert labels == {"nic": "recurrent", "nnc": "seen", "new": "unseen",
                      "fsm": None}
    assert (isoforms["nic"].pub_samples, isoforms["nic"].pub_studies,
            isoforms["nic"].pub_platforms) == (10, 4, 2)
    assert isoforms["new"].pub_samples == 0
    row = isoforms["fsm"].as_dict()
    assert row["pantranscriptome"] == "NA" and row["pub_samples"] == "NA"
    assert "ujc" not in row and row["pubLR_prob"] == "NA"
    for name in ("nic", "nnc", "new"):
        assert 0.0 <= isoforms[name].pubLR_prob <= 1.0


def test_annotation_without_catalog() -> None:
    iso = myQueryTranscripts(isoform="nic",
                             structural_category="novel_in_catalog",
                             ujc="chr1+150-200_300-400")
    pantranscriptome_annotation({"nic": iso}, None, "unused.gtf",
                                "unused.fa")
    assert iso.pantranscriptome is None and iso.pub_samples is None
    assert iso.pubLR_prob is None
