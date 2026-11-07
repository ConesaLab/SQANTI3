"""
test_pantranscriptome.py

Unit tests of the pantranscriptome annotation of SQANTI3 QC: catalog
packing, intron chain keys, catalog building and lookup, and the
columns added to NIC/NNC.

Version: 1.0
Author: Adrian Berenguer-Agustí (I2SysBio-CSIC)
Date: 2026-11-07
"""

import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.classification_steps import intron_chain_key  # noqa: E402
from src.qc_classes import (  # noqa: E402
    Pantranscriptome, genePredRecord, myQueryTranscripts
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


@pytest.fixture
def catalog_path(tmp_path) -> str:
    """Build a catalog of CHAINS and return its path."""
    keys = list(CHAINS)
    counts = np.array([CHAINS[key] for key in keys])
    table = pa.table({"ujc": keys, "n_samples": counts[:, 0],
                      "n_studies": counts[:, 1],
                      "n_platforms": counts[:, 2]})
    parquet_path = str(tmp_path / "evidence.parquet")
    pq.write_table(table, parquet_path)
    result = str(tmp_path / "catalog.npz")
    build_catalog(parquet_path, result, "GRCh38")

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


def test_not_a_catalog(tmp_path) -> None:
    path = tmp_path / "not_a_catalog.npz"
    path.write_text("ujc\tn_samples\n")
    with pytest.raises(ValueError):
        Pantranscriptome(str(path))


def test_annotation(catalog_path: str) -> None:
    isoforms = {
        "nic": myQueryTranscripts(isoform="nic",
                                  structural_category="novel_in_catalog",
                                  ujc="chr1+150-200_300-400"),
        "nnc": myQueryTranscripts(isoform="nnc",
                                  structural_category="novel_not_in_catalog",
                                  ujc="chr3+5-9"),
        "new": myQueryTranscripts(isoform="new",
                                  structural_category="novel_not_in_catalog",
                                  ujc="chr9+1-2"),
        "fsm": myQueryTranscripts(isoform="fsm",
                                  structural_category="full-splice_match",
                                  genes=["G1"], transcripts=["T1"]),
    }
    pantranscriptome_annotation(isoforms, catalog_path)
    labels = {name: iso.pantranscriptome for name, iso in isoforms.items()}
    assert labels == {"nic": "recurrent", "nnc": "seen", "new": "unseen",
                      "fsm": None}
    assert (isoforms["nic"].pub_samples, isoforms["nic"].pub_studies,
            isoforms["nic"].pub_platforms) == (10, 4, 2)
    assert isoforms["new"].pub_samples == 0
    row = isoforms["fsm"].as_dict()
    assert row["pantranscriptome"] == "NA" and row["pub_samples"] == "NA"
    assert "ujc" not in row


def test_annotation_without_catalog() -> None:
    iso = myQueryTranscripts(isoform="nic",
                             structural_category="novel_in_catalog",
                             ujc="chr1+150-200_300-400")
    pantranscriptome_annotation({"nic": iso}, None)
    assert iso.pantranscriptome is None and iso.pub_samples is None
