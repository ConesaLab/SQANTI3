"""
pantranscriptome_build.py

Build the pantranscriptome catalog that SQANTI3 QC reads with
--pantranscriptome. The input is a Parquet table with one row per intron
chain seen in public long-read samples (columns ujc, n_samples,
n_studies and n_platforms; ujc as in intron_chain_key()). The catalog
keeps a 44-bit hash of each chain instead of the chain itself, in two
sorted groups: chains seen in a single sample, whose counts are 1, and
the rest, stored with their counts. Hashes and counts are compressed
with pack_array(). The catalog also carries the pantranscriptome model
that scores the NIC/NNC chains (src/pantranscriptome_model.py): the
LightGBM model text, the donor and acceptor PWMs, the variables in
order, the PU constant c and the probability threshold, read from the
directory written by the training of the model (model.txt,
pwm_donor.npy, pwm_acceptor.npy, meta.json).

Version: 1.2
Author: Adrian Berenguer-Agustí (I2SysBio-CSIC)
Date: 2026-11-08
"""

import argparse
import json
import os
import sys

import numpy as np
import pyarrow.parquet as pq

# Run as a script from any directory: put the SQANTI3 root on the path.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", ".."))
from src.config import PANTRANSCRIPTOME_FORMAT  # noqa: E402
from src.config import PANTRANSCRIPTOME_HASH_BITS  # noqa: E402
from src.utils import hash_intron_chains  # noqa: E402
from src.utils import pack_array  # noqa: E402

# Columns of the input table.
COLUMNS = ["ujc", "n_samples", "n_studies", "n_platforms"]
# Rows read from the input table at a time.
BATCH_SIZE = 1_000_000
# Files of the pantranscriptome model directory.
MODEL_FILES = ("model.txt", "pwm_donor.npy", "pwm_acceptor.npy",
               "meta.json")
# Largest count that fits in the catalog (uint16 samples and studies).
MAX_COUNT = np.iinfo(np.uint16).max


def read_evidence(parquet_path: str) -> dict[str, np.ndarray]:
    """Hash the intron chains of an evidence table.

    Args:
        parquet_path: Parquet table with the columns in COLUMNS.

    Returns:
        Arrays "hash" (uint64), "samples", "studies" and "platforms"
        (int64), one value per row of the table.

    Raises:
        FileNotFoundError: If parquet_path does not exist.
        ValueError: If a column is missing or a count is out of range.
    """
    if not os.path.isfile(parquet_path):
        raise FileNotFoundError(parquet_path)
    table = pq.ParquetFile(parquet_path)
    missing = set(COLUMNS) - set(table.schema_arrow.names)
    if missing:
        raise ValueError(f"{parquet_path} lacks the columns "
                         f"{sorted(missing)}.")
    parts = {"hash": [], "samples": [], "studies": [], "platforms": []}
    for batch in table.iter_batches(batch_size=BATCH_SIZE, columns=COLUMNS):
        parts["hash"].append(hash_intron_chains(
            batch.column(0).to_pylist(), PANTRANSCRIPTOME_HASH_BITS))
        for i, name in enumerate(("samples", "studies", "platforms"), 1):
            parts[name].append(batch.column(i).to_numpy().astype(np.int64))
    result = {name: np.concatenate(values) for name, values in parts.items()}
    counts = np.column_stack([result["samples"], result["studies"],
                              result["platforms"]])
    if counts.min() < 1 or counts.max() > MAX_COUNT:
        raise ValueError(f"Counts must be between 1 and {MAX_COUNT}.")
    if result["platforms"].max() > np.iinfo(np.uint8).max:
        raise ValueError("Too many platforms for the catalog format.")

    return result


def drop_collisions(evidence: dict[str, np.ndarray]
                    ) -> dict[str, np.ndarray]:
    """Keep one row per hash, sorted by hash.

    Two chains with the same hash are indistinguishable in the catalog;
    the one seen in more samples is kept.

    Args:
        evidence: Output of read_evidence().

    Returns:
        The same arrays, sorted by hash, without repeated hashes.
    """
    order = np.lexsort((-evidence["samples"], evidence["hash"]))
    hashes = evidence["hash"][order]
    first = np.ones(len(hashes), dtype=bool)
    first[1:] = hashes[1:] != hashes[:-1]
    result = {name: values[order][first] for name, values in evidence.items()}
    print(f"[pantranscriptome] {len(hashes)} chains, "
          f"{len(hashes) - first.sum()} hash collisions dropped",
          file=sys.stderr)

    return result


def pack_hashes(hashes: np.ndarray) -> np.ndarray:
    """Pack sorted hashes as differences from the previous one.

    Args:
        hashes: Sorted uint64 hashes.

    Returns:
        pack_array() output of the differences (the first one from 0).
    """
    result = pack_array(np.diff(hashes, prepend=np.uint64(0)))

    return result


def read_model(model_dir: str) -> dict[str, np.ndarray]:
    """Read the pantranscriptome model as catalog arrays.

    Args:
        model_dir: Directory with model.txt, pwm_donor.npy,
            pwm_acceptor.npy and meta.json.

    Returns:
        Arrays model_text (packed model text), model_pwm_donor,
        model_pwm_acceptor, model_features, model_c and
        model_threshold.

    Raises:
        FileNotFoundError: If a file of the model is missing.
    """
    paths = {name: os.path.join(model_dir, name)
             for name in MODEL_FILES}
    missing = [path for path in paths.values() if not os.path.isfile(path)]
    if missing:
        raise FileNotFoundError(
            f"Missing pantranscriptome model files: {missing}")
    with open(paths["model.txt"], "rb") as handle:
        text = handle.read()
    with open(paths["meta.json"]) as handle:
        meta = json.load(handle)
    result = {
        "model_text": pack_array(np.frombuffer(text, dtype=np.uint8)),
        "model_pwm_donor": np.load(paths["pwm_donor.npy"]),
        "model_pwm_acceptor": np.load(paths["pwm_acceptor.npy"]),
        "model_features": np.array(meta["features"]),
        "model_c": np.array(float(meta["c"])),
        "model_threshold": np.array(float(meta["threshold"])),
    }

    return result


def build_catalog(parquet_path: str, model_dir: str, out_path: str,
                  assembly: str) -> None:
    """Write the pantranscriptome catalog of an evidence table.

    Args:
        parquet_path: Parquet table with the columns in COLUMNS.
        model_dir: Directory of the pantranscriptome model (see read_model()).
        out_path: Output .npz file.
        assembly: Genome assembly of the coordinates, e.g. GRCh38.

    Raises:
        FileNotFoundError: If parquet_path or a model file does not
            exist.
        ValueError: If the table is not a valid evidence table.
    """
    model = read_model(model_dir)
    evidence = drop_collisions(read_evidence(parquet_path))
    single = ((evidence["samples"] == 1) & (evidence["studies"] == 1)
              & (evidence["platforms"] == 1))
    multi = ~single
    np.savez(out_path,
             format=np.array(PANTRANSCRIPTOME_FORMAT),
             hash_bits=np.array(PANTRANSCRIPTOME_HASH_BITS),
             assembly=np.array(assembly),
             single=pack_hashes(evidence["hash"][single]),
             multi=pack_hashes(evidence["hash"][multi]),
             multi_samples=pack_array(
                 evidence["samples"][multi].astype(np.uint16)),
             multi_studies=pack_array(
                 evidence["studies"][multi].astype(np.uint16)),
             multi_platforms=pack_array(
                 evidence["platforms"][multi].astype(np.uint8)),
             **model)
    print(f"[pantranscriptome] {single.sum()} single-sample and "
          f"{multi.sum()} multi-sample chains -> {out_path}",
          file=sys.stderr)


def main() -> None:
    """Parse the command line and build the catalog."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    parser.add_argument("evidence", help="Parquet table with the columns "
                        + ", ".join(COLUMNS))
    parser.add_argument("model_dir", help="Directory of the pantranscriptome "
                        "model (" + ", ".join(MODEL_FILES) + ")")
    parser.add_argument("output", help="Output catalog (.npz)")
    parser.add_argument("--assembly", default="GRCh38",
                        help="Genome assembly of the coordinates "
                        "(default: %(default)s)")
    args = parser.parse_args()
    build_catalog(args.evidence, args.model_dir, args.output,
                  args.assembly)


if __name__ == "__main__":
    main()
