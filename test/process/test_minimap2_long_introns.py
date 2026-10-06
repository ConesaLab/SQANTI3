"""
Process test for the minimap2 command used in mapping mode (issue #216).
Reference transcripts with introns longer than minimap2's default limit (200 kb) must align as a
single primary alignment, not be split into supplementary alignments.
"""
import os
import re
import shutil
import subprocess
from collections import defaultdict

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from src.commands import MINIMAP2_CMD

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
GENOME = os.path.join(main_path, "test/test_data/genome/genome_test.fasta")
REF_GTF = os.path.join(main_path, "test/test_data/reference/test_reference.gtf")

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not installed")

@pytest.fixture
def long_intron_fasta(tmp_path):
    """Spliced sequences of the reference transcripts with an intron longer than 200 kb."""
    genome = {r.id: r.seq for r in SeqIO.parse(GENOME, "fasta")}
    exons, strands = defaultdict(list), {}
    for line in open(REF_GTF):
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon" or f[0] not in genome:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons[(tid, f[0])].append((int(f[3]), int(f[4])))
        strands[tid] = f[6]

    fasta = tmp_path / "long_introns.fa"
    n = 0
    with open(fasta, "w") as out:
        for (tid, chrom), ex in exons.items():
            ex.sort()
            if max((b[0] - a[1] - 1 for a, b in zip(ex, ex[1:])), default=0) <= 200000:
                continue
            seq = Seq("".join(str(genome[chrom][s - 1:e]) for s, e in ex))
            if strands[tid] == "-":
                seq = seq.reverse_complement()
            out.write(f">{tid}\n{seq}\n")
            n += 1
    assert n > 0, "test annotation has no transcripts with introns > 200 kb"
    return str(fasta), n

def test_long_introns_not_split(long_intron_fasta, tmp_path):
    fasta, n_transcripts = long_intron_fasta
    sam = tmp_path / "long_introns.sam"
    cmd = MINIMAP2_CMD.format(cpus=2, g=GENOME, i=fasta, o=sam)
    subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)

    records = defaultdict(list)
    for line in open(sam):
        if not line.startswith("@"):
            f = line.split("\t")
            records[f[0]].append(int(f[1]))

    assert len(records) == n_transcripts
    split = [tid for tid, flags in records.items() if any(flag & 2048 for flag in flags)]
    assert split == [], f"Transcripts split into supplementary alignments: {split}"
