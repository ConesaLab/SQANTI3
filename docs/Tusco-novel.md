# TUSCO‑novel (Novel Isoform Stress Test)

## Overview
TUSCO‑novel benchmarks a pipeline’s ability to detect truly novel isoforms by intentionally “lying” to the pipeline: the reference GTF is altered only for curated single‑isoform TUSCO genes so that the real expressed isoform appears novel relative to the supplied annotation.

Note: The curated TUSCO gene list is derived from the TSV panel (`src/utilities/report_qc/tusco_<species>.tsv`). TUSCO itself uses TSV panels only; no TUSCO GTF is required.

---

## Concept
- Alters single‑isoform, multi‑exon TUSCO genes to hide one real internal junction and replace it with a plausible synthetic junction (canonical GT–AG; reasonable intron length; still multi‑exon).
- Tools that depend on the given annotation tend to degrade; reference‑agnostic discovery plus robust filtering (e.g., Iso‑Seq + SQANTI3 ML) better controls false positives.

Design constraints (for controlled difficulty and fairness):
- Modify only internal splice junctions; preserve TSS/TTS and multi‑exon structure.
- Use canonical motifs (GT–AG) and intron lengths within empirical bounds for the species/tissue.
- Leave all non‑TUSCO genes unchanged to localize the perturbation.

---

## How It Works
1. Inputs: native reference GTF, genome FASTA, curated TUSCO single‑isoform gene list (human/mouse).
2. For each multi‑exon TUSCO gene:
   - Remove one true internal splice junction in the transcript.
   - Insert a plausible synthetic junction (canonical; reasonable distance; preserves multi‑exon structure).
3. Use this altered GTF as the only “reference” for reconstruction and for downstream classification/evaluation.

Simulator code: https://github.com/TianYuan-Liu/tusco-paper/blob/main/src/tusco_novel_simulator/tusco_novel_sim.py

---

## When To Use
- Quantify dependence on annotation vs sequencing data.
- Assess novel discovery under under‑annotated or misleading references (species/tissues).

---

## End‑to‑End Workflow
### 1) Prepare inputs
- Genome FASTA: `hg38.fa` or `mm10.fa`.
- Native reference GTF: e.g., GENCODE.
- TUSCO gene list: derive from `src/utilities/report_qc/tusco_<species>.tsv` (human or mouse).
  - Extract the TUSCO gene identifiers to a text file, one per line (e.g., `tusco_genes.txt`).
- Your reads/alignments as required by each pipeline (BAMs or CCS/FASTQ for Iso‑Seq).

### 2) Build the TUSCO‑novel GTF
Use the simulator to modify only TUSCO genes:
```bash
python tusco_novel_sim.py \
  --refGTF native.gtf \
  --genome hg38.fa \
  --tusco-list tusco_genes.txt \
  --out tusco_novel.gtf \
  --seed 42
```
Sanity checks:
- Modify only TUSCO single‑isoform genes and only internal junctions.
- Ensure synthetic junctions are canonical and yield valid multi‑exon transcripts.
- Log which junction was modified per gene for reproducibility.

Provenance and determinism:
- Record the commit of the simulator and configuration used.
- Fix the RNG seed (e.g., `--seed 42`) and preserve the simulator log.

### 3) Run reconstruction with the TUSCO‑novel reference
- StringTie2:
  ```bash
  stringtie aligned.bam -G tusco_novel.gtf -o stringtie.gtf -L
  ```
- FLAIR (guide with altered annotation):
  ```bash
  flair correct -q reads.fastq -g hg38.fa -f tusco_novel.gtf -o flair_correct
  flair collapse -g hg38.fa -r reads.fastq -q flair_correct.bed -f tusco_novel.gtf -o flair_collapse
  ```
- Bambu (R): provide `tusco_novel.gtf` as the annotation object.
- Iso‑Seq + SQANTI3 ML: run Iso‑Seq reference‑free; use `tusco_novel.gtf` only for SQANTI3 classification.

### 4) Evaluate with SQANTI3 (TUSCO metrics)
- TUSCO‑novel classification:
  ```bash
  python sqanti3_qc.py \
    --isoforms <tool_output.gtf> \
    --refGTF tusco_novel.gtf \
    --refFasta hg38.fa \
    --tusco human|mouse \
    --report html -o <prefix_novel> -d <outdir>
  ```
- Baseline (native reference):
  ```bash
  python sqanti3_qc.py \
    --isoforms <tool_output.gtf> \
    --refGTF native.gtf \
    --refFasta hg38.fa \
    --tusco human|mouse \
    --report html -o <prefix_native> -d <outdir>
  ```
Outputs are identical to the standard TUSCO report. See [TUSCO Quick Start](Tusco-quick-start.md) for report paths, metric definitions, and interpretation guidance.

---

## Interpreting Results
- Expect larger drops under TUSCO‑novel for reference‑guided tools (e.g., StringTie2, Bambu) if they rely on annotation for novel splice discovery.
- Reference‑free discovery with stringent filtering (Iso‑Seq + SQANTI3 ML) typically:
  - Minimizes false positives (FDR), especially when junction support is enforced.
  - May show lower sensitivity/precision for exact TSS/TTS due to read length and end adjustments.
- Compare native vs TUSCO‑novel: the gap indicates reliance on annotation vs data.

Report context:
- TUSCO‑novel uses the same TSV panel and metrics as TUSCO. Interpreting Sn, nrPre, rPre, 1−FDR, PDR, and 1/red follows the same guidance as in the Quick Start.

---

## Practical Tips
- Restrict modification to multi‑exon TUSCO genes with well‑supported junctions; leave the rest of the annotation unchanged.
- Use a fixed RNG seed and write a per‑gene change log.
- Verify synthetic junctions against the genome (canonical motifs, reasonable intron sizes).
- Keep pipeline parameters identical between native and TUSCO‑novel runs for fair comparison.

Limitations:
- The stress test focuses on splice‑junction novelty; it does not simulate alternative TSS/TTS or structural variants.
- Results can depend on the choice of TUSCO panel (species/tissue) and the simulator constraints; report both explicitly.

---

## See Also
- [TUSCO Quick Start](Tusco-quick-start.md)
- [SQANTI3 Quality Control](Running-SQANTI3-Quality-Control.md)

## Reproducibility Checklist
- SQANTI3 version/commit and command lines for both native and TUSCO‑novel runs.
- Provenance of the TUSCO panel TSV (filename, species, checksum).
- Simulator commit, configuration, RNG seed, and logs.
- All generated logs (including `tusco_report.log`) and HTML reports archived.
