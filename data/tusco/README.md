# TUSCO Example Dataset

## Overview
Minimal dataset for testing SQANTI3 TUSCO benchmarking (~8.8 MB total).

## Files

| File | Size | Description |
|------|------|-------------|
| `tusco_genome.fa` | 4.7 MB | Genomic regions around 46 TUSCO genes (±50kb flanking) |
| `tusco_genome.fa.fai` | 2 KB | FASTA index |
| `tusco_annotation.gtf` | 3.7 MB | GENCODE v49 annotation subset (8,688 features) |
| `tusco_input.gtf` | 13 KB | Example input transcripts (112 entries, 40 genes) |

## Source Data
- Genome: GRCh38.p14 (extracted regions)
- Annotation: GENCODE v49
- Input: WTC11 PacBio cDNA transcripts (subset)

## Coordinate System
Files use region-based chromosome names (e.g., `chr1:1182237-1285041`)
to match the extracted genomic regions.

## Usage

```bash
python sqanti3_qc.py \
    --isoforms data/tusco/tusco_input.gtf \
    --refGTF data/tusco/tusco_annotation.gtf \
    --refFasta data/tusco/tusco_genome.fa \
    --tusco human \
    -d output_dir
```

## TUSCO Genes
Contains all 46 human TUSCO reference genes from `src/utilities/report_qc/tusco_human.tsv`.
