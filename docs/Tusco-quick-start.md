## Overview

**TUSCO** (Transcriptome Universal Single-isoform COntrol) is a curated internal reference set of genes lacking alternative isoforms, designed to benchmark long-read transcriptome sequencing quality without external spike-in controls.

Unlike BUSCO (which can misinterpret alternative splicing as gene duplications) or spike-ins like SIRVs/ERCCs (which oversimplify real-sample complexity and neglect RNA degradation artifacts), TUSCO uses endogenous single-isoform genes to evaluate:

- **Precision**: Identifying transcripts that deviate from reference annotations
- **Sensitivity**: Verifying detection completeness of known transcripts

This module is integrated into SQANTI3 QC and generates an interactive HTML report with benchmarking metrics and IGV-style genome visualization plots.

---

## Example Dataset

The [`data/tusco/`](https://github.com/ConesaLab/SQANTI3/tree/master/data/tusco) directory contains a minimal dataset for testing SQANTI3 TUSCO benchmarking (~8.8 MB total):

| File | Size | Description |
|------|------|-------------|
| `tusco_genome.fa` | 4.7 MB | Genomic regions around 46 TUSCO genes |
| `tusco_genome.fa.fai` | 2 KB | FASTA index |
| `tusco_annotation.gtf` | 3.7 MB | GENCODE v49 annotation subset (8,688 features) |
| `tusco_input.gtf` | 13 KB | Example input transcripts (112 entries, 40 genes) |

---

## Prerequisites

### Create Conda Environment

Navigate to the SQANTI3 directory and create the conda environment:

```bash
cd /path/to/SQANTI3
conda env create -f SQANTI3.conda_env.yml
```

This installs ~70+ packages including Python 3.11, R 4.3+, bioinformatics tools (minimap2, samtools, bedtools), and R packages (ggplot2, plotly, Gviz, rmarkdown).

**Tip:** For faster installation, use mamba instead of conda:
```bash
mamba env create -f SQANTI3.conda_env.yml
```

### Activate Environment

```bash
conda activate sqanti3
```

---

## Basic Usage

Run SQANTI3 with TUSCO benchmarking enabled:

```bash
python sqanti3_qc.py \
    --isoforms data/tusco/tusco_input.gtf \
    --refGTF data/tusco/tusco_annotation.gtf \
    --refFasta data/tusco/tusco_genome.fa \
    --tusco human \
    -d output/tusco_example \
    --skipORF
```

**Parameters:**
- `--isoforms`: Input transcript GTF file to benchmark
- `--refGTF`: Reference annotation GTF
- `--refFasta`: Reference genome FASTA
- `--tusco human`: Enable TUSCO benchmarking (use `human` or `mouse`)
- `-d`: Output directory
- `--skipORF`: Skip ORF prediction (faster for testing)

---

## Outputs

The TUSCO module generates the following output files:

| File | Description |
|------|-------------|
| `<prefix>_TUSCO_report.html` | Interactive HTML benchmarking report with metrics and visualizations |
| `<prefix>_TUSCO_results.tsv` | Transcript categorization (transcript_id, associated_gene, structural_category, subcategory, TUSCO_category) |
| `igv_plots/` | Directory with IGV-style genome visualization PNG plots (one per gene) |
| `logs/tusco_report.log` | Execution log for troubleshooting |

Standard SQANTI3 outputs are also generated:
- `<prefix>_classification.txt` - Full SQANTI3 classification results
- `<prefix>_corrected.gtf` - Corrected GTF file
- `<prefix>_junctions.txt` - Junction information

To view the report, open `<prefix>_TUSCO_report.html` in a web browser.

---

## Bundled Reference Panels

SQANTI3 includes pre-built TUSCO reference panels:

| Species | File | Genes | Location |
|---------|------|-------|----------|
| Human | `tusco_human.tsv` | 46 | `src/utilities/report_qc/` |
| Mouse | `tusco_mouse.tsv` | 33 | `src/utilities/report_qc/` |

Each TSV file contains columns: Ensembl Gene ID, Ensembl Transcript ID, Gene Symbol, Entrez ID, RefSeq mRNA, RefSeq Protein.

---

## Report Contents and Interpretation

### Benchmarking Metrics

The TUSCO report calculates 7 benchmarking metrics:

| Metric | Description |
|--------|-------------|
| **Sensitivity** | Proportion of reference transcripts correctly detected (TP / (TP + FN)) |
| **Non-redundant Precision** | Proportion of unique predicted transcripts that are correct (TP / (TP + FP)) |
| **Redundant Precision** | Precision accounting for redundant predictions |
| **Positive Detection Rate** | Rate of true positive detections among all predictions |
| **False Discovery Rate** | Proportion of predictions that are false positives (FP / (TP + FP)) |
| **False Detection Rate** | Rate of reference transcripts not detected (FN / (TP + FN)) |
| **Redundancy** | Ratio of total predictions to unique predictions |

### Classification Categories

Transcripts are classified into four TUSCO categories:

| Category | Definition |
|----------|------------|
| **TP** (True Positive) | Exact structural match to a TUSCO reference transcript (FSM - Full Splice Match) |
| **PTP** (Partial True Positive) | Partial match to reference (ISM - Incomplete Splice Match, or NIC/NNC with shared junctions) |
| **FN** (False Negative) | TUSCO reference transcript not detected in the input |
| **FP** (False Positive) | Predicted transcript that does not match any TUSCO reference |

---

## Troubleshooting

### Conda Solver Issues
If conda fails to solve the environment:
1. Use mamba instead of conda (faster solver)
2. Remove specific version constraints if packages are unavailable for your platform

### Apple Silicon (ARM64) Notes
Some packages may have limited ARM64 support. The core TUSCO functionality works on Apple Silicon, but you may need to:
- Install packages without strict version pins
- Skip optional dependencies like `parasail` if they fail to build

### Missing R Packages
If R packages fail to load, verify they are installed in the conda environment:
```bash
Rscript -e "library(Gviz); library(ggplot2); library(plotly); library(rmarkdown)"
```

### Empty or Missing IGV Plots
If IGV plots are not generated:
- Check `logs/tusco_report.log` for errors
- Ensure the reference genome FASTA matches the annotation coordinates
- Verify Gviz R package is properly installed

---

## Source Data

- **Genome:** GRCh38.p14 (extracted regions)
- **Annotation:** GENCODE v49
- **Input:** WTC11 PacBio cDNA transcripts (subset)

### Coordinate System

Files use region-based chromosome names (e.g., `chr1:1182237-1285041`) to match the extracted genomic regions.

---

## See Also

- [SQANTI3 Documentation](https://github.com/ConesaLab/SQANTI3)
- [TUSCO Preprint](https://doi.org/10.1101/2025.08.23.671926) - Liu et al., bioRxiv 2025
