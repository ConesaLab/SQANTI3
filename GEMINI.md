# GEMINI.MD: AI Collaboration Guide

This document provides essential context for AI models interacting with the SQANTI3 project. Adhering to these guidelines will ensure consistency and maintain code quality.

## 1. Project Overview & Purpose

* **Primary Goal:** SQANTI3 is a tool for the quality control, characterization, and structural analysis of long-read transcriptomes. It merges features from SQANTI and SQANTI2 with new additions to provide comprehensive isoform classification and quantification.
* **Business Domain:** Bioinformatics, specifically Long-Read Transcriptomics (Oxford Nanopore, PacBio Iso-Seq). It is part of the Functional IsoTranscriptomics (FIT) framework.

## 2. Core Technologies & Stack

* **Languages:** Python 3.11.x, R >= 4.3.3.
* **Frameworks & Runtimes:** CLI-based Python scripts, R for report generation (RMarkdown).
* **Databases:** None (File-based: GTF/GFF3, BAM, SAM, VCF, BED, FASTQ, FASTA).
* **Key Libraries/Dependencies:** 
    * **Python:** `pysam`, `pandas`, `numpy`, `scikit-learn`, `biopython`, `bx-python`, `pyyaml`.
    * **R:** `ggplot2`, `dplyr`, `tidyr`, `rmarkdown`, `plotly`, `caret`.
    * **Bioinformatics Tools:** `minimap2`, `samtools`, `bedtools`, `gffread`, `STAR`, `kallisto`, `TransDecoder`.
* **Package Manager(s):** `conda` (Mamba), `pixi`.

## 3. Architectural Patterns

* **Overall Architecture:** Modular CLI application. The project consists of several main entry-point scripts that utilize a shared core library in `src/`.
* **Directory Structure Philosophy:**
    * `/src`: Core logic, divided into feature-specific modules (classification, qc, filter, rescue, reads).
    * `/src/utilities`: Helper scripts, bundled third-party tools (e.g., Cupcake), and RMarkdown templates for reports.
    * `/test`: Comprehensive test suite including unit, integration, and functional tests.
    * `/docs`: Project documentation and wiki source files.
    * Root directory: Contains the main CLI entry points (`sqanti3_qc.py`, `sqanti3_filter.py`, etc.) and configuration files.

## 4. Coding Conventions & Style Guide

* **Formatting:**
    * **Python:** Adheres to PEP 8. Formatting is checked via `flake8` (configured in `setup.cfg`). Note: `max-line-length` is set to 500.
    * **R:** Follows Tidyverse style preferences. Checked via `.lintr`.
* **Naming Conventions:**
    * `variables`, `functions`: `snake_case` (e.g., `isoform_classification_pipeline`).
    * `classes`: `PascalCase` (e.g., `SAMReader`).
    * `files`: `snake_case` (e.g., `qc_pipeline.py`).
* **API Design:** Primarily a CLI-based tool. Arguments are managed using `argparse`, with centralized validation logic in `src/argparse_utils.py`.
* **Error Handling:** Uses standard Python exceptions and centralized logging via `src/module_logging.py` and `src/logging_config.py`.

## 5. Key Files & Entrypoints

* **Main Entrypoint(s):**
    * `sqanti3_qc.py`: Quality Control pipeline.
    * `sqanti3_filter.py`: Isoform filtering.
    * `sqanti3_rescue.py`: Isoform rescue module.
    * `sqanti3_reads.py`: SQANTI-reads analysis.
    * `sqanti3`: Wrapper script for the above tools.
* **Configuration:** 
    * `sqanti3_config.yaml`: Default configuration parameters.
    * `src/config.py`: Hardcoded configuration and versioning.
* **CI/CD Pipeline:** GitHub Actions workflows in `.github/workflows/` (e.g., `run-tests.yml`, `build-test-conda.yml`).

## 6. Development & Testing Workflow

* **Local Development Environment:** 
    * Managed via Conda/Mamba using `SQANTI3.conda_env.yml` or Pixi using `pixi.toml`.
    * Docker support via `Dockerfile` and `build_docker.sh`.
* **Testing:** 
    * Run Python tests using `pytest`.
    * Tests are organized in `test/unit/`, `test/integration/`, and `test/functional/`.
    * `test/utils.py` contains helper functions for testing.
* **CI/CD Process:** Automated tests run on PRs and commits to `master` via GitHub Actions. Docker images are automatically built and pushed upon release.

## 7. Specific Instructions for AI Collaboration

* **Contribution Guidelines:** (From `CONTRIBUTING.md`)
    * Fork the repo and create branches from `master`.
    * Add tests for new features or bug fixes.
    * Ensure all tests pass and code lints correctly.
    * Submissions are under the GPL-3.0 License.
* **Security:** Handle genomic data securely; do not hardcode paths or assume specific local file structures. Use the provided validation helpers in `src/argparse_utils.py`.
* **Dependencies:** New dependencies should be added to `pyproject.toml`, `pixi.toml`, and `SQANTI3.conda_env.yml`.
* **Commit Messages:** Follow existing patterns in the commit history. Messages are generally concise and descriptive (e.g., `Fixed qc_pipeline full_length_quantification function`, `Refactor parallel.py to remove unused import`).
* **Performance:** Prioritize memory efficiency and parallelization, as SQANTI3 is often used with large genomic datasets. Use `src/parallel.py` for implementing parallel processing.

## 8. Documentation

* **System & Framework:** Documentation is managed using **MkDocs** with the **Material theme** (configured in `mkdocs.yml`). It is also mirrored/maintained as a GitHub Wiki.
* **Format:** All documentation source files are written in **Markdown** (`.md`) and located in the `/docs` directory.
* **Documentation Style:**
    * **Technical & Comprehensive:** Provides in-depth guides for each CLI module, including detailed argument explanations, input requirements, and algorithmic logic.
    * **Visual & Diagrammatic:** Frequently uses diagrams and flowcharts to explain complex bioinformatic workflows (e.g., rescue strategy, isoform classification).
    * **Structured Content:** Most pages include a "Table of Contents" with anchor links for easy navigation. Use of `⚠️WARNING` and `> Note` callouts is standard for highlighting critical information.
    * **Domain-Specific Guidance:** Includes theoretical sections explaining biological concepts (e.g., Isoform categories, Intra-priming, RT-switching) to provide context for the tool's output.
    * **Educational Resources:** Features both step-by-step written tutorials and video tutorial links.
* **Maintenance:** Navigation structure is centrally defined in `mkdocs.yml`. When adding new features or modules, corresponding documentation must be added to `/docs` and updated in the navigation menu.
