# Introduction

This page contains a detailed time-aware (DD/MM/YYYY) account of **SQANTI3 releases (>5.0)** and the new features and changes introduced in each of them.

Please be aware that v5.0 represented a major release of the SQANTI3 software. Versions of SQANTI3 >= 5.0 do not have backward compatibility with previous releases and their output (v4.3 and earlier). Users that wish to apply any of the new functionalities in v5.0 to output files from older versions will herefore need to re-run SQANTI3 QC.

See the dedicated site for [installation instructions](Dependencies-and-installation.md).

| [See all SQANTI3 releases](https://github.com/ConesaLab/SQANTI3/tags)



## SQANTI3 v6.1 [LATEST, 20/10/2026]

#### Major changes:
- **Supplementary alignments** ([issue #216](https://github.com/ConesaLab/SQANTI3/issues/216)): minimap2 now allows introns up to 2 Mb (`-G 2000000`, the same limit as GMAP). Before, transcripts with introns longer than 200 kb were split into several alignments. Only the primary alignment of each isoform is classified. Supplementary alignments are listed in `_supplementary_alignments.txt`, and isoforms that could not be aligned in `_unmapped.txt`.
- **SQANTI-reads report**: reworked plots, an interactive HTML report that matches the PDF, new UJC-level, junction, ISM and RTS metrics, upset and saturation plots, and a `summary.json` with the run parameters. SQANTI-reads also runs faster.
- **TSS ratio with short reads**: reads are now counted with pysam through the BAM index instead of `bedtools coverage`. The results are the same, and memory use drops from more than 200 GB to about 1 GB on a 72 GB BAM. Unindexed `--SR_bam` files are indexed once into `<dir>/bam_index`.
- **ML filter and rescue**: numeric columns are no longer turned into text after NA replacement ([PR #633](https://github.com/ConesaLab/SQANTI3/pull/633)), and `ratio_exp` is no longer truncated to an integer ([issue #634](https://github.com/ConesaLab/SQANTI3/issues/634)). ⚠️**Warning**: because of these two fixes, the ML filter can classify some isoforms differently than v6.0.2 on the same data.
- **Prevalence per experimental group** ([PR #636](https://github.com/ConesaLab/SQANTI3/pull/636)): with a multi-sample `--fl_count` file, the new `--counts_design` option assigns the samples to experimental groups, and QC adds one `prevalence_<group>` column per group. `--min_expression` sets the minimum count for a transcript to be expressed in a sample (by default, any count above 0, so fractional counts are kept). The rules filter can use `prevalence` and `prevalence_<group>` like any other column, for example to keep transcripts detected in two samples of any group, or of every group.
- **Rescue requantifies by default** ([PR #636](https://github.com/ConesaLab/SQANTI3/pull/636)): requantification now runs unless `--skip_requant` is given, and it needs the `--counts` file. The `-q/--requant` option was removed. ⚠️**Warning**: scripts and wrapper configuration files that use `-q`, `--requant` or `requant: true` must be updated.
- **Evidence check of rescued reference transcripts** ([PR #636](https://github.com/ConesaLab/SQANTI3/pull/636)): reference transcripts have no long-read counts of their own, so the rules filter applied to the reference transcriptome no longer uses the requisites on `FL`, `FL.<sample>`, `prevalence` or `prevalence_<group>`. Instead, rescue computes the counts that each reference transcript receives from its artifacts and applies those requisites to them before reintroducing it. Reference transcripts that fail are not reintroduced, and their counts go to a long read-defined isoform or to the gene residual. `--skip_evidence_check` keeps the behaviour of previous versions.

#### Minor fixes/enhancements:
- New `min_intron_length` column in the classification (the shortest intron of each transcript, NA for mono-exonic transcripts), which can be used in the rules filter ([PR #636](https://github.com/ConesaLab/SQANTI3/pull/636)).
- After requantification, `prevalence` and `prevalence_<group>` are recomputed in the rescued classification, including the reintroduced reference transcripts ([PR #636](https://github.com/ConesaLab/SQANTI3/pull/636)).
- Rescue stops at the start, before any mapping, when `--counts` is missing, or when the rules need a `--counts_design` file that was not given. QC and rescue also check that every sample of the `--counts_design` file is in the count file before any long step runs.
- Rescue no longer reuses the mapping hits of an earlier run in the same output folder, which crashed rescue when it was run again with a different filter result.
- New `prevalence` column in the classification: the number of samples in which each transcript is detected ([PR #627](https://github.com/ConesaLab/SQANTI3/pull/627)).
- Fractional FL counts are kept, and rescue distributes them proportionally ([issue #624](https://github.com/ConesaLab/SQANTI3/issues/624)).
- With several samples, isoforms missing from one sample's counts get 0 for that sample ([issue #613](https://github.com/ConesaLab/SQANTI3/issues/613)).
- Rescue-by-mapping no longer produces duplicated transcripts ([issue #622](https://github.com/ConesaLab/SQANTI3/issues/622)).
- The ML filter reuses an existing `_randomforest.RData` instead of retraining every run. ML rescue with a classifier trained on multi-sample counts now stops with a clear message instead of crashing.
- GFF/GFF3 reference annotations are no longer rejected at input validation ([PR #635](https://github.com/ConesaLab/SQANTI3/pull/635)).
- QC accepts gzipped FASTA/FASTQ as `--isoforms`.
- FASTQ input can be run in parallel (`-n` > 1), and the renamed FASTA is written to `--dir` instead of the input folder.
- With FASTA/FASTQ input and `-n` > 1, the sequences are aligned once, with all the CPUs given with `-t`, instead of once per chunk. Before, each chunk ran its own aligner and loaded its own copy of the genome index at the same time, so runs with many chunks could run out of memory ([PR #643](https://github.com/ConesaLab/SQANTI3/pull/643)).
- QC stops with an error when two isoform IDs become identical after renaming.
- Fixed `--SR_bam` input paths ([issue #611](https://github.com/ConesaLab/SQANTI3/issues/611)), the `--ratio_TSS_metric 3quartile` crash, `RTS_junction` written as `????` in `_junctions.txt`, and transcript ID parsing in TD2 ORF headers ([PR #626](https://github.com/ConesaLab/SQANTI3/pull/626)).

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v6.1)

_______

## SQANTI3 v6.0.2 (22/07/2026)

#### Minor fixes/enhancements:
- Fixes in Cupcake SAM to GFF3 conversion and short-reads utilities.
- Environment and build configuration updates (`pixi.toml`, `pyproject.toml`).

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v6.0.2)

_______

## SQANTI3 v6.0.1 (24/04/2026)

#### Minor fixes/enhancements:
- Bugfixes in parallel processing and QC computations.
- Added failsafe handling for low-memory environments.
- Updated rules filter helper functions and test suites.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v6.0.1)

_______

## SQANTI3 v6.0 (30/03/2026)

#### Major changes:
- **Rescue Module Overhaul**: Fully revised transcript recovery and reintroduction logic with stricter criteria.
- **CDS Predictions & Classification**: Updated classification table attributes and report format for CDS predictions.
- **TUSCO Integration**: Included TUSCO (Transcriptome Universal Single-isoform COntrol) reference datasets for benchmarking transcriptome reconstruction.
- **Parallel Execution**: Refactored parallel execution framework (`src/parallel.py`).

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v6.0)

_______

## SQANTI3 v5.5 (23/06/2025)

#### Major changes:
- Replaced GeneMarkS-T with **TransDecoder2 (TD2)** for improved ORF predictions.
- Conda environment updated to Python 3.11 to prevent memory freezes on large datasets.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.5)

_______

## SQANTI3 v5.4 (23/04/2025)

#### Major changes:
- Standardized command-line argument naming across Quality Control, Filter, and Rescue modules.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.4)

_______

## SQANTI3 v5.3 (04/12/2024)

#### Minor fixes/enhancements:
- QC class updates, bugfixes in test suites, and general performance improvements.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.3.0)

_______

## SQANTI3 v5.2 (04/10/2023)

#### Major changes:
**QC**
- Added read coverage threshold (>3 reads) for TSS ratio computation.
- Added `--ratio_TSS_metric` argument for metric selection for TSS ratio computation.

**ML filter**
- Expanded report with performance metrics: resulting probability distribution, confusion matrix, test set statistics.
- Minimum number of TP/TN isoforms to run ML filter is now 250.
- Implemented automatic exclusion of dist_* cols when RM are TP
- 'ML' arg to select ML filter mode is now **lowercase**, i.e. 'ml'.

**Rescue**
- Added `--mode` argument to rescue to allow selection of `automatic` or `full`, modes (default: `automatic`)

#### Minor fixes/enhancements:
- Added output directory creation and prefix checks in Rescue.
- Updated filter report syntax to match ggplot2 and tidyverse updates (ggplot>=3.4.0 now required for report).
- Update requirements to R>=4.3.0.
- RColorConesa now installed via CRAN.
- SQ3 version now correctly reported by all scripts.
- Fixed bugs in filter during QC output file reading (fasta, etc.).
- Minor bugs/typos.


> [Download release](Dependencies-and-installation.md#1-downloading-sqanti3)

_______

### Patch to SQANTI3 v5.1.2 (20/07/2023)

**Changes:**
* Speed improvements in **rescue** when finding best match IDs for rescue candidates.
* Output all reference transcripts associated to FSM/ISM found during **QC** if the `--isoform_hits` flag is supplied.
* Fixed bug in **ML filter** leading to errors when input data had no mono-exonic transcripts.
* Fixed saturation curve bug in pigeon report.
* Fixes conda environment installation of bcbio-gff.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.1.2)

### Patch to SQANTI3 v5.1.1 (19/01/2023)

**Changes:**
* Adapt for pigeon compatibility.
* Bug fixes in STM function.
* Minor fixes in ML filter documentation and output file handling.
* Fixed bug leading to incorrect classification of **genic intron** transcripts.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.1.1)


## SQANTI3 v5.1 (22/07/2022)

#### Major changes:
* Implemented new **rescue strategy** to recover transcriptome diversity lost after filtering (see details at the [SQ rescue wiki](Running-SQANTI3-rescue.md)).
* Updated **conda environment** to include rescue dependencies. We recommend creating the environment again in order for SQANTI3 to run without error.
* Fixed behavior of **mono-exon transcripts** during **ML filter**:
  - FSM now undergo intra-primming evaluation if they are mono-exons.
  - Corrected ML filter output when `--force_multi_exon` option is supplied: mono-exon transcripts will now be labeled as Artifacts.
* Fixed reasons file output by **rules filter**: the table now includes correct filtering reasons for **mono-exon transcripts**.
* Added an option to rules filter to control for mono-exon transcripts (previously available in ML filter).
* Modified the **output of SQANTI3 QC** to incorporate the creation of a complete `params.txt` file, i.e. including all arguments and the full paths of all supplied files.

 #### Minor fixes/enhancements:
   - Fixed output path for IsoAnnotLite GFF3 that prevented writing the file to the correct output directory when -gff3 option was not used.
   - Set temporary file dir for HTML report creation (fixes Singularity container error).
   
> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.1)

___________


## SQANTI3 v5.0 (01/06/2022)

#### Major changes:
* Implemented new **machine learning-based filter**.
* Updated **rules filter**: users can now define their own set of rules using a JSON file. By default, the rules filter applies the same set of rules that were implemented in the old `sqanti3_RulesFilter.py` script.
 * The `sqanti3_RulesFilter.py` script is now deprecated and has been replaced by `sqanti3_filter.py`, which works a wrapper for both filters (see details in the [documentation](Running-SQANTI3-filter.md)).
 * IsoAnnotLite updated to version 2.7.3.
 * Substantial modification of the SQANTI3 directory structure, with `utilities` folder now being divided into subfolders that group the scripts by their function.
 * Added a column in the classification file to indicate whether a polyA motif was found, which adds to the existing column detailing the detected motif (details [here](https://github.com/ConesaLab/SQANTI3/issues/138)).
* Changed CAGE argument and CAGE/polyA columns to capital letters (for consistency across columns and arguments).
* The `example` folder now includes sample commands and output files for SQANTI3 QC, rules filter and machine learning filter.
* Added new supported transcript model (STM) plots to the SQANTI3 QC report.

#### Minor fixes/enhancements:
   * Included cython (cDNA_cupcake dependency) as a dependency in the SQANTI3 conda environment.
   * pip installed in conda environment.
   * When supplied, the new `sqanti3_filter.py` filters the `sqanti3_qc.py` output files using the filter result (rules or ML). This was not previously done by `sqanti3_RulesFilter.py`.
   * Antisense vs intergenic bug: fixed inconsistencies in classification of isoforms across the two categories.
   * Fixed deprecation warnings in calculation of ratioTSS.
   * Minor report updates.

> [Download release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.0)