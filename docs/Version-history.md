# Introduction

This page contains a detailed time-aware (DD/MM/YYYY) account of **SQANTI3 releases (>5.0)** and the new features and changes introduced in each of them.

Please be aware that v5.0 represented a major release of the SQANTI3 software. Versions of SQANTI3 >= 5.0 do not have backward compatibility with previous releases and their output (v4.3 and earlier). Users that wish to apply any of the new functionalities in v5.0 to output files from older versions will herefore need to re-run SQANTI3 QC.

See the dedicated site for [installation instructions](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation).

| [See all SQANTI3 releases](https://github.com/ConesaLab/SQANTI3/tags)

## SQANTI3 v5.2 [LATEST, 04/10/2023]

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


> [Download release](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation#1-downloading-sqanti3)

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
* Implemented new **rescue strategy** to recover transcriptome diversity lost after filtering (see details at the [SQ rescue wiki](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-rescue)).
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
 * The `sqanti3_RulesFilter.py` script is now deprecated and has been replaced by `sqanti3_filter.py`, which works a wrapper for both filters (see details in the [documentation](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter)).
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