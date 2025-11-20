![SQANTI3 logo](https://github.com/ConesaLab/blob/master/src/utilities/data/sqanti3_logo.png?raw=true)

# SQANTI3

![GitHub Release](https://img.shields.io/github/v/release/ConesaLab/SQANTI3)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/ConesaLab/SQANTI3)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues-closed/ConesaLab/SQANTI3)
![GitHub Repo stars](https://img.shields.io/github/stars/ConesaLab/SQANTI3)


SQANTI3 is the newest version of the [SQANTI tool](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) that merges features from [SQANTI](https://github.com/ConesaLab/SQANTI) and  [SQANTI2](https://github.com/Magdoll/SQANTI2), together with new additions. SQANTI3 will continue as an integrated development aiming to provide the best characterization for your new long read-defined transcriptome. 

SQANTI3 is the first module of the [Functional IsoTranscriptomics (FIT)](https://tappas.org/) framework, which also includes [IsoAnnot](https://github.com/ConesaLab/IsoAnnot) and [tappAS](https://github.com/ConesaLab/tappAS).

SQANTI3 also includes TUSCO (Transcriptome Universal Single-isoform COntrol), a curated single-isoform reference for benchmarking transcriptome reconstruction from long-read sequencing; see preprint https://doi.org/10.1101/2025.08.23.671926.

## Installation
The [latest SQANTI3 release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.5.3) (17/11/2025) is **version 5.5.3**. See our wiki for [installation instructions](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation).

For information about previous releases and features introduced in them, see the [version history](https://github.com/ConesaLab/SQANTI3/wiki/Version-history).

**⚠️WARNING:** v5.0 represented a major release of the SQANTI3 software. **Versions of SQANTI3 >= 5.0 will not have backward compatibility** with previous releases and their output (v4.3 and earlier). Users that wish to apply any of the new functionalities in v5.0 to output files from older versions will therefore need to re-run SQANTI3 QC. See below for a full list of changes implemented in SQANTI3 v5.0.

**⚠️WARNING:** v5.4 has changed the naming of some of the arguments in the command line input of [quality control](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control#arguments-and-parameters-in-sqanti3-qc), [filter](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter#running-the-rules-filter) and [rescue](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-rescue#running-sq3-rescue). Please, check their respective documentation for more information.

**⚠️WARNING:** v5.5 has replaced GeneMarkS-T for TransDecoder2. This has improved the results at the cost of longer runtimes. However, the changes are expected to be minimal. **The conda environment has to be recreated to add TD2 distribution.** Release 5.5.1 updated the python version to 3.11.13, as on the previous version freezes might happen in some datasets.

## Documentation

For detailed documentation, please visit [the SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki).

### Wiki contents:
* [Introduction to SQANTI3](https://github.com/ConesaLab/SQANTI3/wiki/Introduction-to-SQANTI3)

* [Dependencies and installation](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation)

* [Version history](https://github.com/ConesaLab/SQANTI3/wiki/Version-history)

* [Isoform classification: categories and subcategories](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories)

* [Running SQANTI3 from the wrapper](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-from-the-wrapper)

* [Running SQANTI3 quality control](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control)
* [TUSCO quick start (SQANTI3 QC)](https://github.com/ConesaLab/SQANTI3/wiki/TUSCO-quick-start)

* [Understanding the output of SQANTI3 QC](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC)

* [IsoAnnotLite](https://github.com/ConesaLab/SQANTI3/wiki/IsoAnnotLite)

* [Running SQANTI3 filter](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter)

* [Running SQANTI3 rescue](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-rescue)

* [Running SQANTI-reads](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI%E2%80%90reads)

* [Tutorial: running SQANTI3 on an example dataset](https://github.com/ConesaLab/SQANTI3/wiki/Tutorial:-running-SQANTI3-on-an-example-dataset)

Please, note that we are currently updating and expanding the wiki to provide as much information as possible and 
enhance the SQANTI3 user experience. Pages under construction -or where information is still missing- will be indicated where appropriate. 
Thank you for your patience!


## How to cite SQANTI3
If you are using SQANTI3 in your research, please cite the following paper in addition to this repository:

- Pardo-Palacios, F.J., Arzalluz-Luque, A. et al. **SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms**. *Nat Methods* (2024). https://doi.org/10.1038/s41592-024-02229-2

- Keil, N., Monzó, C., McIntyre, L., Conesa, A. **Quality assessment of long read data in multisample lrRNA-seq experiments with SQANTI-reads**. *Genome Res* (2025). https://doi.org/10.1101/gr.280021.124

- Liu, T., Paniagua, A., Jetzinger, F., Ferrández-Peral, L., Frankish, A., Conesa, A. **Transcriptome Universal Single-isoform COntrol: A Framework for Evaluating Transcriptome reconstruction Quality**. *bioRxiv* (2025). https://doi.org/10.1101/2025.08.23.671926
