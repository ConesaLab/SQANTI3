![SQANTI3 logo](https://github.com/FJPardoPalacios/public_figures/blob/master/sq3-logo.png)

# SQANTI3

SQANTI3 is the newest version of the [SQANTI tool](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) that merges features from [SQANTI](https://github.com/ConesaLab/SQANTI) and  [SQANTI2](https://github.com/Magdoll/SQANTI2), together with new additions. SQANTI3 will continue as an integrated development aiming to provide the best characterization for your new long read-defined transcriptome. 

SQANTI3 is the first module of the [Functional IsoTranscriptomics (FIT)](https://tappas.org/) framework, which also includes IsoAnnot and tappAS.

## Latest updates
Latest SQANTI3 release (11/05/2022) is **version 5.0**.

**WARNING:** v5.0 constitutes a major release of the SQANTI3 software. **Versions of SQANTI3 >= 5.0 will not have backward compatibility** with previous releases and their output (v4.3 and earlier). Users that wish to apply any of the new functionalities in v5.0 to output files from older versions will herefore need to re-run SQANTI3 QC.

New features implemented in SQANTI3 v5.0:

* Implemented new **machine learning-based filter**.
* Updated **rules filter**: users can now define their own set of rules using a JSON file. By default, the rules filter applies the same set of rules that were implemented in the old `sqanti3_RulesFilter.py` script.
 * The `sqanti3_RulesFilter.py` script is now deprecated and has been replaced by `sqanti3_filter.py`, which works a wrapper for both filters (see details in the [documentation](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter)).
 * IsoAnnotLite updated to version 2.7.3.
 * Substantial modification of the SQANTI3 directory structure, with `utilities` folder now being divided into subfolders that group the scripts by their function.
 * Added a column in the classification file to indicate whether a polyA motif was found, which adds to the existing column detailing the detected motif (details [here](https://github.com/ConesaLab/SQANTI3/issues/138)).
* Changed CAGE argument and CAGE/polyA columns to capital letters (for consistency across columns and arguments).
* The `example` folder now includes sample commands and output files for SQANTI3 QC, rules filter and machine learning filter.
* Added new supported transcript model (STM) plots to the SQANTI3 QC report.
* Minor fixes/enhancements:
   * Included cython (cDNA_cupcake dependency) as a dependency in the SQANTI3 conda environment.
   * pip installed in conda environment.
   * When supplied, the new `sqanti3_filter.py` filters the `sqanti3_qc.py` output files using the filter result (rules or ML). This was not previously done by `sqanti3_RulesFilter.py`.
   * Antisense vs intergenic bug: fixed inconsistencies in classification of isoforms across the two categories.
   * Fixed deprecation warnings in calculation of ratioTSS.
   * Minor report updates.


## Documentation

For detailed documentation, please visit [the SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki).

### Wiki contents:
* [Introduction to SQANTI3](https://github.com/ConesaLab/SQANTI3/wiki/Introduction-to-SQANTI3)
* [Dependencies and installation](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation)
* [Isoform classification: categories and subcategories](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories)
* [Running SQANTI3 quality control](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control)
* [Understanding the output of SQANTI3 QC](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC)
* [Running SQANTI3 filter](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter)
* [Tutorial: running SQANTI3 on an example dataset](https://github.com/ConesaLab/SQANTI3/wiki/Tutorial:-running-SQANTI3-on-an-example-dataset)

Please, note that we are currently updating and expanding the wiki to provide as much information as possible and 
enhance the SQANTI3 user experience. Pages under construction -or where information is still missing- will be indicated where appropriate. 
Thank you for your patience!


## How to cite SQANTI3

SQANTI3 paper is currently in preparation. In the meantime, when using SQANTI3 in your research, please cite the [original SQANTI paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) as well as this repository:

- Tardaguila M, de la Fuente L, Marti C, et al. SQANTI: extensive characterization of long-read transcript sequences for quality control in full-length transcriptome identification and quantification. *Genome Res*, 2018. **28**(3):396-411. doi:10.1101/gr.222976.117

