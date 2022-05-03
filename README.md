![SQANTI3 logo](https://github.com/FJPardoPalacios/public_figures/blob/master/sq3-logo.png)

#		 SQANTI3

SQANTI3 is the newest version of the SQANTI tool ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/)) that merges features from SQANTI, ([code repository](https://github.com/ConesaLab/SQANTI)) and SQANTI2 ([code repository](https://github.com/Magdoll/SQANTI2)), together with new additions. SQANTI3 will continue as an integrated development aiming to provide the best characterization for your new long read-defined transcriptome. 

SQANTI3 is the first module of the [Functional IsoTranscriptomics (FIT)](https://tappas.org/) framework, that also includes IsoAnnot and tappAS.

## Lastest updates
Current version (12/04/2022): SQANTI3 version 4.3

New features implemented in SQANTI3

**version 4.3:**
* Fixed bug printing the ORF in the classification file
* Fixed bug using PolyA peaks
* IsoAnnot Lite updated.

**version 4.2:**
* Extra aligner for mapping Long-Read sequences: **uLTRA** . To use it, set the option `--aligner_choice uLTRA`.
* NMD prediction now takes into account the 50bp rule, meaning that isoforms will be tagged as `TRUE` for `prediction_NMD` if there's a predicted ORF and CDS ends at least 50bp before the last junction.
* Fixed minor bugs

**version 4.1:**
* Change for `--gtf` argument. Now, by default, SQANTI3 will expect a GTF file as input. It will still be possible to provide the sequences of your isoforms as a FASTA or FASTQ file by activatingt the `--fasta` option.
* New plots and minor bugs of the report fixed.

**version 4.0:**
* Creation of HTML report. Using the `--report` argument, it is possible to choose which type of report: `html` (default), `pdf`, `both` or `skip`.
* Short-reads processing pipeline included. Now, providing directly your short-read data (FASTA/FASTQ format), SQANTI3 will:
    * Map them against the genome to identify SJ and calculate their coverage.
    * Calculate the "ratio_TSS" value for each isoform of your transcriptome.
    * If pair-end data is provided, isoform expression will be computed using kallisto.
* New subcategories for FSM: Reference match, Alternative 3' UTR, Alternative 5' UTR and Alternative 5' and 3' UTRs.
* IsoAnnotLite implemented to generate tappAS compatible GFF3 files. GFF3 output may incorporate functional annotation labels for model species supported by tappAS.
* New plots:
    *  Saturation curves only plot when `--saturation` option activated.
* Installation provided as a conda yml environment file  

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

SQANTI3 paper is in preparation, in the meantime it is possible to cite the [original SQANTI paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/).

*Tardaguila M, de la Fuente L, Marti C, et al. SQANTI: extensive characterization of long-read transcript sequences for quality control in full-length transcriptome identification and quantification. Genome Res. 2018;28(3):396-411. doi:10.1101/gr.222976.117*

