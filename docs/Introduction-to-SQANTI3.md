### Table of contents:

- <a href="#intro">Introduction</a>

    * <a href="#fit">SQANTI3 and the Functional IsoTranscriptomics (FIT) pipeline</a>
    * <a href="#workflow">Before running SQANTI3: recommended long read processing workflow</a>

- <a href="#sqanti">How does SQANTI3 work?</a>

- <a href="#structure">SQANTI3 structure</a>

    * <a href="#source">Source folder
    * <a href="#utilities">Utilities folder</a>
    * <a href="#data">Data folder</a>
    * <a href="#example">Example folder</a>

__________________________________

<a name="intro"/>

## Introduction

SQANTI3 is the newest version of the [SQANTI tool](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/). 
SQANTI3 combines features from the original [SQANTI](https://github.com/ConesaLab/SQANTI) and from
[SQANTI2](https://github.com/Magdoll/SQANTI2), as well as newly implemented functionalities and 
transcript features.

**Disclaimer:** please note that, although still available, both SQANTI and SQANTI2 are deprecated and will no longer be 
maintained/updated. All development efforts will continue in SQANTI3, aiming to providing the most
comprehensive characterization of long read-defined transcriptomes for the community.


<a name="fit"/>

### SQANTI3 and the Functional IsoTranscriptomics (FIT) pipeline

SQANTI3 constitutes the first module of the [Functional IsoTranscriptomics (FIT)](https://tappas.org/) 
pipeline, which is an end-to-end strategy to perform isoform-level bioinformatics analyses. 
The SQANTI3 tool is designed to enable **quality control and filtering** of long read-defined transcriptomes, which are often rich in artifacts and false-positive isoforms. Therefore, a good curation of the transcriptome is indispensable to proceed with FIT analysis and produce valid, biologically sound conclusions/hypothesis. 

After generating a high-quality transcriptome with SQANTI3, downstream steps include:

- **Functional annotation** of isoform models, including positionally-defined functional features such as motifs,
domains, etc. IsoAnnot, a tool for *de novo* annotation of isoforms, is currently under development, however,
users can run [IsoAnnotLite](https://github.com/ConesaLab/SQANTI3/wiki/IsoAnnotLite) within or outside of SQANTI3 to impute functional features from other already-annotated
transcriptomes.
- **Expression-based functional analysis** using [tappAS](https://app.tappas.org/). tappAS is a Java GUI application 
that leverages both expression and domain/motif annotation information to gain insight into the functional implications 
of alternative isoform expression.


<a name="workflow"/>

### Before running SQANTI3: recommended long read processing workflow

Here is our recommended workflow, including the best way to generate the SQANTI3 inputs and how to proceed after QC and filtering:

1. **Sample pooling**: while we are aware that some users may have long read data from several replicates and/or samples,
we recommend pooling all long read samples to build a single transcriptome per experiment.
2. **Long read data processing** using your preferred transcriptome-building tool. We do not recommend using SQANTI3 on 
raw long reads, as it is NOT designed as a tool for long read data QC.
3. **Collapse of isoform models**. Typically, long read data processing pipelines generate a large number of highly redundant
isoform models. We recommend collapsing these using tools such as [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) or 
[TAMACollapse](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) to reduce the number of isoforms and create unique isoform models prior to running SQANTI3.
4. **Quality control and filtering**: we strongly encourage users to do as careful an inspection of their long read-defined
transcriptomes as possible, including filtering their transcriptome to remove potential false positive isoforms, which are abundant
in long read-generated transcriptomes.
5. **Quantification** of the filtered transcriptome using short/long reads and your preferred tool. We do not recommend using the
expression estimates input into SQANTI3 for downstream analysis: these are used for quality control purposes only. Once all artifacts
are removed from the transcriptome, the reads can be used to obtain a more accurate quantification.



<a name="sqanti"/>

## How does SQANTI3 work?

SQANTI3 is a tool for in-depth characterization of isoforms obtained by full-length transcript sequencing, 
which are commonly returned in a fasta or GTF file format. SQANTI3 combines the long read-defined transcripts with the **reference annotation** as well as with other **orthogonal data** to provide a wide range of descriptors of transcript quality. SQ3 generates a comprehensive report to facilitate quality control and filtering of the isoform models.

SQANTI3 is mainly designed to perform two 
different tasks, both of them equally important:

1. **Isoform classification and quality control** ([SQANTI3 QC](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control)) for long read-defined transcriptomes. The SQ3 [categories and subcategories](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories), 
together with a long list of transcript-level attributes and descriptors, allow users to carefully inspect the properties of their isoform models, as well as identify potential problems generated during library preparation and raw data processing.
2. **Artifact filtering** ([SQANTI3 filter](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter)) for long read-defined transcriptomes. Using the large number of descriptors calculated by SQANTI3, users can make informed decisions to remove potential false positive isoforms from their transcriptomes. This is particularly relevant considering the biases and pitfalls of current long read sequencing protocols.

To gain insight into these two steps, we encourage reading the original [SQANTI publication](https://genome.cshlp.org/content/early/2018/02/09/gr.222976.117). Recently, however, we have implemented a final step in the SQANTI3 workflow:

3. **Reference transcript rescue** ([SQANTI3 rescue](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-rescue)) to find matching reference transcript for discarded artifacts. This module intends to keep the diversity in the transcriptome, i.e. preventing loss of transcripts with (mostly) valid junction chains for which a long-read defined isoform could not be validated or even the removal of genes for which all isoforms were catalogued as artifacts by the filter, in spite of having long read-based evidence of expression. By running the rescue, SQANTI3 will select a set of reference transcripts that the discarded artifacts can be confidently assigned to, and add them to expand the filtered transcriptome.

![Sqanti3 workflow](https://github.com/FJPardoPalacios/public_figures/blob/master/SQ3_qc.png)



<a name="structure"/>

## SQANTI3 structure

At the moment, the SQANTI3 tool is directly called using the wrapper `sqanti3`. This main wrapper orchestrates the execution of the **three primary modules**. Each module is designed to perform a specific function within the SQANTI3 workflow, implemented as separate Python scripts:

- `sqanti3_qc.py` for quality control.
- `sqanti3_filter.py` for transcriptome filtering/curation.
- `sqanti3_rescue.py` to replace artifacts by their closest matching reference transcript.

To maintain a clean and organized codebase, the detailed implementation for each core module and the auxiliary tools used in the workflow are housed within the `src` (source) directory. _As of the [latest release](https://github.com/ConesaLab/SQANTI3/releases/latest), only sqanti3_qc.py has been fully modularized._

<a name="source"/>

### Source folder

The `src` directory contains all the minor modules that each one of the sqanti blocks and the parser use to function properly. Each script has functions or classes according to their function within SQANTI3: classification, utilities, parsers, QC classes, etc.


<a name="utilities"/>
#### Utilities folder

The `utilities` folder contains all the auxiliary scripts and functions required to run `sqanti3_qc.py` and `sqanti3_filter.py`.
Under the main directory, you will find:

- `IsoAnnotLite_SQ3.py`: [IsoAnnotLite](https://isoannot.tappas.org/isoannot-lite/) script. Used to generate functional annotations and tappAS-compatible outputs.
- gtfToGenePred: UCSC conversion tool, obtained (see details [here](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)).
Used during QC.
- `indels_annot.py`: when a fasta file is provided as input, corrects indels in transcripts (note that SQANTI3 QC requires a GTF by default, so this function is rarely used).
- `rt_switching.py`: function to check splice junctions for potential RT switching.
- `short_reads.py`: used to process raw short-read data when supplied to `sqanti3_qc.py`.

##### Filter

The `filter` subfolder contains two R scripts, each for one of the [filters available in SQANTI3](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter) (i.e. rules and machine learning-based). 
In addition, a JSON file containing the default set of rules for the rules filter (which are equivalent to running the old RulesFilter in SQANTI2 and SQANTI3 for versions < 5.0).

##### Rescue 

The `rescue` subfolder contains a series of R scripts implementing the different [rescue steps](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-rescue). These are run internally by the rescue wrapper (`sqanti3_rescue.py`) to perform the rescue.

##### Reports

- The `report_qc` folder contains the main R script and auxiliary plotting functions/scripts used to generate the **SQANTI3 QC report**.
- The `report_filter` folder contains the main R script and auxiliary plotting functions used to generate the **SQANTI3 filter report**.

##### GMST

GeneMarkST tool. Used internally during QC for ORF prediction.

<a name="data"/>

### Data folder

The `data` folder contains pre-computed TSS/TTS orthogonal data that is ready to be used to run SQANTI3 QC (see [additional inputs](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control#minimal) section), avoiding the need for users to download them from reference databases.

- `ref_TSS_annotation` subfolder: human and mouse CAGE data from the [refTSS database](http://reftss.clst.riken.jp/reftss/Main_Page).
- `polyA_motif` subfolder: contains a list of human and mouse canonical polyadenylation motifs. 

<a name="example"/>

### Example folder

This folder contains all the necessary inputs to run SQANTI3 on an example dataset. This data was obtained from **public PacBio data**, in particular, it consists in a set of transcripts generated using long reads from the **Universal Human Reference RNA (UHR)** (see dataset details [here](https://github.com/Magdoll/cDNA_Cupcake/wiki/Iso-Seq-Public-Datasets#4-universal-human-reference-rna-uhrr-2021-release)). For convenience, we have only selected those in chromosome 22. 
This small dataset can be used to test SQANTI3 and/or to follow the 
[SQANTI3 example tutorial](https://github.com/ConesaLab/SQANTI3/wiki/Tutorial:-running-SQANTI3-on-an-example-dataset) in the Wiki. 

In addition, we provide the outputs generated after running QC and filters (both machine learning and rules) on the example dataset under the corresponding `*_output` subfolders.

### Test folder

Newly implemented, and in development, test suite for SQANTI3. This test suite is being developed as SQANTI3 is being modularized, with the aim to keep SQANTI3 as stable as possible in future releases. If you find any issue with the test suite, add it as an issue and it will be updated as soon as possible.
