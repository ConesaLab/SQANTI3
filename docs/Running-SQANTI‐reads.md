#  Running SQANTI-reads  

You can find out more about SQANTI-reads here: Keil N, Monzó C, McIntyre L, Conesa A (2025). Quality assessment of long read data in multisample lrRNA-seq experiments with SQANTI-reads. Genome Res. DOI: 10.1101/gr.280021.124

SQANTI-reads leverages SQANTI3, a tool for the analysis of the quality of transcript models, to develop a quality control protocol for replicated long-read RNA-seq experiments. The number/distribution of reads, as well as the number/distribution of unique junction chains (transcript splicing patterns), in SQANTI3 structural categories are compiled. Multi-sample visualizations of QC metrics can also be separated by experimental design factors. We introduce new metrics for 1) the identification of potentially under-annotated genes and putative novel transcripts and 2) variation in junction donors and acceptors.
  
## Introduction to SQANTI-reads  
  
SQANTI-reads has two modes of running: _fast_ and _simple_.  
  
### SQANTI-reads _fast_ mode (RECOMMENDED). 
In the _fast_ mode, SQANTI-reads takes as (minimum) input i) a design CSV file with a column **sampleID** that includes the names of the samples (as we want them to be represented in the output plots), and a column **file_acc** with the name of the directory where the output SQANTI3 is stored. And ii) the reference annotation.gtf file.

An example design file is:

| sampleID | file_acc |
| --------- | ------- |
| wtc11_PBcDNA | ENCFF003QZT |
| wtc11_PBCapTrap | ENCFF023EXJ |
| wtc11_ONTR2C2 | ENCFF063ASB |
| wtc11_ONTdRNA | ENCFF104BNW |
| wtc11_cDNA | ENCFF105WIJ |

To run SQANTI-reads in _fast_ mode, you must first, run SQANTI3-QC on all your samples individually (great opportunity for you to parallelize this using your system's slurm, gnu parallel or your favourite parallelization method). Using the parameter ```--input_dir```, you can give sqanti3_reads.py the path to the parent directory where all the SQANTI3-QC output directories are stored.

We recommend using _fast_ mode, this ensures you can easily parallelize your work and have more control over mapping parameters. The usual reads pre-processing pipelines for ONT and PacBio to run SQANTI-reads include:  
- For ONT: Running pychopper to strand the reads, remove primers and polyA tails.  
- For PacBio: Running lima and refine from IsoSeq pipeline, and ```bamtools convert -format fastq -in fl.bam``` to transform PacBio FL bams into fastq files for mapping.  
- Common next steps: mapping with minimap2, transforming to gtf using spliced_bam2gff and running SQANTI3-QC.

Directories with SQANTI3-QC output are named as {file_acc} (e.g. ./ENCFF003QZT/) and SQANTI3-QC output files are named as {sampleID} (e.g. ./ENCFF003QZT/wtc11_PBcDNA_classification.txt)

Example run starting after running spliced_bam2gff:
```
python sqanti3_qc.py ENCFF003QZT.gff --annotation hg38.ensGene.gtf --genome hg38.fa --skipORF --min_ref_len 0 --aligner_choice minimap2 -t 8 -d ./ENCFF003QZT -o wtc11_PBcDNA
python sqanti3_reads.py --design design.csv --annotation hg38.ensGene.gtf
```
  
### SQANTI-reads (minimum) _simple_ mode.  
In the _simple_ mode, SQANTI-reads takes as (minimum) input a design CSV file with a column **sampleID** that includes the names of the samples (as we want them to be represented in the output), and a column **file_acc** with the name of the FASTQ or GTF/GFF files. It also requires the reference genome fasta, and reference annotation gtf.

When running in _simple_ mode, it will call SQANTI3-QC pipeline sequentially for each sample to:  
- If FASTQ files were inputted, call SQANTI3-QC and first map using minimap2 (by default, but can be changed to uLTRA), and second transform to GTF for isoform classification.  
- If GTF files were inputted, it will call SQANTI3-QC directly for isoform classification.  
  
Once the SQANTI3-QC pipeline is run within SQANTI-reads, it will automatically fill two columns in the design CSV file, **classification_file** and **junction_file**. These columns will include the names of the files generated for each **sampleID**. Each of these generated classification and junction files, will be stored within one sub-directory per sample, named as the **file_acc**. Next, it will calculate all the SQANTI-reads metrics. Finally, it will generate the plots and tables summarizing the metrics results.

Where design.csv is:

| sampleID | file_acc |
| --------- | ------- |
| wtc11_PBcDNA | ENCFF003QZT |
| wtc11_PBCapTrap | ENCFF023EXJ |
| wtc11_ONTR2C2 | ENCFF063ASB |
| wtc11_ONTdRNA | ENCFF104BNW |
| wtc11_cDNA | ENCFF105WIJ |

Example run:
```
python sqanti3_reads.py --design design.csv --annotation hg38.ensGene.gtf --genome hg38.fa
```

Fastq files are named {file_acc}\*.fastq (e.g. ENCFF003QZT_PB.fastq) and are stored in the current directory.

## Getting ready

Before running SQANTI-reads, you will need to:

**Activate the SQANTI3 conda environment:**

```
(base)-bash-4.1$ conda activate SQANTI3.env
(SQANTI3.env)-bash-4.1$
```

## Arguments and parameters in SQANTI-reads

The SQANTI-reads script accepts the following arguments:

```
usage: sqanti3_reads.py [-h] [--genome GENOME] --annotation ANNOTATION -de INDESIGN [-i INPUT_DIR] [-f INFACTOR]
                       [-p PREFIX] [-d DIR] [--min_ref_len MIN_REF_LEN] [--force_id_ignore]
                       [--aligner_choice {minimap2,uLTRA}] [-t CPUS] [-n CHUNKS] [-s SITES] [-ge ANNOTEXP]
                       [-je JXNEXP] [-pc PERCCOV] [-pj PERCMAXJXN] [-fl FACTORLVL] [--all_tables] [--report {pdf,html,both}] [--pca_tables]
                       [--skip_hash] [--verbose] [-v]

Structural and Quality Annotation of Novel Transcript Isoforms

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       Reference genome (Fasta format).
  --annotation ANNOTATION
                        Reference annotation file (GTF format).
  -de INDESIGN, --design INDESIGN
                        Path to design file, must have sampleID and file_acc column.
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Path to directory where fastq/GTF files are stored. Or path to parent directory with
                        children directories of SQANTI3 runs. Default: Directory where the script was run.
  -f INFACTOR, --factor INFACTOR
                        This is the column name that plots are to be faceted by. Default: None
  -p PREFIX, --prefix PREFIX
                        SQANTI-reads output filename prefix. Default: sqantiReads
  -d DIR, --dir DIR     Directory for output sqanti_reads files. Default: Directory where the script was run.
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length. Default: 0 bp
  --force_id_ignore     Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)
  --aligner_choice {minimap2,uLTRA}
                        Default: minimap2
  -t CPUS, --cpus CPUS  Number of threads used during alignment by aligners. Default: 10
  -n CHUNKS, --chunks CHUNKS
                        Number of chunks to split SQANTI3 analysis in for speed up. Default: 1
  -s SITES, --sites SITES
                        Set of splice sites to be considered as canonical (comma-separated list of splice sites).
                        Default: GTAG,GCAG,ATAC.
  -ge ANNOTEXP, --gene_expression ANNOTEXP
                        Expression cut off level for determining underannotated genes. Default = 100
  -je JXNEXP, --jxn_expression JXNEXP
                        Coverage threshold for detected reference donors and acceptor. Default = 10
  -pc PERCCOV, --perc_coverage PERCCOV
                        Percent gene coverage of UJC for determining well-covered unannotated transcripts.
                        Default = 20
  -pj PERCMAXJXN, --perc_junctions PERCMAXJXN
                        Percent of the max junctions in gene for determining near full-length putative novel
                        transcripts. Default = 80
  -fl FACTORLVL, --factor_level FACTORLVL
                        Factor level to evaluate for underannotation
  --all_tables          Export all output tables. Default tables are gene counts, ujc counts, length_summary, cv
                        and and underannotated gene tables
  --pca_tables          Export table for making PCA plots  
  --skip_hash           Skip the hashing step
  --report {pdf,html,both}
                        Default: pdf
  --verbose             If verbose is run, it will print all steps, by default it is FALSE
  -v, --version         Display program version number.
```

## SQANTI-reads output

The sqanti-reads output is written to the path specified in ```--dir``` or ```-d``` argument, also appending the prefix provided via the ```--prefix``` or -p argument.

The following output files are generated after running it:

### Modified input files

New reads_classification.txt files will be generated, including all columns from the original files plus two additional columns:
- jxn_string: includes the unique junction chain string as ```chromosome_strand_junction1_junction2_junctionN```. Monoexons are coded as ```chromosome_strand_monoexon_readID```.
- jxnHash: includes the hashing of the jxn_string for easier work in next steps.

The design file will be modified to include all initial columns (e.g. sampleID, file_acc, platform, age, sex) and two additional columns:

- classification_file: will include the full paths to the reads_classification.txt files (output of SQANTI-reads).
- junction_file: will include the full paths to the junctions.txt files (output of SQANTI3).

### Results default .csv files

- sqantiReads_cv.csv: Provides metrics on the coefficient of variance of reference junctions, per sample.
- sqantiReads_gene_counts.csv: Provides the number of reads in each structural category, per gene, per sample.
- sqantiReads_length_summary.csv: Provides the number and percentage of reads in length categories per sample.
- sqantiReads_ujc_counts.csv: Provides a list of junction strings in each sample and the number of reads in each sample associated with each junction string.
- sqantiReads_gene_classification.csv: For genes with coverage meeting a user defined threshold (-ge), provides the annotation category of each gene.
- sqantiReads_putative_novel_transcripts.csv: Provides metrics on NIC and NNC UJCs and flags putative novel transcripts.

### Results optional .csv files

- sqantiReads_cv_acc/don_counts.csv: Provides the number of detected annotated donors and acceptors in each junction variation category.
- sqantiReads_cv_acc/don_counts.csv: Provides the number of detected annotated donors and acceptors in each junction variation category.
- sqantiReads_FSM/ISM/NIC_NNC_counts.csv: Provides the number of reads in each subcategory for FSMs, ISMs, NICs and NNCs.
- sqantiReads_err_counts.csv: Provides the number and percentage of reads with. evidence of intrapriming, RT-switching and non-canonical junctions per samples.
- sqantiReads_pca_loadings.csv: Gives loadings of PC1 and PC2 for the PCA analaysis.
- sqantiReads_pca_variance.csv: Gives the proportion variance explained by each PC.

### Results plots

- sqantiReads_plots.pdf: Pdf file with one output plot visualizing QC metrics per page.
- sqantiReads_annotation_plots.pdf: Pdf file with one output plot visualizing underannotation metrics per page.

# Running example SQANTI-reads

There is an example dataset with three small samples on the SQANTI3/example/sqanti_reads_test/ directory. If you run it from the main SQANTI3 directory, you can run it in simple mode with the following command:

```
python sqanti3_reads.py --design ./example/sqanti_reads_test/sqR_design_file.csv -i ./example/sqanti_reads_test/ -p SQ_R --annotation ./example/gencode.v38.basic_chr22.gtf --genome ./example/GRCh38.p13_chr22.fasta --dir ./example/sqanti_reads_test/ --report both
```

## Citing SQANTI-reads
If you are using SQANTI-reads, please cite:  
Keil N, Monzó C, McIntyre L, Conesa A (2025). Quality assessment of long read data in multisample lrRNA-seq experiments with SQANTI-reads. Genome Research, 35 (4), 987. DOI: 10.1101/gr.280021.124
  
SQANTI-reads is based, and uses SQANTI3, please also cite:  
Pardo-Palacios FJ, Arzalluz-Luque A, Kondratova L, Salguero P, Mestre-Tomás J, Amorín R, Estevan-Morió E, Liu T, Nanni A, McIntyre L, Tseng E, Conesa A (2024). SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. Nature Methods, 21, 793-797. DOI: 10.1038/s41592-024-02229-2