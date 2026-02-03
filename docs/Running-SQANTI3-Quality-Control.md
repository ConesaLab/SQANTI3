### Table of contents:

* <a href="#ready">Getting ready</a>

* <a href="#args">Arguments and parameters in SQANTI3 QC</a>

* <a href="#running">Running SQANTI3 QC</a>

    * <a href="#minimal">Minimal input (mandatory)</a>

    * <a href="#optional">Additional inputs (optional)</a>

* <a href="#SR">Providing short reads</a>

    * <a href="#fasta">Raw short read data (FASTA/FASTQ files)</a>

    * <a href="#exp">Short read isoform expression (Kallisto/RSEM expression files)</a>

    * <a href="#cov">Splice junction coverage by short reads (SJ.out.tab files)</a>

    * <a href="#bam">Short read mapping results (BAM files)</a>

* <a href="#cage">Incorporating CAGE peak data</a>
        
* <a href="#polya">Incorporating polyA information</a>

    * <a href="#polyamotif">PolyA motif data</a>

    * <a href="#polyasite">PolyA site data</a>
        
* <a href="#flcount">Supplying single or multi-sample full-length (FL) counts</a>

* <a href="#tssratio">Selecting the metric used to aggregate TSS ratio across samples</a>
        
* <a href="#gff3">IsoAnnotLite and SQANTI3: using an existing tappAS GFF3 annotation file</a>

* <a href="#parallelization">Parallelization</a>


______________________________________________________________________

<a name="ready"/>

## Getting ready
</a>
Before running SQANTI3, you will need to:

**Activate the SQANTI3 conda environment:**

```
(base)-bash-4.1$ conda activate SQANTI3.env
(SQANTI3.env)-bash-4.1$
```

<a name="args"/>

## Arguments and parameters in SQANTI3 QC
</a>
The SQANTI3 quality control script accepts the following arguments:

```bash
usage: sqanti3_qc.py [-h] --isoforms ISOFORMS --refGTF REFGTF --refFasta
                     REFFASTA [--min_ref_len MIN_REF_LEN] [--force_id_ignore]
                     [--fasta] [--genename] [--short_reads SHORT_READS]
                     [--SR_bam SR_BAM] [--novel_gene_prefix NOVEL_GENE_PREFIX]
                     [--aligner_choice {minimap2,deSALT,gmap,uLTRA}]
                     [-x GMAP_INDEX] [-s SITES] [--skipORF]
                     [--orf_input ORF_INPUT] [--CAGE_peak CAGE_PEAK]
                     [--polyA_motif_list POLYA_MOTIF_LIST]
                     [--polyA_peak POLYA_PEAK] [--phyloP_bed PHYLOP_BED]
                     [-o OUTPUT] [-d DIR] [--saturation]
                     [--report {html,pdf,both,skip}] [--isoform_hits]
                     [--ratio_TSS_metric {max,mean,median,3quartile}]
                     [-t CPUS] [-n CHUNKS] [-l {ERROR,WARNING,INFO,DEBUG}]
                     [--is_fusion] [-e EXPRESSION] [-c COVERAGE] [-w WINDOW]
                     [-fl FL_COUNT] [-v] [--isoAnnotLite] [--gff3 GFF3]
```

<details>
  <summary>Arguments explanation</summary>

```bash
Structural and Quality Annotation of Novel Transcript Isoforms

options:
  -h, --help            show this help message and exit

Required arguments:
  --isoforms ISOFORMS   Isoforms (FASTA/FASTQ) or GTF format. It is
                        recommended to provide them in GTF format, but if it
                        is needed to map the sequences to the genome use a
                        FASTA/FASTQ file with the --fasta option.
  --refGTF REFGTF       Reference annotation file (GTF format)
  --refFasta REFFASTA   Reference genome (Fasta format)

Customization and filtering:
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length (default: 0 bp)
  --force_id_ignore     Allow the usage of transcript IDs non related with
                        PacBio's nomenclature (PB.X.Y)
  --fasta               Use when running SQANTI by using as input a
                        FASTA/FASTQ with the sequences of isoforms
  --genename            Use gene_name tag from GTF to define genes. Default:
                        gene_id used to define genes
  --short_reads SHORT_READS
                        File Of File Names (fofn, space separated) with paths
                        to FASTA or FASTQ from Short-Read RNA-Seq. If
                        expression or coverage files are not provided,
                        Kallisto (just for pair-end data) and STAR,
                        respectively, will be run to calculate them.
  --SR_bam SR_BAM       Directory or fofn file with the sorted bam files of
                        Short Reads RNA-Seq mapped against the genome
  --novel_gene_prefix NOVEL_GENE_PREFIX
                        Prefix for novel isoforms (default: None)

Aligner and mapping options:
  --aligner_choice {minimap2,deSALT,gmap,uLTRA}
                        Select your aligner of choice: minimap2, deSALT, gmap,
                        uLTRA (default: minimap2)
  -x GMAP_INDEX, --gmap_index GMAP_INDEX
                        Path and prefix of the reference index created by
                        gmap_build. Mandatory if using GMAP unless -g option
                        is specified.
  -s SITES, --sites SITES
                        Set of splice sites to be considered as canonical, in
                        a comma separated list. (default: ATAC,GCAG,GTAG)

ORF prediction:
  --skipORF             Skip ORF prediction (to save time)
  --orf_input ORF_INPUT
                        Input fasta to run ORF on. By default, ORF is run on
                        genome-corrected fasta - this overrides it. If input
                        is fusion (--is_fusion), this must be provided for ORF
                        prediction.

Functional annotation:
  --CAGE_peak CAGE_PEAK
                        FANTOM5 Cage Peak (BED format, optional)
  --polyA_motif_list POLYA_MOTIF_LIST
                        Ranked list of polyA motifs (text, optional)
  --polyA_peak POLYA_PEAK
                        PolyA Peak (BED format, optional)
  --phyloP_bed PHYLOP_BED
                        PhyloP BED for conservation score (BED, optional)

Output options:
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files. (Default: Directory where
                        the script was run.)
  --saturation          Include saturation curves into report
  --report {html,pdf,both,skip}
                        Select report format: html, pdf, both, skip (default:
                        html)
  --isoform_hits        Report all FSM/ISM isoform hits in a separate file
  --ratio_TSS_metric {max,mean,median,3quartile}
                        Define which statistic metric should be reported in
                        the ratio_TSS column (default: max)

Performance options:
  -t CPUS, --cpus CPUS  Number of threads used during alignment by aligners.
                        (default: 10)
  -n CHUNKS, --chunks CHUNKS
                        Number of chunks to split SQANTI3 analysis in for
                        speed up (default: 1).
  -l {ERROR,WARNING,INFO,DEBUG}, --log_level {ERROR,WARNING,INFO,DEBUG}
                        Set the logging level INFO

Optional arguments:
  --is_fusion           Input are fusion isoforms, must supply GTF as input
  -e EXPRESSION, --expression EXPRESSION
                        Expression matrix (supported: Kallisto tsv)
  -c COVERAGE, --coverage COVERAGE
                        Junction coverage files (provide a single file, comma-
                        delmited filenames, or a file pattern, ex:
                        "mydir/*.junctions").
  -w WINDOW, --window WINDOW
                        Size of the window in the genomic DNA screened for
                        Adenine content downstream of TTS (default: 20)
  -fl FL_COUNT, --fl_count FL_COUNT
                        Full-length PacBio abundance file
  -v, --version         Display program version number.
  --isoAnnotLite        Run isoAnnot Lite to output a tappAS-compatible gff3
                        file
  --gff3 GFF3           Precomputed tappAS species specific GFF3 file. It will
                        serve as reference to transfer functional attributes

```
</details><br>

<a name="running"/>

## Running SQANTI3 QC
</a>
<a name="minimal"/>

### Minimal input (mandatory)
</a>
This are the minimal files that you will need to run SQANTI3 QC:

*   **Long read-defined transcriptome**. It can be obtained after processing data from any of the available Third Generation Sequencing techonologies like Iso-Seq (PacBio) or Oxford Nanopore (ONT). SQANTI3 is compatible with the output of any long read-based transcriptome building pipeline, such as IsoSeq3, TALON, or FLAIR. SQANTI3 accepts it in several formats such as FASTA, FASTQ and GTF:

    - **GTF (default)**: by default, SQANTI3 expects the transcriptome to be provided as a GTF file, and we recommend to stick to this format if your transcriptome construction pipeline allows it.
    - **FASTA/FASTQ**: if you provide your transcript sequences in FASTA or FASTQ format, you will also need to supply the `--fasta` option. In this case, a mapping step will be initially performed. SQANTI3 currently supports minimap2 (default), uLTRA, deSALT and GMAP for mapping. To select a mapper, just supply the name to the `--aligner_choice` argument.

*   **Reference transcriptome annotation** in GTF format. This file will be used as reference to describe the amount of novelty across long read transcripts. Some examples of reference transcriptomes for different species can be found in [GENCODE](https://www.gencodegenes.org/) or [CHESS](http://ccb.jhu.edu/chess/).

*   **Reference genome** in FASTA format, for example, hg38. Before running SQANTI3, make sure of the following:
    - The transcriptome annotation GTF must match the provided reference genome. 
    - The chromosome/scaffolds names must be same in the reference annotation and the reference genome.

**⚠️WARNINGS:** 
- Before running SQANTI3, it is strongly recommended to **collapse redundant transcript sequences** using [isoseq collapse](https://isoseq.how/classification/isoseq-collapse.html) or [TAMA](https://github.com/GenomeRIK/tama/wiki). To learn more about this, we suggest taking a look at our [recommended long read processing workflow](https://github.com/ConesaLab/SQANTI3/wiki/Introduction-to-SQANTI3#workflow).

- Note that, by default, SQANTI3 expects PacBio-formatted IDs (i.e. PB.XX.XX). Users whose data was not processed using IsoSeq3 should add the `--force_id_ignore` flag to override this behavior.


<a name="optional"/>

### Additional inputs (optional)
</a>
* **Short read data**: SQANTI3 QC uses short-read data to validate several aspects of long read defined transcripts, i.e. to compute gene/isoform expression, junction coverage and TSS ratio. To learn more about the different ways in which short reads can be supplied to SQANTI3 QC, see <a href="#SR">Provding short reads</a> section.

* **CAGE Peak data:** In SQANTI2, a [CAGE peak database for the hg38 genome](https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz) was provided, originally retrieved from [FANTOM5](http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/). Now, we recommend to use data from [refTSS](http://reftss.clst.riken.jp/reftss/Main_Page). To ease usage, some versions of this annotation for several species have been uploaded to the `data/ref_TSS_annotation` folder, and we will try to keep them updated for new genome releases.

* **polyA motif list:** A ranked list of polyA motifs to find upstream of the 3' end site. See [human polyA list](https://raw.githubusercontent.com/Magdoll/images_public/master/SQANTI2_support_data/human.polyA.list.txt) for an example. We have included a list of motifs that can be used for mouse and human data at the `data/polyA_motifs/` folder. Please, if you know which are the most likely polyA motifs for other type of organism/clade, we would highly appreciate if you let us know.

* **FL count information:** See <a href="#flcount">FL count section</a> to include Iso-Seq FL count information for each isoform.

* **tappAS-annotation file:** when the `--isoAnnotLite` flag is activated, a GFF3 file containing isoform-level functional annotations for a reference transcriptome (e.g. Ensembl) can also be supplied via the `--gff3` flag.

* **[Intropolis](https://github.com/nellore/intropolis/blob/master/README.md) junction BED file:** In previous versions of SQANTI, it was provided a version of [Intropolis for hg38 genome, modified into STAR junction format](https://github.com/Magdoll/images_public/tree/master/SQANTI2_support_data) which is still valid.


If you don't feel like running the ORF prediction part, use `--skipORF`. Just know that all your transcripts will be annotated as non-coding.

If you have short read data, you can run STAR to get the junction file (usually called `SJ.out.tab`, see [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)) and supply it to SQANTI3 as is. If you have several samples, it is recommended to provide them as separated `*SJ.out.tab` files. 

By way of example, the following command was used to run SQANTI3 with the `example` data. The input files are a subdivision of the [Universal Human Reference (UHR) IsoSeq data](https://github.com/PacificBiosciences/DevNet/wiki/Sequel-II-System-Data-Release:-Universal-Human-Reference-(UHR)-Iso-Seq) that corresponds just to those polished and collapsed sequences located at chromosome 22. 

```bash
python sqanti3_qc.py --isoforms example/UHR_chr22.gtf \
                     --refGTF example/gencode.v38.basic_chr22.gtf \
                     --refFasta example/GRCh38.p13_chr22.fasta \
                     --CAGE_peak data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
                     --polyA_motif_list data/polyA_motifs/mouse_and_human.polyA_motif.txt    \
                     -o UHR_chr22 -d example/SQANTI3_output -fl example/UHR_abundance.tsv    \
                     --short_reads example/UHR_chr22_short_reads.fofn --cpus 4 --report both                     
```


It is *highly recommended* to run SQANTI3 using a GTF file with the Long-Read defined isoforms instead of input FASTA/FASTQ file with their sequences. If you ran [isoseq collapse](https://isoseq.how/classification/isoseq-collapse.html), the obtained `*collapsed.gff` file is actually in GFF2 format (equivalent to GTF), so you can use it as is. If you provide the sequences, a mapping step will be needed and the decisions took by the mapping algorithm will affect the final results. 


<a name="SR"/>

## Providing short reads
</a>
Short read data can be supplied to SQANTI3 in several ways:

* **Raw short-read data:** for users who may want SQANTI3 to compute all short read-based attributes internally without the need to perform any additional pre-processing, short read FASTA/FASTQ files can be provided using the dedicated argument (see <a href="#fasta">Raw short read data</a> section below).

Alternatively, users may supply pre-processed short read data in one or more of these formats:

* **Short read-based transcript expression:** see <a href="#exp">Short read isoform expression</a> section below for details on how to include pre-computed, short read-based expression files into SQANTI3 QC (RSEM or Kallisto output files or a user-defined expression matrix).

* **Splice junction coverage by short-reads:** short reads can be mapped to the genome using STAR to obtain junction coverage information, included in the the `SJ.out.tab` file (see <a href="#cov">Splice junction coverage section</a> for details).

* **Short read BAMs:** should be obtained by mapping short reads against the reference genome using a splice-aware mapper, such as STAR. SQANTI3 QC will use this information to compute the TSS ratio.


<a name="fasta"/>

### Raw short read data (FASTA/FASTQ files)
</a>
Matching RNA-Seq data (i.e. obtained from the same samples as your long-read data) can be supplied to SQANTI3 in the form of FASTA/FASTQ files via the `--short_reads` argument. The expected file format is a **File Of File Names** (.fofn), which is a text file that contains the paths of the Short-Read RNA-Seq data with one line per replicate and separated by one space in the case of paired-end data. It should look like this:

```text
/path/to/replicate1.R1.fastq /path/to/replicate1.R2.fastq
/path/to/replicate2.R1.fastq /path/to/replicate2.R2.fastq
```
[STAR](https://github.com/alexdobin/STAR) and [kallisto](https://pachterlab.github.io/kallisto/) will be run to automatically calculate splice-junction coverage, isoform expression and the `ratio_TSS` attribute, a metric that quantifies the differences of mean coverage 100bp upstream and downstream the reported TSS. 

Please note that, if you are going to use the `--short_reads` option, you should provide enough memory and computational resources for mapping and quantifying.

⚠️**Warning**: Kallisto will only be run if paired-end data is provided. If you wish to run Kallisto on single-end data, you may do it before running SQANTI3 and supply the isoform expression matrix as indicated below.



<a name="exp"/> 

### Short read isoform expression
</a>

Users can supply their pre-computed isoform expression via the `--expression` argument. Depending on whether this information is provided as one or multiple files, two formats are supported: 
- If you input **one expression file per sample**, you can provide several expression data files as a chain of comma-separated paths or by providing a directory where ONLY expression data is present. Accepted formats include the output of Kallisto and RSEM (see below).
- Alternatively, it is possible to provide a **pre-computed expression matrix**. The matrix should be a tab-separated file with transcripts as rows and samples/replicates as columns. Transcript identifiers must match the ones included in the input long read-defined transcriptome, and the name for that column must be `ID`. The rest of the columns in the header can be named as desired, and should correspond to the sample/replicate identifier.

#### Kallisto expression files

Kallisto expression files have the following format:

```
target_id   length  eff_length  est_counts  tpm
PB.1.1  1730    1447.8  0   0
PB.1.3  1958    1675.8  0   0
PB.2.1  213 54.454  0   0
PB.2.2  352 126.515 0   0
PB.3.1  153 40.3918 0   0
PB.4.1  1660    1377.8  0   0
PB.5.1  2767    2484.8  0   0
```

#### RSEM expression files

RSEM expression files have the following format:

```
transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct  posterior_mean_count    posterior_standard
_deviation_of_count     pme_TPM pme_FPKM        IsoPct_from_pme_TPM     TPM_ci_lower_bound      TPM_ci_upper_bound      FPKM_ci_lower_boun
d       FPKM_ci_upper_bound
PB.1.1  PB.1    1516    1369.11 8.00    0.05    0.22    100.00  8.00    0.00    0.06    0.25    100.00  0.0233365       0.0984532       0.
0981584 0.413561
PB.10.1 PB.10   1101    954.11  0.00    0.00    0.00    0.00    0.00    0.00    0.01    0.04    100.00  4.80946e-08     0.0283585       2.
01809e-07       0.11905
PB.100.1        PB.100  979     832.11  0.00    0.00    0.00    0.00    6.62    6.60    0.08    0.35    4.00    3.51054e-06     0.241555 1.47313e-05      1.01397
PB.100.2        PB.100  226     81.11   20.18   2.26    9.47    100.00  16.84   5.20    1.99    8.36    96.00   0.559201        3.45703 2.
34572   14.5141
```

<a name="cov"/>

### Splice junction coverage by short reads
</a>
If you have short read data, you can run STAR to get the **junction file** (usually called `SJ.out.tab`) and supply it to SQANTI3 via the `-c` argument. It is possible to use a different mapper, however, SQANTI3 requires that the format of the junction is similar to that of STAR's `SJ.out.tab` file (see below). 

As described in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf):

> SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that
STAR defines the junction start/end as intronic bases, while many other software define them as
exonic bases. 
> The columns in `SJ.out.tab` have the following meaning:
> - column 1: chromosome
> - column 2: first base of the intron (1-based)
> - column 3: last base of the intron (1-based)
> - column 4: strand (0: undefined, 1: +, 2: -)
> - column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
> - column 6: 0: unannotated, 1: annotated in the splice junctions database. Note that in 2-pass mode, junctions detected in the 1st pass are reported as annotated, in addition to annotated junctions from GTF.
> - column 7: number of uniquely mapping reads crossing the junction
> - column 8: number of multi-mapping reads crossing the junction
> - column 9: maximum spliced alignment overhang
> The filtering for this output file is controlled by the `--outSJfilter` parameters.

If you have several samples, we recommend providing them as separated `*SJ.out.tab` files instead of combined into a single file. To provide several files, you can follow different strategies: 
- Supply the path of the directory where all the `SJ.out.tab` are situated.
- Provide each of their individual paths separated by a comma.
- Use a wildcard to provide all paths at the same time, e.g. `/path/to/*my_suffix.tsv`.

<a name="bam"/>

### Short read mapping results (BAM files)
</a>
If you used the `-c` argument to provide coverage information, but you want to calculate `ratio_TSS` values for each isoform, SQANTI3 QC will still require you to supply short read mapping information. Therefore, we recommend users to supply the BAM files obtained as a result of mapping short reads to the genome via the `--SR_bam` option. Otherwise, you may supply short reads as raw data to SQANTI3 and have all pre-processing steps done internally, however, keep in mind that this will take longer as well as require more memory/computational resources.


<a name="cage"/>

## Incorporating CAGE peak data (`--CAGE_peak`)
</a>
To perform quality control of the **Transcription Start Site (TSS)**, users may provide CAGE peak data. This is particularly relevant given that 5' RNA degradation can generate variability in the TSS that can be mistaken for a novel TSS.

CAGE peak data must be provided as a BED file, where:

- Column 1 contains the **chromosome** (chr1, chr2...).
- Column 2 contains the zero-based index **start** coordinate.
- Column 3 contains one-based index **end** coordinate.
- Column 6 contains **strand** information.

The BED file needs to be supplied to `sqanti3_qc.py` via the `--CAGE_peak` argument.

Please, note that **public CAGE peak data** for human and mouse from the [refTSS](http://reftss.clst.riken.jp/reftss/Main_Page) database is supplied in SQANTI3's `data` folder.

SQANTI3 QC will use this data to validate the transcript's TSS. First, it will evaluate whether the TSS of the different long read-defined transcripts fall inside a CAGE peak (`within_CAGE_peak` column). Then, it will calculate the distance to the nearest detected CAGE peak (`dist_to_CAGE_peak` column in the output `*_classification.txt` file). A detailed glossary of classification file columns can be found [here](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC).


<a name="polya"/>

## Incorporating polyA information
</a>
To perform quality control of the **Transcription Termination Site (TTS)**, users may provide polyA motif or polyA site/peak information. SQANTI3 QC will use this data to create several descriptors regarding the presence of the TTS within polyA sites/motifs, as well as the distance to them, if detected.


<a name="polyamotif"/>

### PolyA motif data (`--polyA_motif_list`)
</a>
The polyA motif list should be supplied in the form of a text file. The most common polyA motifs for human and mouse are supplied in the SQANTI3 `data` folder (see [here](https://github.com/ConesaLab/SQANTI3/blob/master/data/polyA_motifs/mouse_and_human.polyA_motif.txt)).

The list of polyA motifs with the `--polyA_motif_list` parameter in `sqanti3_qc.py`.

If a polyA motif is identified at the 3' end of the transcript, the `polyA_motif_found` column in the classification file will be `TRUE`. In addition, SQANTI3 returns the sequence (`polyA_motif`) and the distance (`polyA_dist`) to the closest polyA motif detected. A detailed glossary of classification file columns can be found [here](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC).


<a name="polyasite"/>

### PolyA site data (`--polyA_peak`)
</a>
Complementary to polyA motif information, polyA site data can be supplied. It must be supplied as a BED file, where:

- Column 1 contains the **chromosome** (chr1, chr2...).
- Column 2 contains the zero-based index **start** coordinate.
- Column 3 contains one-based index **end** coordinate.
- Column 6 contains **strand** information.

For human, mouse, and worm, you can download **public polyA site data** from the [PolyASite atlas](https://polyasite.unibas.ch/atlas#2). Note, however, that the chromosomes in the BED files are missing the "chr" prefix. You will need to modify the downloaded BED file as follows:

```
sed -i 's/^/chr/' atlas.clusters.2.0.GRCh38.96.bed
```
 
The polyA site BED file must be supplied to `sqanti3_qc.py` via the `--polyA_peak` argument.

Using this information, SQANTI3 will write out the distance to the closest polyA site/peak (`dist_to_polyA_site` column) and whether the transcript's TTS was found inside a polyA site/peak (`within_polyA_site` column) to the output `*_classification.txt` file. A detailed glossary of classification file columns can be found [here](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC).


<a name="flcount"/>

## Supplying single or multi-sample full-length (FL) counts (`--fl_count`)
</a> 
`sqanti3_qc.py` supports single or multi-sample FL counts from PacBio Iso-Seq pipeline. There are three acceptable formats.

### Single Sample FL Count

A single sample FL Count file is automatically produced by the Iso-Seq With Mapping pipeline in [SMRTLink/SMRTAnalysis](https://www.pacb.com/products-and-services/analytical-software/) with the following format:

```
#
# -----------------
# Field explanation
# -----------------
# count_fl: Number of associated FL reads
# norm_fl: count_fl / total number of FL reads, mapped or unmapped
# Total Number of FL reads: 1065
#
pbid    count_fl        norm_fl
PB.1.1  2       1.8779e-03
PB.1.2  6       5.6338e-03
PB.1.3  3       2.8169e-03
PB.2.1  3       2.8169e-03
PB.3.1  2       1.8779e-03
PB.4.1  8       7.5117e-03

```


### Multi Sample Chained FL Count

A multi-sample FL Count file produced by the [chain_samples.py](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#chain) script in Cupcake will have the following format:

| superPBID | sample1 | sample2 | 
| --------- | ------- | ------- |
| PB.1.2 | 3 | NA |
| PB.2.1 | 2 | NA |
| PB.3.1 | 2 | 2 |
| PB.3.2 | 5 | 3 |
| PB.3.3 | 5 | 2 | 

This is a tab-delimited file.


### Multi Sample Demux FL Count

A multi-sample Demux FL Count file produced by the [demux_isoseq_with_genome.py](https://github.com/Magdoll/cDNA_Cupcake/wiki/Tutorial:-Demultiplexing-SMRT-Link-Iso-Seq-Jobs) script in Cupcake will have the following format:

```
id,sample1,sample2
PB.1.1,3,10
PB.1.2,0,11
PB.1.3,4,4
```

This is a comma-delimited file.


### FL count information in the SQANTI3 QC output

For each sample provided through the `--fl_count` option, `sqanti3_qc.py` will create a column in the `*_classification.txt` output file that is `FL.<sample>`. Note that this is the raw FL count provided. The sum of all the FL reads accross the samples associated to one transcript will be recorded in th `FL` column of the `*_classification.txt` output file.

When plotted, the script [SQANTI3_report.R](ConesaLab/SQANTI3/blob/tree/master/utilities/SQANTI3_report.R) will convert the FL counts to TPM using the formula:

```
                           raw FL count for PB.X.Y in sample1
FL_TPM(PB.X.Y,sample1) = ------------------------------------- x 10^6
                               total FL count in sample1
```

Two additional columns, `FL_TPM.<sample>` and `FL_TPM.<sample>_log10` will be added and output to a new classification file with the suffix `*_classification_TPM.txt`. Please, do not mix up `*_classification_TPM.txt`  and `*_classification.txt` files. The one used in subsequent scripts (filtering, future isoAnnot, etc.) will be the `_classification.txt` one. A detailed glossary of classification file columns can be found [here](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC).


<a name="tssratio"/>

## Selecting the metric used to aggregate TSS ratio across samples
</a>
As of v5.2, SQANTI3 includes the `--ratio_TSS_metric` argument, which can be used to tweak the results in the `TSS_ratio` column of the
classification.txt file. Briefly, SQANTI3 calculates the TSS ratio metric for all supplied short-read samples and replicates, and then summarizes the results to provide a single TSS ratio value.

The following options are available:
- `max` (default): uses the maximum value among the short-read samples and replicates supplied. 
- `mean`: mean value of samples and replicates is used. 
- `median`: median value of samples and replicates is used.
- `3rd quartile`: the 3rd quartile value of samplesl and replicates is used.

When aiming to discover novel and/or rare TSS, this metric will ensure that underrepresented sites will not be disregarded just because they are only captured in one sample. `median` and `3rd quartile`, on the other hand, allow the user to enforce different levels of robustness, preventing the metric from being driven by outliers and favoring widely-detected TSS to be preserved after QC. `mean` values, on the other hand, provide a balance between both scenarios in situations in which the degree of TSS novelty of the sample is unknown.


<a name="gff3"/>

## IsoAnnotLite and SQANTI3: using an existing tappAS annotation file (`--gff3`)
</a>
When the `--isoAnnotLite` flag is supplied to SQANTI3 QC, the tool will run [IsoAnnotLite](https://isoannot.tappas.org/isoannot-lite/) internally. IsoAnnotLite is a Python script designed to use the SQANTI3 QC output to generate the input for [tappAS](https://app.tappas.org/), i.e. tappAS-compatible  GFF3 file (see "Annotation features file format" section in the [tappAS overview](https://app.tappas.org/overview/)).

If you want IsoAnnotLite to perform **positional transfer of functional features**, you will also need to provide a pre-annotated tappAS GFF3 file via the `--gff3` argument. By doing this, functional feature information will be added to your long read-defined transcriptome, unlocking tappAS functional and feature-level analyses. You can find all functional annotation files available for tappAS [here](http://app.tappas.org/resources/downloads/gffs/). 

Note that, if **no pre-annotated tappAS GFF3 file is provided** (e.g. because there is no file available for your species), IsoAnnotLite will only perform structural information transfer, but a tappAS-like GFF3 file containing this structural information will still be generated. Therefore, in this case, users will still be able to load their data into tappAS for isoform-level expression analyses.

**⚠️Warning:** we are aware that some of the tappAS GFF3 files correspond to old transcriptome releases of reference databases (ENSEMBL, RefSeq, etc.). We are currently working on updating these annotations to later transcriptome versions. 


<a name="parallelization"/>

## Parallelization
</a>
There are two options related to parallelization:

- The `-t` (`--cpus`) parameter designates the number of CPUs used by the aligners for long and short reads. If your have supplied your transcriptome a GTF file and you do not provide short-read FASTA/FASTQ files, the `-t` option has no effect.
- The `-n` (`--chunks`) parameter divides the input transcriptome (GTF or FASTA/FASTQ) into chunks and runs SQANTI3 in parallel for each fragment, then combining the results into one classification and junctions file. 

Note that both of this parameters can be combined, for instance, if you have set `-t 30 -n 10`, then each chunk will get (30/10=3) CPUs.

