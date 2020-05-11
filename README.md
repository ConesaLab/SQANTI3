# SQANTI2
# This is the SQANTI2 version as it is right now. New commits and features will be added soon until getting SQANTI3


Last Updated: 02/07/2020 (v7.3.2)   Now works with Python 3.7 exclusively!!

## What is SQANTI2

SQANTI2 is a developer's version of SQANTI ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/), [code repository](https://bitbucket.org/ConesaLab/sqanti)). 


New features implemented in SQANTI2 not available in SQANTI:

* NMD detection -- new field `predicted_NMD` in classification output.
* Intron retention --- marked with `intron_retention` in `subcategory` field in classification output.
* CAGE peak --- new fields `dist_peak` and `within_peak` in classification output. Must provide CAGE peak data.
* polyA motif --- new field `polyA_motif` in classification output. Must provide polyA motif list.
* CDS-annotated GFF --- SQANTI2 outputs a `xxxx.cds.gff` GFF file that annotates CDS regions.
* PacBio Iso-Seq FL count multi-sample plotting --- use the `--fl_count` option for single or multi-sample FL counts


## SQANTI2 HowTos

* <a href="#install">Setting up SQANTI2</a>
* <a href="#input">Running SQANTI2 Classification</a>
   * <a href="#flcount">Single or Multi-Sample FL Count Plotting</a>
   * <a href="#exp">Short Read Expression</a>
* <a href="#filter">Filtering Isoforms using SQANTI2</a>
* <a href="#explain">SQANTI2 Output Explanation</a>
   * <a href="#class">Classification Output Explanation</a>
   * <a href="#junction">Junction Output Explanation</a>
   

![sqanti2workflow](https://github.com/Magdoll/images_public/blob/master/SQANTI2_figures/sqanti2_workflow.png)

## Updates

2020.02.07 updated to version 7.3.2. fixed default param for intrapriming to 0.6.

2020.02.06 updated to version 7.3.1. `sqanti_filter2.py` supports runA length filtering (`--runAlength`), monoexonic (`--filter_mono_exonic`) and skipping GTF/Isoform output (`--skipGTF` and `--skipFaFq`).

2020.02.04 updated to version 7.2.0. Fixed classification bug where multi-exon transcripts were listed as having exon count of 1.

2020.02.04 updated to version 7.1.0. Fixed classification `_TPM.txt` output format bug. Added `seq_A_downstream_TTS` field to classification output.

2020.01.31 updated to version 7.0.0. Minor re-classification changes for edge case mono-exonic transcripts. Supporting multithreading with `-n` option!

2019.12.12 updated to version 6.0.0. MAJOR re-classification changes! See [notes](https://www.dropbox.com/s/d1ya562ib00qhzc/Dec2019_SQANTI2_6.0.0_release_Notes.pptx?dl=0)

2019.11.24 updated to version 5.1.0. `--is_fusion` must be used together with `--gtf`. Added addtional tables if `--fl_count` supplied.

2019.11.12 updated to version 5.0.0. Now works exclusively with Python 3.7!

2019.09.24 updated to version 4.1. `sqanti_qc2.py` now supports `--fl_count` multi-sample with plotting.

2019.08.22 updated to version 4.0. `sqanti_filter2.py` now outputs GTF.

<details>
   <summary>Click here to see older update logs.</summary>
    
    2019.08.19 updated to version 3.9. Fixed minor bug in removing superPBID (not always the case) in `sqanti_qc2.py`
    
    2019.08.13 updated to version 3.8. Fixed `SQANTI_report2.R` printing bug for when polyA info is NA.
    
    2019.08.01 updated to version 3.7. `--expression` supports RSEM and Kallisto output. Expression and polyA motifs added to PDF plots.
    
    2019.07.26 updated to version 3.6. `--fl_count` is supported again!
    
    2019.07.25 updated to version 3.5. Checks for Cupcake version (8.1+)
    
    2019.07.24 updated to versoin 3.4. Now `--gtf` input option works with collapsed GFF format.
    
    2019.07.23 updated to version 3.3. Added CDS for GFF support and IR fix.
    
    2019.07.19 updated to version 3.2. `sqanti_qc2.py` fusion mode surpressed ORF prediction for the meantime. Minor mod to gmst to work on small input.
    
    2019.07.17 updated to version 3.1. `sqanti_qc2.py` now supports GTF input `--gtf` again. Minor change to NIC subtype categorization naming.
    
    2019.07.16 updated to version 3.0. now use Bioconda install of `gtfToGenePred` and `gffread`.
    
    2019.07.12 updated to version 2.9. `sqanti_qc2.py` now annotates NMD prediction.
    
    2019.06.19 updated to version 2.8. `sqanti_qc2.py` now works with fusion transcripts (must have ID `PBfusion.X`) using the `--is_fusion` option
    
    2019.06.02 updated to version 2.7. Fixed GMAP option bug + added distance to closest annotated start/end for the gene (not ref isoform) and filtering afterwards.
    
    2019.05.02 updated to version 2.6. Added deSALT aligner support using `--aligner=deSALT`. 
    
    2019.03.18 minor typo fixed for version 2.5. updated doc for `sqanti_filter2.py`
    
    2019.03.10 updated to version 2.5. Fixed `sqanti_filter2.py` missing fusion category also using polyA_motif as part of filtering.
    
    2019.02.27 updated to version 2.4.  Added polyA motif finding.
    
    2019.02.27 updated to version 2.3. `junction_category` fixed to check for (ss5,ss3) pairs in provided GTF.
    
    2019.02.26 updated to version 2.2. added support for CAGE peak (FANTOM5) and Intropolis junction BED. 
    
    2018.10.15 updated to version 1.1. modified use of SAM to GFF with added `source` parameter.

</details>



<a name="install"/>

## Prerequisite

* Perl
* Minimap2 

### Python-related libraries

* Python (3.7)
* pysam
* psutil
* bx-python
* BioPython
* BCBioGFF
* [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#install)

### R-related libraries

* R (>= 3.4.0)
* R packages (for `sqanti_qc2.py`): ggplot2, scales, reshape, gridExtra, grid, dplyr

## Installing Python dependencies

I recommend using Anaconda which makes installing all the Python packages much easier. If you already have Anaconda installed because you use [Iso-Seq3](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda) or [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step), you can activate your current environment and directly go to step (4).

(1)  Here's the generic Anaconda installation for [Linux environment](http://docs.continuum.io/anaconda/install/#linux-install). Currently only Linux environment supported.

```
export PATH=$HOME/anacondaPy37/bin:$PATH
conda -V
conda update conda
```

(2) Create a virutal environment. I will call it `anaCogent3`. Type `y` to agree to the interactive questions.

```
conda create -n anaCogent3 python=3.7 anaconda
source activate anaCogent3
```

(3) Once you have activated the virtualenv, you should see your prompt changing to something like this:

```
(anaCogent3)-bash-4.1$
```

(4) Install additional required libraries:

```
conda install -n anaCogent3 -c bioconda pysam
conda install -n anaCogent3 psutil
conda install -n anaCogent3 biopython
conda install -n anaCogent3 -c bioconda bx-python
conda install -n anaCogent3 -c bioconda bcbiogff
conda install -n anaCogent3 -c bioconda gffread
```

We also need to install [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html) that seems to have some issues with Python 3.7 (or openssl). At this point, the easiest solution is to download it from [UCSC Download Page](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) and add the binary to your `$PATH` variable:

```
export PATH=$PATH:<path_to>/gtfToGenePred
```


If you don't already have [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#install) installed, you can do that now:

```
$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
$ cd cDNA_Cupcake
$ python setup.py build
$ python setup.py install
```

No installation for SQANTI2 itself is required. The scripts can be run directly.


## Running SQANTI2

Activate the Anaconda environment. Make sure minimap2 works. Add `cDNA_Cupcake/sequence` to `$PYTHONPATH`.

```
$ source activate anaCogent3
(anaCogent3)-bash-4.1$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(anaCogent3)-bash-4.1$ minimap2 --version
2.15-r905
```

<a name="input"/>

#### Input to SQANTI2 Classification

* *Iso-Seq output*. Preferably already mapped to the genome and [collapsed to unique transcripts](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse). (FASTA/FASTQ/GTF)
* *Reference annotation* in GTF format. For example [GENCODE](https://www.gencodegenes.org/releases/current.html) or [CHESS](http://ccb.jhu.edu/chess/).
* *Reference genome*, in FASTA format. For example hg38. *Make sure your annotation GTF is based on the correct ref genome version!*

Optionally:

* CAGE Peak data (from FANTOM5). I've provided a version of [CAGE Peak for hg38 genome](https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz) which was originally from [FANTOM5](http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/). 

* [Intropolis](https://github.com/nellore/intropolis/blob/master/README.md) Junction BED file. I've provided a version of [Intropolis for hg38 genome, modified into STAR junction format](https://github.com/Magdoll/images_public/tree/master/SQANTI2_support_data).

* polyA motif list. A ranked  list of polyA motifs to find upstream of the 3' end site. See [human polyA list](https://raw.githubusercontent.com/Magdoll/images_public/master/SQANTI2_support_data/human.polyA.list.txt) for an example.

* FL count information. See <a href="#flcount">FL count section</a> to include Iso-Seq FL count information for each isoform.

* Short read expression. See <a href="#exp">Short Read Expression section</a> to include short read expression (RSEM or Kallisto).

### Running SQANTI2 Classification

The script usage is:

```
python sqanti_qc2.py [-t cpus] [-n chunks]
     [--gtf] [--skipORF] 
     [-c shortread_STAR_junction_out] 
     [--cage_peak CAGE_PEAK_BED]
     [--polyA_motif_list POLYA_LIST]
     [--fl_count FL_COUNT]
     [--expression EXPRESSION]
     [--aligner_choice=minimap2,deSALT]
     [--is_fusion]
     <input_fasta> <annotation_gtf> <genome_fasta>
```

If you don't feel like running the ORF prediction part, use `--skipORF`. Just know that all your transcripts will be annotated as non-coding.
If you have short read data, you can run STAR to get the junction file (usually called `SJ.out.tab`, see [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)) and supply it to SQANTI2.

If `--aligner_choice=minimap2`, the minimap2 parameter used currently is: `minimap2 -ax splice --secondary=no -C5 -O6,24 -B4 -uf`
If `--aligner_choice=deSALT`, the deSALT parameter used currently is: `deSALT aln -x ccs`. 

You can look at the [`MINIMAP2_CMD` and `DESALT_CMD` in `sqanti_qc2.py` for the full command format](https://github.com/Magdoll/SQANTI2/blob/master/sqanti_qc2.py#L61).


There are two options related to parallelization. The first is `-t` (`--cpus`) that designates the number of CPUs used by the aliger. 
If your input is GTF (using `--gtf` option), the `-t` option has no effect.
The second is `-n` (`--chunks`) that chunks the input (GTF or fasta) into chunks and run SQANTI2 in parallel before combining them. 
Note that if you have `-t 30 -n 10`, then each chunk gets (30/10=3) CPUs.

For example:

```
python sqanti_qc2.py -t 30 -n 10 --gtf example/test.gtf \
     gencode.v29.annotation.gtf hg38.fa \
     --fl_count example/test.chained_count.txt \
     --polyA_motif_list example/polyA.list \
     --cage_peak hg38.cage_peak_phase1and2combined_coord.bed \
     -c "Public_Intronpolis/*10.modified"
```


For fusion transcripts, you must use the `--is_fusion` option for `sqanti_qc2.py` to work properly. Furthermore, the IDs in the input FASTA/FASTQ *must* have the format `PBfusion.X`, as is output by [`fusion_finder.py` in Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#fusion).


### SQANTI2 Classification Output

You can look at the [example](https://github.com/Magdoll/SQANTI2/tree/master/example) subfolder for a sample output. The PDF file shows all the figures drawn using R script [SQANTI_report2.R](https://github.com/Magdoll/SQANTI2/blob/master/utilities/SQANTI_report2.R), taking the `_classification.txt` and `_junctions.txt` as the two input. If you know R well, you are free to modify the R script to add new figures! I will be constantly adding new figures as well.

Detailed explanation of `_classification.txt` and `_junctions.txt` <a href="#explain">below</a>.


<a name="flcount"/>

### Single or Multi-Sample FL Count Plotting

Supported since: version 4.1

`sqanti_qc2.py` supports single or multi-sample FL counts from PacBio Iso-Seq pipeline. There are three acceptable formats.

#### Single Sample FL Count

A single sample FL Count file is automatically produced by the Iso-Seq With Mapping pipeline in [SMRTLink/SMRTAnalysis](https://www.pacb.com/products-and-services/analytical-software/) with the following format:

#### Multi Sample Chained FL Count

A multi-sample FL Count file produced by the [chain_samples.py](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#chain) script in Cupcake will have the following format:

| superPBID | sample1 | sample2 | 
| --------- | ------- | ------- |
| PB.1.2 | 3 | NA |
| PB.2.1 | 2 | NA |
| PB.3.1 | 2 | 2 |
| PB.3.2 | 5 | 3 |
| PB.3.3 | 5 | 2 | 

This is a tab-delimited file.


#### Multi Sample Demux FL Count

A multi-sample Demux FL Count file produced by the [demux_isoseq_with_genome.py](https://github.com/Magdoll/cDNA_Cupcake/wiki/Tutorial:-Demultiplexing-SMRT-Link-Iso-Seq-Jobs) script in Cupcake will have the following format:

```
id,sample1,sample2
PB.1.1,3,10
PB.1.2,0,11
PB.1.3,4,4
```

This is a comma-delimited file.


#### SQANTI2 output of FL Count information

For each sample provided through the `--fl_count` option, `sqanti_qc2.py` will create a column in the `.classification.txt` output file that is `FL.<sample>`. Note that this is the raw FL count provided.

When plotted, the script [SQANTI_report2.R](https://github.com/Magdoll/SQANTI2/blob/master/utilities/SQANTI_report2.R) will convert the FL counts to TPM using the formula:

```
                           raw FL count for PB.X.Y in sample1
FL_TPM(PB.X.Y,sample1) = ------------------------------------- x 10^6
                               total FL count in sample1
```

Two additional columns, `FL_TPM.<sample>` and `FL_TPM.<sample>_log10` will be added and output to a new classification file with the suffix `.classification_TPM.txt`.


<a name="exp"/>

### Short Read Expression

Use `--expression` to optionally provide short read expression data. Two formats are supported.

#### Kallisto Expression Input

Kallisto expression files have the format:

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

#### RSEM Expression Input

RSEM expression files have the format:

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

<a name="filter"/>

### Filtering Isoforms using SQANTI2


I've made a lightweight filtering script based on SQANTI2 output that filters for two things: (a) intra-priming and (b) short read junction support.  

The script usage is:

```
usage: sqanti_filter2.py sqanti_class isoforms gtf_file
                         [--sam SAM] [--faa FAA] 
                         [-r RUNALENGTH] [-a INTRAPRIMING]
                         [-c MIN_COV] [-m MAX_DIST_TO_KNOWN_END]
                         [--filter_mono_exonic] 
                         [--skipGTF] [--skipFaFq] 
                         sqanti_class isoforms gtf_file
sqanti_filter2.py: error: the following arguments are required: sqanti_class, isoforms, gtf_file

python sqanti_filter2.py [classification] [fasta] [sam] [gtf]
         [-a INTRAPRIMING] [-c MIN_COV] [-m MAX_DIST_TO_KNOWN_END]
```

where `-a` determines the fraction of genomic 'A's above which the isoform will be filtered. The default is `-a 0.6`. 
`-r` is another option for looking at genomic 'A's that looks at the immediate run-A length. The default is `-r 6`. 

`-m` sets the maximum distance to an annotated 3' end (the `diff_to_gene_TTS` field in classification output) to offset the intrapriming rule.

`-c` is the filter for the minimum short read junction support (looking at the `min_cov` field in `.classification.txt`), and can only be used if you have short read data.


For example:

```
python sqanti_filter2.py test_classification.txt \
                         test.renamed_corrected.fasta \
                         test.gtf
```

The current filtering rules are as follow:

* If a transcript is FSM, then it is kept unless the 3' end is unreliable (intrapriming).
* If a transcript is not FSM, then it is kept only if all of below are true: (1) 3' end is reliable (2) does not have a junction that is labeled as RTSwitching (3) all junctions are either canonical or has short read coverage above `-c` threshold


<a name="explain"/>

### SQANTI2 Output Explanation


SQANTI/SQANTI2 categorizes each isoform by finding the best matching reference transcript in the following order:

* FSM (*Full Splice Match*): meaning the reference and query isoform have the same number of exons and each internal junction agree. The exact 5' start and 3' end can differ by any amount.

* ISM (*Incomplete Splice Match*): the query isoform has fewer 5' exons than the reference, but each internal junction agree. The exact 5' start and 3' end can differ by any amount.

* NIC (*Novel In Catalog*): the query isoform does not have a FSM or ISM match, but is using a combination of known donor/acceptor sites.

* NNC (*Novel Not in Catalog*): the query isoform does not have a FSM or ISM match, and has at least one donor or acceptor site that is not annotated.

* *Antisense*: the query isoform does not have overlap a same-strand reference gene but is anti-sense to an annotated gene. 

* *Genic Intron*: the query isoform is completely contained within an annotated intron.

* *Genic Genomic*: the query isoform overlaps with introns and exons.

* *Intergenic*: the query isoform is in the intergenic region.

![sqanti_explain](https://github.com/Magdoll/images_public/blob/master/SQANTI2_figures/sqanti2_classification.png)


Some of the classifications have further subtypes (the `subtype`) field in SQANTI2 classification output. They are explained below.

![ISM_subtype](https://github.com/Magdoll/images_public/blob/master/SQANTI2_figures/sqanti2_ISM_subtype.png)

Novel isoforms are subtyped based on whether they use a combination of known junctions (junctions are pairs of <donor>,<acceptor> sites), a combination of known splice sites (the individual donor and acceptor sites are known, but at least combination is novel), or at least one splice site (donor or acceptor) is novel. 

![NIC_subtype](https://github.com/Magdoll/images_public/blob/master/SQANTI2_figures/sqanti2_NIC_subtype.png)

<a name="class"/>

#### Classification Output Explanation

The output `.classification.txt` has the following fields:

1. `isoform`: the isoform ID. Usually in `PB.X.Y` format.
2. `chrom`: chromosome.
3. `strand`: strand.
4. `length`: isoform length.
5. `exons`: number of exons.
6. `structural_category`: one of the categories ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron"]
7. `associated_gene`: the reference gene name.
8. `associated_transcript`: the reference transcript name.
9. `ref_length`: reference transcript length.
10. `ref_exons`: reference transcript number of exons.
11. `diff_to_TSS`: distance of query isoform 5' start to reference transcript start end. Negative value means query starts downstream of reference.
12. `diff_to_TTS`: distance of query isoform 3' end to reference annotated end site. Negative value means query ends upstream of reference.
13. `diff_to_gene_TSS`: distance of query isoform 5' start to the closest start end of *any* transcripts of the matching gene. This field is different from `diff_to_TSS` since it's looking at all annotated starts of a gene. Negative value means query starts downstream of reference.
14. `diff_to_gene_TTS`: distance of query isoform 3' end to the closest end of *any* transcripts of the matching gene. Negative value means query ends upstream of reference.
13. `subcategory`: additional splicing categorization, separated by semi-colons. Categories include: `mono-exon`, `multi-exon`. Intron rentention is marked with `intron_retention`. 
14. `RTS_stage`: TRUE if one of the junctions could be a RT switching artifact.
15. `all_canonical`: TRUE if all junctions have canonical splice sites.
16. `min_sample_cov`: sample with minimum coverage.
17. `min_cov`: minimum junction coverage based on short read STAR junction output file. NA if no short read given.
18. `min_cov_pos`: the junction that had the fewest coverage. NA if no short read data given.
19. `sd_cov`: standard deviation of junction coverage counts from short read data. NA if no short read data given.
20. `FL` or `FL.<sample>`: FL count associated with this isoform per sample if `--fl_count` is provided, otherwise NA.
21. `n_indels`: total number of indels based on alignment.
22. `n_indels_junc`: number of junctions in this isoform that have alignment indels near the junction site (indicating potentially unreliable junctions).
23. `bite`: TRUE if all junctions match reference junctions completely.
24. `iso_exp`: short read expression for this isoform if `--expression` is provided, otherwise NA.
25. `gene_exp`: short read expression for the gene associated with this isoform (summing over all isoforms) if `--expression` is provided, otherwise NA.
26. `ratio_exp`: ratio of `iso_exp` to `gene_exp` if `--expression` is provided, otherwise NA.
27. `FSM_class`: ignore this field for now.
28. `ORF_length`: predicted ORF length.
29. `CDS_length`: predicted CDS length. 
30. `CDS_start`: CDS start.
31. `CDS_end`: CDS end.
32. `CDS_genomic_start`: genomic coordinate of the CDS start. If on - strand, this coord will be greater than the end.
33. `CDS_genomic_end`: genomic coordinate of the CDS end. If on - strand, this coord will be smaller than the start.
34. `predicted_NMD`: TRUE if there's a predicted ORF and CDS ends before the last junction; FALSE if otherwise. NA if non-coding.
35. `perc_A_downstreamTTS`: percent of genomic "A"s in the downstream 20 bp window. If this number if high (say > 0.8), the 3' end site of this isoform is probably not reliable.
36. `seq_A_downstream_TTS`: sequence of the downstream 20 bp window.
37. `dist_peak`: distance to closest TSS based on CAGE Peak data. Negative means upstream of TSS and positive means downstream of TSS. Strand-specific. SQANTI2 only searches for nearby CAGE Peaks within 10000 bp of the PacBio transcript start site. Will be `NA` if none are found within 10000 bp.
38. `within_peak`: TRUE if the PacBio transcript start site is within a CAGE Peak. 
39. `polyA_motif`: if `--polyA_motif_list` is given, shows the top ranking polyA motif found within 50 bp upstream of end.
40. `polyA_dist`: if `--polyA_motif_list` is given, shows the location of the  last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream. 


<a name="junction"/>

### Junction Output Explanation


THe `.junctions.txt` file shows every junction for every PB isoform. NOTE because of this the *same* junction might appear multiple times if they are shared by multiple PB isoforms. 

1. `isoform`: Isoform ID
2. `junction_number`: The i-th junction of the isoform
3. `chrom`: Chromosome 
4. `strand`: Strand
5. `genomic_start_coord`: Start of the junction (1-based), note that if on - strand, this would be the junction acceptor site instead.
6. `genomic_end_coord`: End of the junction (1-based), note that if on - strand, this would be the junction donor site instead.
7. `transcript_coord`: Currently not implemented. Ignore.
8. `junction_category`: `known` if the (donor-acceptor) combination is annotated in the GTF file, `novel` otherwise. Note that it is possible to have a `novel` junction even though both the donor and acceptor site are known, since the combination might be novel.
9. `start_site_category`: `known` if the junction start site is annotated. If on - strand, this is actually the donor site.
10. `end_site_category`: `known` if the junction end site is annotated. If on - strand, this is actually the acceptor site.
11. `diff_to_Ref_start_site`: distance to closest annotated junction start site. If on - strand, this is actually the donor site.
12. `diff_to_Ref_end_site`: distance to closest annotated junction end site. If on - strand, this is actually the acceptor site.
13. `bite_junction`: TRUE if either or both the junction start/end site matches annotation.
14. `splice_site`: Splice motif.
15. `RTS_junction`: TRUE if junction is predicted to a template switching artifact.
16. `indel_near_junct`: TRUE if there is alignment indel error near the junction site, indicating potential junction incorrectness.
17. `sample_with_cov`: If `--coverage` (short read junction coverage info) is provided, shows the number of samples (cov files) that have short read that support this junction.
18. `total_coverage`: Total number of short read support from all samples that cover this junction.

