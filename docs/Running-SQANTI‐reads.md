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

Only **sampleID** and **file_acc** are required. `sampleID` is the label used for the sample everywhere in the output; `file_acc` links the row to that sample's SQANTI3-QC output on disk (see the directory layout below).

**Adding experimental factors (for faceted plots).** The design file can carry any number of extra columns describing your samples — e.g. `platform`, `condition`, `age`, `sex`, `tissue`, `RIN_group`. Give one of those column names to `--factor` and every multi-sample plot is split by that factor (one panel/colour per level), and the replicate-concordance and transcript-divergency comparisons are computed within and between its groups. For example:

| sampleID | file_acc | platform |
| --------- | ------- | -------- |
| wtc11_PBcDNA | ENCFF003QZT | PacBio |
| wtc11_ONTR2C2 | ENCFF063ASB | ONT |
| wtc11_ONTdRNA | ENCFF104BNW | ONT |

then run with `--factor platform`. The extra columns are optional and only used when you name one via `--factor`; if you omit `--factor`, all samples are shown together.

To run SQANTI-reads in _fast_ mode, you must first, run SQANTI3-QC on all your samples individually (great opportunity for you to parallelize this using your system's slurm, gnu parallel or your favourite parallelization method). Using the parameter ```-d/--sqanti_dirs```, you can give sqanti3_reads.py the path to the parent directory where all the SQANTI3-QC output directories are stored.

We recommend using _fast_ mode, this ensures you can easily parallelize your work and have more control over mapping parameters. The usual reads pre-processing pipelines for ONT and PacBio to run SQANTI-reads include:  
- For ONT: Running pychopper to strand the reads, remove primers and polyA tails.  
- For PacBio: Running lima and refine from IsoSeq pipeline, and ```bamtools convert -format fastq -in fl.bam``` to transform PacBio FL bams into fastq files for mapping.  
- Common next steps: mapping with minimap2, transforming to gtf using spliced_bam2gff and running SQANTI3-QC.

**Directory layout SQANTI-reads expects.** Under the `--sqanti_dirs` folder, each sample's SQANTI3-QC output must live in a sub-directory named exactly as its `file_acc`, and the files inside must be prefixed with that sample's `sampleID`. For the first design row (`wtc11_PBcDNA` / `ENCFF003QZT`):

```
<sqanti_dirs>/                            # the folder you pass to --sqanti_dirs
└── ENCFF003QZT/                          # = file_acc  (one sub-dir per sample)
    ├── wtc11_PBcDNA_classification.txt   # = {sampleID}_classification.txt
    ├── wtc11_PBcDNA_junctions.txt        # = {sampleID}_junctions.txt
    └── wtc11_PBcDNA_corrected.gtf        # = {sampleID}_corrected.gtf
```

You control these two names in the SQANTI3-QC run: `-d ./ENCFF003QZT` sets the sub-directory (= `file_acc`) and `-o wtc11_PBcDNA` sets the file prefix (= `sampleID`). Making these match your design file is the most common thing to get wrong — if a run reports missing classification/junction files, check that the sub-directory name and file prefix line up with the `file_acc` and `sampleID` columns.

Example run starting after running spliced_bam2gff (one QC call per sample, then one SQANTI-reads call):
```
python sqanti3_qc.py --isoforms ENCFF003QZT.gff --refGTF hg38.ensGene.gtf --refFasta hg38.fa --min_ref_len 0 --aligner_choice minimap2 -t 8 -d ./ENCFF003QZT -o wtc11_PBcDNA
# ...repeat sqanti3_qc.py for every sample...

# then aggregate all of them (add --factor <column> to split plots by a design factor):
python sqanti3_reads.py --design design.csv --refGTF hg38.ensGene.gtf --sqanti_dirs ./ --report both
python sqanti3_reads.py --design design.csv --refGTF hg38.ensGene.gtf --sqanti_dirs ./ --factor platform --report both
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
python sqanti3_reads.py --design design.csv --refGTF hg38.ensGene.gtf --refFasta hg38.fa
```

Fastq files are named {file_acc}\*.fastq (e.g. ENCFF003QZT_PB.fastq) and are stored in the current directory.

As in _fast_ mode, you can add extra design columns (e.g. `platform`, `age`) and pass `--factor <column>` to facet the plots, and `--report {pdf,html,both}` to choose the report format.

## Getting ready

Before running SQANTI-reads, you will need to:

**Create the SQANTI3 conda environment** (once, from the YAML shipped in the repo — it creates an environment named `sqanti3`) **and activate it:**

```
(base)-bash-4.1$ conda env create -f SQANTI3.conda_env.yml   # first time only
(base)-bash-4.1$ conda activate sqanti3
(sqanti3)-bash-4.1$
```

## Arguments and parameters in SQANTI-reads

The SQANTI-reads script accepts the following arguments:

```
usage: sqanti3_reads.py [-h] [--refFasta REFFASTA] --refGTF REFGTF -de
                        INDESIGN [-i INPUT_DIR] [-p PREFIX] [-d SQANTI_DIRS]
                        [-o OUTPUT] [--report {pdf,html,both}]
                        [--config CONFIG] [--all_tables] [--pca_tables]
                        [--min_ref_len MIN_REF_LEN] [-ge ANNOTEXP]
                        [-je JXNEXP] [-pc PERCCOV] [-pj PERCMAXJXN]
                        [--aligner_choice {minimap2,uLTRA}] [-s SITES]
                        [--skip_hash] [-f INFACTOR] [-fl FACTORLVL]
                        [--skip_plots] [-t CPUS] [-n CHUNKS] [-j JOBS]
                        [--force_id_ignore] [--verbose] [-v]

Structural and Quality Annotation of Novel Transcript Isoforms

options:
  -h, --help            show this help message and exit
  --force_id_ignore     Allow the usage of transcript IDs non related with
                        PacBio's nomenclature (PB.X.Y)

Required arguments:
  --refFasta REFFASTA   Reference genome (Fasta format). Required unless
                        running in fast mode (--sqanti_dirs).
  --refGTF REFGTF       Reference annotation file (GTF format)
  -de INDESIGN, --design INDESIGN
                        Path to design file, must have sampleID and file_acc
                        column.

Input/Output options:
  -i INPUT_DIR, --raw_data_dir INPUT_DIR
                        Path to directory where fastq/GTF files are stored
                        (for running SQANTI3 from scratch). Default: Directory
                        where the script was run.
  -p PREFIX, --prefix PREFIX
                        SQANTI-reads output filename prefix. Default:
                        sqantiReads
  -d SQANTI_DIRS, --sqanti_dirs SQANTI_DIRS
                        Directory containing existing SQANTI3 output folders.
                        Use this to skip re-running QC and proceed directly to
                        aggregation.
  -o OUTPUT, --output OUTPUT
                        Directory for output sqanti_reads files (plots,
                        tables, design file). Default: Directory where the
                        script was run.
  --report {pdf,html,both}
                        Default: pdf
  --config CONFIG       YAML file overriding QC-flag thresholds and other cut-
                        offs (intra-priming %A, PCA variance, length bins).
                        Defaults reproduce the standard behaviour.
  --all_tables          Export all output tables. Default tables are gene
                        counts, ujc counts, length_summary, cv and
                        underannotated gene tables
  --pca_tables          Export table for making PCA plots

Filtering options:
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length. Default: 0 bp
  -ge ANNOTEXP, --gene_expression ANNOTEXP
                        Expression cut off level for determining
                        underannotated genes. Default = 100
  -je JXNEXP, --jxn_expression JXNEXP
                        Coverage threshold for detected reference donors and
                        acceptor. Default = 10
  -pc PERCCOV, --perc_coverage PERCCOV
                        Percent gene coverage of UJC for determining well-
                        covered unannotated transcripts. Default = 20
  -pj PERCMAXJXN, --perc_junctions PERCMAXJXN
                        Percent of the max junctions in gene for determining
                        near full-length putative novel transcripts. Default =
                        80

Analysis options:
  --aligner_choice {minimap2,uLTRA}
                        Default: minimap2
  -s SITES, --sites SITES
                        Set of splice sites to be considered as canonical
                        (comma-separated list of splice sites). Default:
                        GTAG,GCAG,ATAC.
  --skip_hash           Skip UJC hashing and read the already-hashed
                        reads_classification straight from --sqanti_dirs. Use
                        when the classification files already carry a jxnHash
                        column (avoids re-hashing).

Visualization options:
  -f INFACTOR, --factor INFACTOR
                        This is the column name that plots are to be faceted
                        by. Default: None
  -fl FACTORLVL, --factor_level FACTORLVL
                        Factor level to evaluate for underannotation
  --skip_plots          Skip the plotting step

Performance options:
  -t CPUS, --cpus CPUS  Number of threads used during alignment by aligners.
                        Default: 10
  -n CHUNKS, --chunks CHUNKS
                        Number of chunks to split SQANTI3 analysis in for
                        speed up. Default: 1
  -j JOBS, --jobs JOBS  Number of samples to process in parallel (UJC hashing,
                        and per-sample SQANTI3-QC in simple mode). Note: each
                        parallel QC job still uses --cpus alignment threads, so
                        peak cores ~= jobs*cpus. Default: 1 (serial)

Optional arguments:
  --verbose             If verbose is run, it will print all steps, by default
                        it is FALSE
  -v, --version         Display program version number.
```

## SQANTI-reads output

The sqanti-reads output is written to the path specified in ```--output``` or ```-o``` argument, also appending the prefix provided via the ```--prefix``` or -p argument.

The following output files are generated after running it:

### Modified input files

New reads_classification.txt files will be generated, including all columns from the original files plus two additional columns:
- jxn_string: includes the unique junction chain string as ```chromosome_strand_junction1_junction2_junctionN```. Monoexons are coded as ```chromosome_strand_monoexon_readID```.
- jxnHash: includes the hashing of the jxn_string for easier work in next steps.

The design file will be modified to include all initial columns (e.g. sampleID, file_acc, platform, age, sex) and two additional columns:

- classification_file: will include the full paths to the reads_classification.txt files (output of SQANTI-reads).
- junction_file: will include the full paths to the junctions.txt files (output of SQANTI3).

### Results default .csv files

- sqantiReads_cv.csv: Provides metrics on the coefficient of variation of reference junctions, per sample.
- sqantiReads_gene_counts.csv: Provides the number of reads in each structural category, per gene, per sample.
- sqantiReads_length_summary.csv: Provides the number and percentage of reads in length categories per sample.
- sqantiReads_ujc_counts.csv: Provides a list of junction strings in each sample and the number of reads in each sample associated with each junction string.
- sqantiReads_gene_classification.csv: For genes with coverage meeting a user defined threshold (-ge), provides the annotation category of each gene.
- sqantiReads_putative_novel_transcripts.csv: Provides metrics on NIC and NNC UJCs and flags putative novel transcripts.
- sqantiReads_jxn_offsets.csv: Splice-site fuzziness table. For each sample, the count of donor/acceptor site observations at each signed offset (in bp) from the nearest reference site, within the `jxn_offset_window`, split by canonical vs non-canonical. Feeds the offset-spectrum, precision-profile and canonical-split plots and the `perc_sites_imprecise` QC flag.
- sqantiReads_completeness.csv: 5'/3' read-end completeness table. For each sample and end (5prime/3prime), the count of reads at each absolute distance (in bp) from the annotated gene end. Feeds the completeness-profile plots and the `perc_5p_within_window` / `perc_3p_within_window` QC flags.
- sqantiReads_site_offsets.csv: Per reference splice site, each sample's median signed offset and observation count (reference coordinate reconstructed as `genomic_coord + diff_to_Ref`, donor/acceptor labelled strand-aware). Feeds the replicate splice-site-precision concordance view (F7).
- sqantiReads_fuzz_by_depth.csv: Per sample and junction read-depth bin (from `total_coverage_unique`), the number of splice-site observations and how many are imprecise (offset ≠ 0). Feeds the imprecision-vs-depth view (F4). Empty when the junction files carry no per-junction coverage.
- sqantiReads_support.csv: Per sample, the orthogonal-support scalars derived from SQANTI3's optional evidence columns — CAGE 5'-end support (`perc_within_cage`, `median_dist_cage`), polyA 3'-end support (`perc_within_polya`, `perc_polya_motif`, `median_dist_polya`), and short-read support (`median_ratio_TSS`, `perc_jxn_SR_supported`). Feeds the CAGE / polyA / short-read cohort views and their scorecard metrics. Every value is empty (NaN) for samples whose SQANTI3 run had no CAGE / polyA / short-read data, in which case the corresponding views simply do not appear.

### Results optional .csv files

- sqantiReads_cv_acc/don_counts.csv: Provides the number of detected annotated donors and acceptors in each junction variation category.
- sqantiReads_cv_acc/don_counts.csv: Provides the number of detected annotated donors and acceptors in each junction variation category.
- sqantiReads_FSM/ISM/NIC/NNC_counts.csv: Provides the number of reads in each sub-category for FSMs, ISMs, NICs and NNCs (one file each; NICs and NNCs are reported separately, and the NNC sub-category plot is coloured in reds to match the NNC structural-category colour).
- sqantiReads_err_counts.csv: Provides the number and percentage of reads with evidence of intrapriming, RT-switching and non-canonical junctions per sample.
- sqantiReads_pca_loadings.csv: Gives loadings of PC1 and PC2 for the PCA analysis.
- sqantiReads_pca_variance.csv: Gives the proportion variance explained by each PC.

### Results plots and reports

**At a glance**, the report is organised into the sections below (each is described in detail further down; a section is skipped when its input is absent):

| Report section | What it tells you | When it appears |
| --- | --- | --- |
| Summary & QC flags | Per-sample pass / warn / fail on the key QC metrics | HTML report |
| Sample-outlier scorecard | Which samples diverge from the cohort (robust z-scores) | ≥ 4 samples |
| Structural-category composition | Reads & UJCs per SQANTI3 category (FSM/ISM/NIC/NNC/…), per sample and per gene | Always |
| Genes & UJCs detected | How many genes / junction chains, split by read support | Always |
| Read length distribution | Length composition + full per-read distribution | Always |
| FSM/ISM/NIC/NNC sub-categories | Which sub-types (alt 5'/3' ends, fragments, intron retention, …) dominate | Always |
| Junction classification & artefacts | Known/novel × canonical junctions; RT-switching, intra-priming, novel non-canonical burden | Always |
| Splice-site detection & fuzziness | Donor/acceptor detection variation, offset spectrum, precision profile, tandem (NAGNAG), imprecision vs depth | Always (some parts need per-junction coverage) |
| 5'/3' completeness | How close read ends sit to the annotated gene boundaries | Always |
| UJC saturation, concordance, UpSet, overlap | Library saturation and how junction chains are shared across samples | Multi-sample |
| Sample similarity (PCA) | Overall clustering of samples by QC profile | ≥ 2 samples |
| Between-sample comparison | Composition drift (JS distance) and depth-normalised junction yield | Multi-sample |
| Replicate concordance (+ MAD) | How tightly replicates agree on composition / length / precision | `--factor` with ≥ 2 replicates |
| Transcript divergency | Per-gene TD = (NIC+NNC)/(FSM+NIC+NNC): per-sample, distribution, and between-group differences | TD differences need `--factor` |
| Orthogonal support (CAGE / polyA / short-read) | Independent confirmation of transcript ends & junctions | Only if that evidence was given to SQANTI3 |
| Under-annotation analysis | Genes likely missing from the annotation / novel-transcript candidates | Always |

The report format is controlled by `--report {pdf,html,both}` (default `pdf`):

- **`pdf`** — `sqantiReads_report.pdf`: single static report, one plot per page. Contains the QC metrics plots followed by an "Under-annotation analysis" section (gene-category bar chart plus per-category scatterplots of gene coverage vs. junctions).
- **`html`** — `sqantiReads_report.html`: a single **self-contained, interactive** report (Plotly, works offline — no browser plugins or internet needed). It opens with a **Summary & QC-flags** panel (per-sample metrics with pass/warn/fail badges) and then interactive versions of the QC figures (hover for values, zoom, toggle samples in the legend), ending with the under-annotation section. No PDF is written in this mode.
- **`both`** — writes both `sqantiReads_report.pdf` and `sqantiReads_report.html`.

When an HTML report is produced (`html` or `both`), a machine-readable **`sqantiReads_qc_summary.json`** is also written, with per-sample metrics and QC flags for pipeline integration. When the cohort has enough samples (`sample_scorecard.min_samples`, default 4), it also carries a **cohort-relative sample-outlier scorecard**: each read-QC metric is turned into a robust z-score against the cohort (median / MAD across samples), and a sample is flagged when it diverges from its peers on several metrics at once. The scorecard is dataset-agnostic — it uses no absolute cut-offs, so it stays quiet when all samples agree and is skipped (with a stated reason) below `min_samples`. Both reports render it as a heatmap (positive / red = worse than peers). The metrics it draws on span read length, structural composition (%ISM, ISM fragment fraction, composition drift), junction quality (novel-junction burden, novel non-canonical burden, splice-site imprecision, depth-normalised junction yield), read-end completeness (5'/3'), the artefact rates (RT-switching, intra-priming), and — when the corresponding evidence was supplied to SQANTI3 — the orthogonal-support fractions (`perc_within_cage`, `perc_within_polya`, `perc_jxn_SR_supported`); any metric whose source column was not computed is silently dropped.

Both reports also include UJC-level views: a **saturation / rarefaction** curve (expected unique junction chains vs. sequencing depth — a plateau indicates a saturated library), a **replicate-concordance** heatmap (per-UJC read-count correlation between samples), an **UpSet** plot of UJCs shared across samples (intersection and set-size bars stacked by structural category), and a **pairwise UJC overlap-index** heatmap — an asymmetric matrix whose entry (row A, column B) is `|A ∩ B| / |A|`, the fraction of sample A's junction-chain repertoire also detected in B (diagonal 1). Reading across a row shows how much of that sample's repertoire each other sample recovers; the asymmetry distinguishes small, deeply-shared samples from large ones.

A group of **splice-site fuzziness** views probes how precisely read junctions land on the reference: a signed **offset spectrum** (distance of each detected donor/acceptor to its nearest reference site, exact matches excluded), a cumulative **precision profile** (% of observations within ±k bp of a reference site — a lagging curve means more sites placed a few bp off the annotated position), and a **canonical split** contrasting canonical vs non-canonical sites. Together these expose imprecision that the older three-bucket CV summary collapsed away, and drive the `perc_sites_imprecise` flag. A **5'/3' read-end completeness profile** shows the cumulative fraction of reads whose ends reach the annotated gene boundary; a sample whose 5' curve lags its peers has systematically shorter 5' ends (consistent with truncation/degradation, but also with alternative TSS or incomplete annotation, so it is read as a comparative signal rather than a verdict).

Several **per-sample cohort-context pages** plot one QC rate per sample against the cohort median and, where defined, the configured warn/fail thresholds, with a marker on whichever samples the outlier scorecard flagged: **RT-switching rate** and **intra-priming rate** (from the sequence-based artefact heuristics), the **ISM fragment fraction** (share of a sample's incomplete-splice-match reads that are 3'/5'/internal fragments — a *shape* signal, read cohort-relatively, with no absolute threshold because fragment-shaped ISMs arise from alternative isoforms as well as degradation), and the **novel non-canonical junction burden** (the junction class most enriched for alignment/calling artefacts, though not exclusively artefactual). All of these degrade gracefully — a page is skipped when its underlying column was not computed.

A further set of **between-sample / fuzziness comparison views** reads *across* samples rather than judging any one in isolation. All are cohort-relative or descriptive (never an absolute pass/fail), and each is skipped when its input is absent:

- **Composition drift** — a pairwise Jensen–Shannon distance heatmap of samples' structural-category composition. A sample far from all others has an atypical composition, worth checking whether it is a genuine biological difference (a distinct condition) or a technical one. Feeds the scorecard (`composition_drift`).
- **Junction yield vs depth** — distinct junctions recovered per 1000 reads. Because it is depth-normalised, a low value means genuinely lower junction complexity (degradation or a simpler transcriptome), not merely fewer reads. Feeds the scorecard (`jxn_per_1k_reads`).
- **Tandem splice sites (NAGNAG)** — excess mass at ±3 (and ±4/±6) bp in the splice-site offset distribution, the signature of tandem donor/acceptor usage: reads selecting a nearby alternative site, a real biological phenomenon distinct from random imprecision. Purely descriptive — a higher tandem fraction is a property of the sample, not a defect.
- **Replicate concordance (multi-axis)** — for samples sharing a design-factor level (replicates), how closely each agrees with its group on structural composition, read-length profile, and splice-site precision (1 = matches its replicates). A more sensitive check than UJC overlap. Each view also reports the group's **robust within-group spread (MAD)** on every axis (hover in HTML, footnote in PDF) — the median absolute deviation of a per-sample scalar (composition Jensen–Shannon distance to the group mean, median read length, % imprecise sites), where 0 means the replicates are identical; this stays informative even for a group of exactly two replicates, where the pairwise agreement bars are identical by construction. Shown only with a `--factor` that has a ≥2-replicate level.
- **Replicate concordance of splice-site precision** — whether replicates place the *same* reference splice sites at the *same* sub-bp offsets. High agreement means the imprecision pattern is reproducible; a replicate that disagrees has sample-specific boundary placement. Also factor + ≥2-replicate only.
- **Splice-site imprecision vs junction depth** — the share of imprecise observations as a function of a junction's read depth. Imprecision concentrated at low depth is consistent with limited support; imprecision that persists at high depth points to a systematic offset. Requires per-junction coverage (`total_coverage_unique`); the view is omitted when that column is not populated.

A final **transcript-divergency** section reports TD (transcript divergency; Monzó, Frankish & Conesa 2025, *Genome Research*) — the stochastic population of RNA molecules that diverge from a sample's condition-defining transcriptional state, distinct from transcriptional noise and from alternative splicing. As an operational proxy (the paper proposes no closed formula), for each annotated gene it computes **TD = (NIC + NNC) / (FSM + NIC + NNC)** — the fraction of the gene's high-confidence reads coming from novel isoforms rather than the reference-matching (FSM) isoform — over genes with at least `transcript_divergency.min_gene_reads` (default 10) such reads. The section shows the **mean per-gene TD per sample** in cohort context, the **per-gene TD distribution** per sample, and, **when faceted, one per-gene TD-difference histogram per group pair** (reads pooled within each group per gene; genes higher in each group shown either side of zero, à la the reference study). It is descriptive only — TD feeds no QC flag or scorecard, and elevated TD may reflect genuine biology (e.g. cellular stress) or residual technical artefact, which the proxy alone cannot separate.

When the matching evidence was supplied to SQANTI3, a set of **orthogonal-support views** puts each sample's independent end/junction confirmation in cohort context (per-sample bar vs. cohort median; each omitted when its evidence is absent): **CAGE 5'-end support** — % of reads whose 5' end falls inside a CAGE peak (needs `--CAGE_peak`); **polyA 3'-end support** — % of reads whose 3' end sits at a polyA site (needs `--polyA_peak` / `--polyA_motif_list`); and **short-read support** — median short-read TSS ratio and % of junctions with short-read coverage (needs `--short_reads` / `--coverage`). These are cohort-relative comparisons, not absolute pass/fail: a sample below its peers has less independently-supported ends/junctions, worth a look. When any of these evidence types was **not** provided, the report closes with a short, non-judgemental **optional-support note** naming the missing type(s) and the SQANTI3 flag(s) that would add them — the report is complete without them, and none are required to interpret the QC results.

### Configurable thresholds

`--config my_config.yaml` overrides the built-in cut-offs; any omitted value keeps its default, so a run without `--config` is unchanged. Overridable keys include the per-sample **QC-flag thresholds** (`qc_flags`, driving the pass/warn/fail panel and `qc_summary.json`), the **intra-priming** `%A`-downstream cut-off (`intrapriming_perc_A_cutoff`, default 60), and the **PCA** cumulative-variance cut-off (`pca_cumulative_variance`, default 0.85). The splice-site fuzziness and completeness analyses add the **`jxn_offset_window`** (bp half-width for the offset spectrum / precision profile, default 15) and **`completeness_window`** (bp within which a read end counts as complete, default 50). The transcript-divergency section adds **`transcript_divergency.min_gene_reads`** (default 10) — the minimum FSM+NIC+NNC reads a gene needs to enter the per-gene TD calculation, below which TD is too noisy to be meaningful. Two additional artefact/quality flags are configurable in `qc_flags`: **`perc_novel_noncanonical_jxn`** (novel non-canonical junctions as % of all junctions; warn 2 / fail 5) has an absolute threshold because it targets a genuine artefact class, whereas the ISM fragment fraction is deliberately *not* thresholded (it is a shape signal, scored only cohort-relatively). The cohort-relative outlier detector is configured under **`sample_scorecard`** — its z-score warn/fail thresholds (`z_warn` 2.5 / `z_fail` 3.5), how many divergent metrics trip the overall flag (`min_metrics_warn` / `min_metrics_fail`), the minimum cohort size (`min_samples`, default 4), and the `metrics` map (each metric plus whether high or low values are worse). Example:

```yaml
qc_flags:
  perc_reads_intrapriming:
    warn: 25.0
    fail: 45.0
    worse: high
    label: "Intra-priming reads (%)"
intrapriming_perc_A_cutoff: 60
pca_cumulative_variance: 0.85
```

# Running example SQANTI-reads

There is an example dataset with three small samples on the SQANTI3/example/sqanti_reads_test/ directory. If you run it from the main SQANTI3 directory, you can run it in simple mode with the following command:

```
python sqanti3_reads.py --design ./example/sqanti_reads_test/sqR_design_file.csv -i ./example/sqanti_reads_test/ -p SQ_R --refGTF ./example/gencode.v38.basic_chr22.gtf --refFasta ./example/GRCh38.p13_chr22.fasta --output ./example/sqanti_reads_test/ --report both
```

## Citing SQANTI-reads
If you are using SQANTI-reads, please cite:  
Keil N, Monzó C, McIntyre L, Conesa A (2025). Quality assessment of long read data in multisample lrRNA-seq experiments with SQANTI-reads. Genome Research, 35 (4), 987. DOI: 10.1101/gr.280021.124
  
SQANTI-reads is based, and uses SQANTI3, please also cite:  
Pardo-Palacios FJ, Arzalluz-Luque A, Kondratova L, Salguero P, Mestre-Tomás J, Amorín R, Estevan-Morió E, Liu T, Nanni A, McIntyre L, Tseng E, Conesa A (2024). SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. Nature Methods, 21, 793-797. DOI: 10.1038/s41592-024-02229-2