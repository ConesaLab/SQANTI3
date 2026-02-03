
## Invocation within SQANTI3 <a name="invocation"></a>

- You pass `--bugsi human` (or `mouse`) on the SQANTI3 command line.
- In `src/qc_pipeline.py`, SQANTI3 sees `args.bugsi` and calls:
  ```python
  generate_bugsi_report(bugsi, outputClassPath, args.isoforms)
  ```

<a name="driving-script-inputs"></a>
## Driving script & inputs

- `generate_bugsi_report()` (in `src/qc_output.py`) builds and executes:
  ```bash
  Rscript /…/utilities/report_qc/BUGSI_report.R \
    <classification.txt> \
    <bugsi_<human|mouse>.gtf> \
    <your_transcript.gtf> \
    <utilities_path>
  ```
- **Inputs**:
  - `classification.txt`: SQANTI3 classification of each isoform
  - `bugsi_<species>.gtf`: gold-standard GTF of known BUGSI genes (with `ensembl`/`refseq`/`gene_name` fields)
  - `transcript.gtf`: your full transcript GTF
  - path to the utilities directory

<a name="core-r-pipeline"></a>
## Core R pipeline (`BUGSI_report.R`)

1. **Load libraries**: ggplot2, dplyr, rtracklayer, Gviz, rmarkdown, etc.
2. **Import data**  
   - `classification_data` ← read SQANTI3 TSV  
   - `bugsi_gtf` ← `rtracklayer::import()` → extract gene‐level table  
   - `transcript_gtf` ← `rtracklayer::import()`
3. **ID-type classification**  
   - Classify each `associated_gene` as "ensembl", "refseq", "gene_name", or "unknown" via regex  
   - Choose dominant ID; if Ensembl, strip version suffixes from `gene_id` in transcripts
4. **Clean & explode fusion records**  
   - Drop fusion records, then re-add them with split `associated_gene` lists  
   - Strip transcript/gene version suffixes and dedupe
5. **Define benchmarking sets**  
   - **BUGSI_transcripts**: isoforms whose `associated_gene` is in the gold list  
   - **TP (True Positives)**: subset of BUGSI_transcripts with `subcategory == "reference_match"`  
   - **PTP (Partial TP)**: FSM/ISM but not RM  
   - **FP (False Positives)**: novel categories (NIC, NNC, genic, fusion)  
   - **FN (False Negatives)**: gold‐standard genes with no FSM/ISM hit
6. **Compute metrics**  
   - Sensitivity = # unique TP genes / # gold‐standard genes  
   - Non-redundant Precision = TP / total BUGSI_transcripts  
   - Redundant Precision = (TP + PTP) / total BUGSI_transcripts  
   - Positive Detection Rate = # unique (TP+PTP) genes / # gold‐standard genes  
   - False Discovery Rate = (total BUGSI_transcripts – TP) / total BUGSI_transcripts  
   - False Detection Rate = FP / total BUGSI_transcripts  
   - Redundancy = (FSM + ISM) / # unique (TP+PTP) genes
7. **Render report**  
   - Tabulate and round percentages  
   - Assign each isoform to "TP", "PTP", "FP", or "Missing" (for FN)  
   - Call `rmarkdown::render()` on `SQANTI3_BUGSI_Report.Rmd`

<a name="html-css-js"></a>
## HTML/CSS/JS

- `bugsi_style.css` and `bugsi_script.js` accompany the Rmd to style the interactive report.

<a name="output"></a>
## Output

- `<your_prefix>_BUGSI_report.html` in your output directory, containing summary tables, bar/pie charts of TP/PTP/FP/FN, and interactive drill‑downs.

> **In short:**  
> BUGSI cross‑links your SQANTI3 classification against a curated GTF of known single‑isoform genes, segments isoforms into TP/PTP/FP/FN, computes standard metrics, and wraps everything in a self‑contained RMarkdown HTML report.

<a name="gene-selection-pipeline"></a>
## BUGSI Gene Selection Pipeline

### 1. Annotation Curation  
- Retrieved GTFs from MANE Select (Human), GENCODE (Human/Mouse), and NCBI RefSeq (Human/Mouse).  
- **Cross-validation**: kept only genes with a single, perfectly matching isoform across all sources (splice junctions, TSS, TTS).  
- **Initial candidates**: 1,925 human genes; 2,345 mouse genes.

### 2. Expression Filtering  
- Quantified median expression using GTEx (Human) and ENCODE (Mouse) RNA‑seq.  
- **Tissue-specific sets**: ≥ 5 TPM in at least one tissue.  
- **Universal set**: ≥ 1 TPM across every evaluated tissue.  
- Integrated housekeeping genes from HRT Atlas v1.0.

### 3. Alternative Splicing Exclusion  
- **Multi-exon genes**:  
  - Extracted annotated junction coverages from Recount3.  
  - Computed μ = (Σ Cᵢ) / n.  
  - Threshold T = α × μ (α = 0.01).  
  - Excluded any gene with novel junction coverage Cₙₒᵥₑₗ > T.  
  - Used IntroVerse (Human) to remove genes with novel junctions in > 50% of GTEx samples per tissue.  
- **Single-exon genes**:  
  - Overlapped coordinates with refTSS; excluded any with alternative TSS evidence.

### 4. Expert Manual Curation  
- Collaborated with GENCODE annotation experts to verify no plausible alternative isoforms.

### 5. Final Sets  
- **Human**: 53 BUGSI genes  
- **Mouse**: 37 BUGSI genes  
- Tissue‑specific BUGSI gene lists are available at the BUGSI portal (https://bugsi.uv.es). 