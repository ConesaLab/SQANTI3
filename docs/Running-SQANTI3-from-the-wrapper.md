## SQANTI3 v5.4: Simplified Execution with Wrapper Script and Configuration File

As of release v5.4, SQANTI3 introduces a streamlined way to run the pipeline using a unified wrapper script and a YAML configuration file. This update aims to simplify SQANTI3 execution by consolidating all parameters for a given run—from quality control (QC) through filtering and rescue—into a single, easily managed file.

---

## **Recommended Setup**

To ensure a smooth experience with SQANTI3, follow these recommendations:

- **Add SQANTI3 to your PATH:**  
  Make the wrapper script accessible system-wide by either adding the SQANTI3 directory to your PATH or creating a symbolic link to the `sqanti3` script in a directory already in your PATH.

- **Use a configuration file per experiment:**  
  Maintain one YAML configuration file per experiment. This allows you to tweak parameters directly from the command line for each run.

<details>
    <summary>Adding SQANTI3 to the PATH in Ubuntu</summary>

You can add the wrapper script to your PATH in two ways:

- **Option 1: Add the SQANTI3 directory to your PATH**
    ```bash
    export PATH=$PATH:/path/to/SQANTI3
    ```
    Add this line to your `.bashrc` file to make it persistent across sessions.

- **Option 2: Create a symbolic link**
    ```bash
    ln -s /path/to/SQANTI3/sqanti3 /usr/local/bin/sqanti3
    ```
</details><br>

---

## **Configuration File**

The configuration file is written in YAML (⚠️ must use the `.yaml` extension). It contains all arguments and parameters required by the SQANTI3 modules and flags to select which modules (qc, filter, rescue) to run. For the filter and rescue modules, you can specify whether to use rules-based or machine learning approaches.

**Structure:**
- *Common options*: Shared across all modules (e.g., reference genome, maximum threads).
- *Module-specific options*: Parameters unique to QC, filter, or rescue.

**Creating a Configuration Template:**
Generate a template with default parameters using:
```bash
sqanti3 init -c 
```
If no filename is provided, it defaults to `sqanti3_config.yaml`.

**Overriding Parameters:**
You can override configuration values directly from the command line using the `-a` flag:
```bash
sqanti3 init -c sqanti3_config.yaml -a cpus=8 dir=sqanti3_results
```
The `-a` flag can also be used when running SQANTI3 to temporarily override config values for a specific run.

**View Available Arguments:**
Use the `-h` flag after a module name to list all available arguments and their default values.

<details>
    <summary>❓Example Configuration File</summary>

```yaml
main:
  refGTF: ''
  refFasta: ''
  cpus: 4
  dir: sqanti3_results
  output: isoforms
  log_level: INFO
qc:
  enabled: true
  options:
    isoforms: ''
    min_ref_len: 0
    force_id_ignore: false
    fasta: false
    genename: false
    short_reads: ''
    SR_bam: ''
    novel_gene_prefix: ''
    aligner_choice: minimap2
    gmap_index: ''
    sites: ATAC,GCAG,GTAG
    skipORF: false
    orf_input: ''
    CAGE_peak: ''
    polyA_motif_list: ''
    polyA_peak: ''
    phyloP_bed: ''
    saturation: false
    report: html
    isoform_hits: false
    ratio_TSS_metric: max
    chunks: 1
    is_fusion: false
    expression: ''
    coverage: ''
    window: 20
    fl_count: ''
    isoAnnotLite: false
    gff3: ''
filter:
  enabled: true
  options:
    common:
      sqanti_class: sqanti3_results/isoforms_classification.txt
      isoAnnotGFF3: ''
      filter_isoforms: sqanti3_results/isoforms_corrected.fasta
      filter_gtf: sqanti3_results/isoforms_corrected.gtf
      filter_sam: ''
      filter_faa: sqanti3_results/isoforms_corrected.faa
      skip_report: false
      filter_mono_exonic: false
    rules:
      enabled: true
      options:
        json_filter: /home/pabloati/Programs/sqanti3/src/utilities/filter/filter_default.json
    ml:
      enabled: false
      options:
        percent_training: 0.8
        TP: ''
        TN: ''
        threshold: 0.7
        force_fsm_in: false
        intermediate_files: false
        remove_columns: ''
        max_class_size: 3000
        intrapriming: 60
rescue:
  enabled: true
  options:
    common:
      filter_class: sqanti3_results/isoforms_RulesFilter_result_classification.txt
      rescue_isoforms: sqanti3_results/isoforms_corrected.fasta
      rescue_gtf: sqanti3_results/isoforms.filtered.gtf
      refClassif: ''
      rescue_mono_exonic: all
      mode: automatic
    rules:
      enabled: true
      options:
        json_filter: /home/pabloati/Programs/sqanti3/src/utilities/filter/filter_default.json
    ml:
      enabled: false
      options:
        random_forest: ''
        threshold: 0.7
```
</details><br>

---

## **Running SQANTI3**

After creating your configuration file, you can run SQANTI3 in two main ways:

- **Run all activated modules (QC, filter, rescue):**
    ```bash
    sqanti3 all -c sqanti3_config.yaml
    ```

- **Run a specific module:**
    ```bash
    sqanti3 qc -c sqanti3_config.yaml
    sqanti3 filter -c sqanti3_config.yaml
    sqanti3 rescue -c sqanti3_config.yaml
    ```

**Override parameters for a single run:**
```bash
sqanti3 qc -c sqanti3_config.yaml -a short_reads=/path/to/short_reads.fastq dir=/path/to/output
```
This command temporarily changes the short reads file and output directory for this run only.

<details>
    <summary>❓Help message</summary>

```bash
usage: sqanti3 [-c CONFIG] [--dry-run] [-a ARGUMENTS [ARGUMENTS ...]] [-h] [-l {ERROR,WARNING,INFO,DEBUG}] {all,qc,filter,rescue,init}

Python wrapper for SQANTI3 pipeline.

positional arguments:
  {all,qc,filter,rescue,init}  Action to perform

options:
  -c CONFIG, --config CONFIG   Path to the configuration file (default: sqanti3_config.yaml)
  --dry-run                   Print the commands that would be executed
  -a ARGUMENTS [ARGUMENTS ...], --arguments ARGUMENTS [ARGUMENTS ...]
                              Non-default arguments to pass to the SQANTI3 modules
  -h, --help                  Show the help message and exit
  -l {ERROR,WARNING,INFO,DEBUG}, --log_level {ERROR,WARNING,INFO,DEBUG}
                              Set the logging level (default: INFO)
```
</details><br>

