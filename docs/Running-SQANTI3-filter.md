### Table of contents:

- <a href="#rules">Rules filter: removing artifacts using the SQANTI3 output and user-defined rules</a>
    * <a href="#running_rules">Running the rules filter</a>
    * <a href="#default_filter">Default filter</a>
    * <a href="#build_json">Making your own filter</a>

- <a href="#rulesout">Rules filter output</a>

- <a href="#ml">Machine learning filter: removing artifacts using a random forest classifier</a>
    * <a href="#runml">Running the machine learning filter</a>
    * <a href="#sets">True Positive (TP) and True Negative (TN) sets</a>
    * <a href="#partition">Partitioning of the training data (model training/testing)</a>
    * <a href="#prob">Probability threshold to define Isoforms/Artifacts</a>
    * <a href="#cols">Classification file columns excluded from the ML filter</a>
    * <a href="#colrem">Excluding additional data columns from the filter</a>
    * <a href="#intraprim">Intra-priming filter</a>
    * <a href="#retain">Forcing the removal/retention of specific isoform groups</a>
    * <a href="#except">Exceptions to running the SQANTI3 ML filter</a>

- <a href="#mlout">Machine learning filter output</a>
    * <a href="#input">Input-related files</a>
    * <a href="#model">Model-related files</a>
    * <a href="#test">Test set files</a>
    * <a href="#reportml">Filter report</a>
    * <a href="#classifml">Classification file</a>
    * <a href="#inclusionml">Inclusion list and filtering of SQANTI QC output files</a>


__________________________________________

<a id="rules"></a>

## Rules filter: removing artifacts using the SQANTI3 output and user-defined rules
</a>
SQANTI3 filter has been reshaped to allow users to supply their own set of rules. These will be used to accept or discard an isoform based on the attributes obtained through [SQANTI3 QC](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control). To define those rules, the user will need to create a **JSON file** specifying which SQ3 attributes and thresholds are to be used for filtering and, therefore, the properties of the transcripts that are to be included in the final transcriptome. 

<img src = "https://github.com/aarzalluz/figures_public/blob/master/SQANTI3/SQ3_rules-filter.png" width = "626" height = "977">

<a id="running_rules"></a>

### Running the rules filter
</a>
The rules filter only needs a SQANTI3 classification file to work. To run it, just execute `sqanti3_filter.py` providing the `rules` argument followed by the path to the SQANTI3 `*_classification.txt` file:

```
python sqanti3_filter.py rules --sqanti_class path/to/classification.txt 
``` 

Using the optional arguments, the script can be set up to filter several file formats (such as **FASTA**, **GTF** and **BAM**) and remove the isoforms that were catalogued as artifacts based on these user-defined rules. This option is designed to make it easy for users filter the files from the [SQANTI3 QC output](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC), enabling quick transcriptome curation.

These are the arguments accepted by `sqanti3_filter.py rules`:

```
usage: sqanti3_filter.py rules [-h] --sqanti_class SQANTI_CLASS [--isoAnnotGFF3 ISOANNOTGFF3] [--filter_isoforms FILTER_ISOFORMS] [--filter_gtf FILTER_GTF] [--filter_sam FILTER_SAM]
                               [--filter_faa FILTER_FAA] [-o OUTPUT] [-d DIR] [--skip_report] [-e] [-v] [-c CPUS] [-l {ERROR,WARNING,INFO,DEBUG}] [-j JSON_FILTER]
```

<details><summary> Arguments description</summary>

```
Rules filter selected

options:
  -h, --help            show this help message and exit

Required arguments:
  --sqanti_class SQANTI_CLASS
                        SQANTI3 QC classification file.

Input options:
  --isoAnnotGFF3 ISOANNOTGFF3
                        isoAnnotLite GFF3 file to be filtered
  --filter_isoforms FILTER_ISOFORMS
                        fasta/fastq isoform file to be filtered
  --filter_gtf FILTER_GTF
                        GTF file to be filtered
  --filter_sam FILTER_SAM
                        SAM alignment of the input fasta/fastq
  --filter_faa FILTER_FAA
                        ORF prediction faa file to be filtered by SQANTI3

Output options:
  -o OUTPUT, --output OUTPUT
                        Prefix for output files.
  -d DIR, --dir DIR     Directory for output files. Default: ./sqanti3_results
  --skip_report         Skip creation of a report about the filtering

Filtering options:
  -e, --filter_mono_exonic
                        All mono-exonic transcripts are automatically filtered

Extra options:
  -v, --version         Display program version number.
  -c CPUS, --cpus CPUS  Number of CPUs to use. Default: 4
  -l {ERROR,WARNING,INFO,DEBUG}, --log_level {ERROR,WARNING,INFO,DEBUG}
                        Set the logging level.
                        Default: INFO

Rules specific options:
  -j JSON_FILTER, --json_filter JSON_FILTER
                        JSON file where filtering rules are expressed. Rules must be set taking into account that attributes described in the filter will be present in those isoforms that should be kept.
                        Default: /home/pabloati/Programs/sqanti3/src/utilities/filter/filter_default.json
```
</details><br>

<a id="default_filter"></a>

### Default filter
</a>

Following the same criteria as in previous versions of SQANTI3, we hereby provide a default set of rules for the filter. These are equivalent to running the old `sqanti3_RulesFilter.py` script. However, we strongly advise users to closely inspect the properties of their transcriptome and define their own filter based on the results obtained during QC. Keep in mind as well that any isoform that has a missing value in any of the conditions of a rule will be considered an artifact. This means that `NA` the likes will be treated as fails. 

When no JSON file is provided to the rules filter, this is the set of filtering rules used by default:

* If a transcript is a **FSM**, it is kept unless the 3' end is unreliable because of a possible intrapriming event. 
    * We consider that an isoform is likely to be an artifact of intrapriming if in the 20bp downstream the annotated Transcription Terminating Site (TTS) there are 12 or more adenines at the genomic level. That means that following the TTS there is a stretch of nucleotides with at least a 60% A's. 
* If a transcript is **not a FSM**, it is kept only if all of the criteria below are met: 
    * (1) 3' end is not an intrapriming artifact.
    * (2) does not have a junction that is labeled as RT-Switching.
    * (3) all junctions are either canonical or have short read coverage above the `-c` threshold.

This filter is supplied in the `utilities/filter/filter_default.json` file:

```json
{
    "full-splice_match": [
        {
            "perc_A_downstream_TTS":[0,59]
        }
    ],
    "rest": [
        {
            "perc_A_downstream_TTS":[0,59],
            "RTS_stage":"FALSE",
            "all_canonical":"canonical"
        },
        {   
            "perc_A_downstream_TTS":[0,59],
            "RTS_stage":"FALSE",
            "min_cov":3
        }
     ]
}
```


<a id="build_json"></a> 

### Making your own filter
</a>
As of SQANTI3 version 5.0, it is possible to supply a user-defined filter including as many rules as desired. To do it, the user needs to create a JSON file defining which characteristics make an isoform reliable. 

#### Defining rules and requisites

Specifically, for each structural category, the user should build an array of objects (or rules). A **rule** is made of one or more **requisites**, all of which must be fulfilled for an entry to be considered a true **Isoform** (they will be evaluated as AND in terms of logical operators). If different rules (i.e. sets of requisites) are defined for the same structural category, they will be treated independently of one another. In that case, to pass the filter, transcripts will need to pass at least one of these independent rules (they will be evaluated as OR in terms of logical operators). 

Of note, it is also possible to set rules for the rest of structural categories, i.e. those rules will be applied to any structural category for which there isn't a rule established. In that case, users will need to define the rule using the category name `rest`.

Here is a basic scheme of how to define rules and the list of requisites that make up a rule in a JSON file:

```json
{
    "structural_category": [
        {
            "numeric_column_name":[0,100]
            "character_column_name": ["accepted_value1", "accepted_value2"]
        }
     ]
}
```

Rules can be set for any numeric or character column in the classification file:
- **Numeric values**: in this case, it is possible to define an interval (`[X,Y]`) or set just the lower limit (`X`). Please, take into account that the **limit values will be included**. 
- **Character and logical columns**: such as `subcategory`, `RTS_stage` or `all_canonical`. In this case, users can simply establish which terms will be accepted. Of note, users may want to accept several of the values in the column (for instance, several subcategories). If so, the requisite can be defined as an array, and the filter will keep the entries if any of those values are present in the specified column.

#### User-defined rules (JSON file) example

As an example, let's say that we want to define a custom filter that will keep only isoforms that pass these rules:

* **FSM**: Keep if the isoform is:
    * NOT a potential intrapriming product (less than 60% in the `perc_A_downstream_TTS` column).
* **ISM**: Keep if:
    * It is larger than 2kb and shorter than 15kb **AND**
    * It is catalogued within the "3prime_fragment", "5prime_fragment" or "internal_fragment" subcategories.
* **NIC**: Keep if:
    * All the splice junctions are canonical OR covered at least by 10 short reads.
* **NNC**: Keep if:
    * All the splice junctions are canonical OR covered at least by 10 short reads **AND**
    * The absolute distance to the closest annotated TSS and TTS is 50bp or less.
* **Rest of the categories**: Keep if :
    * All the splice junctions are canonical and they are not suspicious of being an RT-Switching (RTS) artifact.
    * It is a coding transcript.
    * It is not an intrapriming product.
    * It is not a mono-exon transcript.

The JSON file for that precise filter will look like this:

```json
{
    "full-splice_match": [
        {
            "perc_A_downstream_TTS":[0,59]
        }
    ],
    "incomplete-splice_match":[
        {
            "length":[2001,14999],
            "subcategory": ["3prime_fragment", "5prime_fragment", "internal_fragment"]
        }
    ],
    "novel_in_catalog":[
        {
            "all_canonical": "canonical"
        },
        {
            "min_cov": 10
        }
    ],
    "novel_not_in_catalog":[
        {
            "all_canonical": "canonical",
            "diff_to_gene_TSS":[-50,50],
            "diff_to_gene_TTS": [-50,50]            
        },
        {
            "min_cov": 10,
            "diff_to_gene_TSS":[-50,50],
            "diff_to_gene_TTS": [-50,50] 
        }
    ],
    "rest": [
        {
            "RTS_stage":"FALSE",
            "all_canonical":"canonical",
            "coding": "coding",
            "perc_A_downstream_TTS":[0,59],
            "exons": 2
        }
     ]
}

```


<a id="rulesout"></a>

## Rules filter output
</a>
The main SQANTI rules filter output files are:

* `*_RulesFilter_result_classification.txt`: New classification file with an extra column called "filter_result" that will classify all the entries as *Isoform* or *Artifact*.
* `*_inclusion-list.txt`: Text file with the IDs of the isoforms that passed the filter.
* `*_filtering_reasons.txt`: TSV file with 3 columns:
    * (1) Isoform ID of the isoform that was. filtered out
    * (2) Structural Category
    * (3) Reason why the isoform was discarded as an artifact. If an isoform is catalogued as *Artifact* because it doesn't fulfill several rules, there will be multiple lines in this file regarding that isoform.
* `*_SQANTI3_filter_report.pdf`: A PDF report with some plots describing the performance of the filtering.

<a id="ml"></a>

## Machine learning filter: removing artifacts using a random forest classifier
</a>
SQANTI3 incorporates an automated filter based on a **random forest classifier**, 
which is designed to discriminate potential artifacts from true isoforms
without the need for user-defined rules or manually-set thresholds. 

Briefly, this classifier learns high and low quality attributes from a True 
Positive (TP) and True Negative (TN) transcript set, building a model that 
discriminates artifacts and isoforms based on TN and TP properties. 
The filter returns a modified  version of the `*_classification.txt` file (named `*_MLresult_classification.txt`). The new table will include a `filter_result` 
column with the label `Isoform` or `Artifact`, which will be the result of the 
random forest classifier.

⚠️ **Disclaimer:** we hereby provide instructions to run the filter and a comprehensive 
description of its parameters. While we have included recommendations on how to best tune 
them to obtain optimal results, users are encouraged to try several configurations of 
the filter to find the best way to run it for their particular dataset. ⚠️

<img src = "https://github.com/aarzalluz/figures_public/blob/master/SQANTI3/SQ3_MLfilter.png" width = "700" height = "904">


<a id="runml"></a>

### Running the machine learning filter
</a>
Similarly to the rules filter, the machine learning filter (ML filter) can be 
run using the `sqanti3_filter.py` script by providing the `ML` argument followed by
the path to the SQANTI `*_classification.txt` file (output by `sqanti3_qc.py`):

```
python sqanti3_filter.py ml path/to/classification.txt 

```

A brief **tutorial** on how to run the SQANTI3 ML filter for an example dataset can 
be found [here](https://github.com/ConesaLab/SQANTI3/wiki/Tutorial:-running-SQANTI3-on-an-example-dataset).

These are the parameters and arguments accepted by `sqanti3_filter.py ML`:

```bash
usage: sqanti3_filter.py ml [-h] --sqanti_class SQANTI_CLASS [--isoAnnotGFF3 ISOANNOTGFF3]
                            [--filter_isoforms FILTER_ISOFORMS] [--filter_gtf FILTER_GTF]
                            [--filter_sam FILTER_SAM] [--filter_faa FILTER_FAA] [-o OUTPUT] [-d DIR]
                            [--skip_report] [-e] [-v] [-c CPUS] [-l {ERROR,WARNING,INFO,DEBUG}]
                            [-t PERCENT_TRAINING] [-p TP] [-n TN] [-j THRESHOLD] [-f] [--intermediate_files]
                            [-r REMOVE_COLUMNS] [-z MAX_CLASS_SIZE] [-i INTRAPRIMING]
```

<details><summary> Arguments description</summary>

```bash
ML filter selected

options:
  -h, --help            show this help message and exit

Required arguments:
  --sqanti_class SQANTI_CLASS
                        SQANTI3 QC classification file.

Input options:
  --isoAnnotGFF3 ISOANNOTGFF3
                        isoAnnotLite GFF3 file to be filtered
  --filter_isoforms FILTER_ISOFORMS
                        fasta/fastq isoform file to be filtered
  --filter_gtf FILTER_GTF
                        GTF file to be filtered
  --filter_sam FILTER_SAM
                        SAM alignment of the input fasta/fastq
  --filter_faa FILTER_FAA
                        ORF prediction faa file to be filtered by SQANTI3

Output options:
  -o OUTPUT, --output OUTPUT
                        Prefix for output files.
  -d DIR, --dir DIR     Directory for output files. Default: ./sqanti3_results
  --skip_report         Skip creation of a report about the filtering

Filtering options:
  -e, --filter_mono_exonic
                        All mono-exonic transcripts are automatically filtered

Extra options:
  -v, --version         Display program version number.
  -c CPUS, --cpus CPUS  Number of CPUs to use. Default: 4
  -l {ERROR,WARNING,INFO,DEBUG}, --log_level {ERROR,WARNING,INFO,DEBUG}
                        Set the logging level.
                        Default: INFO

Machine Learning specific options:
  -t PERCENT_TRAINING, --percent_training PERCENT_TRAINING
                        Proportion of the data that goes to training (parameter p of the function createDataPartition).                     
                        Default: 0.8
  -p TP, --TP TP        Path to file containing the list of the TP transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.
  -n TN, --TN TN        Path to file containing the list of the TN transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.
  -j THRESHOLD, --threshold THRESHOLD
                        Machine Learning probability threshold to classify transcripts as positive isoforms.         
                        Default: 0.7
  -f, --force_fsm_in    Forces retaining FMS transcripts regardless of ML filter result (FSM are threfore automatically classified as isoforms).
  --intermediate_files  Outputs ML filter intermediate files.
  -r REMOVE_COLUMNS, --remove_columns REMOVE_COLUMNS
                        Path to single-column file (no header) containing the names of the columns in SQ3's classification.txt file that are to be excluded during random forest training (optional).
  -z MAX_CLASS_SIZE, --max_class_size MAX_CLASS_SIZE
                        Maximum number of isoforms to include in True Positive and True Negative sets. TP and TN sets will be downsized to this value if they are larger.
                        Default: 3000
  -i INTRAPRIMING, --intrapriming INTRAPRIMING
                        Adenine percentage at genomic 3' end to flag an isoform as intra-priming. Default: 60
```
</details><br>

We next provide a detailed description of each of the ML filter parameters in order
to help users understand how to best apply it to their own data.


<a id="sets"></a>

### True Positive (TP) and True Negative (TN) sets
</a>
As stated above, the random forest classifier model needs to be trained on a subset
of isoforms from the transcriptome. This will include a set of **True Positive isoforms**
-that is, isoforms that are reliable enough for the classifier to learn the properties 
of a *real isoform*- and a set of **True Negative isoforms** -that is, isoforms that
can be considered low-quality or false, from which the classifier will learn what
constitutes an *artifact* or false-positive isoform-.

Although the ML filter includes built-in selection of training data from specific 
SQANTI categories, it also accepts **user-defined TN and TP sets** (`--TN [-n]` 
and `--TP [-p]` parameters, respectively). These should be provided as 
single-column files (with no header) containing the IDs of the selected isoforms, 
and will be used to train and test the random forest classifier.

If not provided, the SQANTI ML filter will use the following as 
**built-in TP and TN sets** by default:

- **TP:** all **Reference Match (RM)** isoforms in the transcriptome (multi-exon only). 
RM are an FSM subcategory in which the TSS and TTS of the isoform are within +/-50bp
of the TSS/TTS of the reference associated transcript (see the category description page
for details). If there are less than 250 RM isoforms, all FSMs will be taken as TP.
- **TN:** all **Novel Not in Catalog (NNC)** isoforms that have at least one non-canonical
junction (multi-exon only). If there are less than 250 isoforms in this group (normally
because there are not enough non-canonical junctions in the transcriptome), all NNC 
isoforms will be taken as TN.

Note that both user-defined and built-in training data will be subset if the number of
isoforms in the TP or TN isoform sets is larger than the other. By default, 
SQANTI3 will **downsize the larger set** to match the smaller of the two sizes.
Downsizing is performed by random sampling. When using the built-in option, and
given this random sampling step, **TP and TN lists will be written out as part of the SQANTI ML filter output**.

Users may change the `--max_class_size [-z]` parameter to set a maximum
number of isoforms to include in the TP and TN sets. By default, TP and TN lists
will be **downsized to 3000 transcripts** if they are larger than this, 
no matter if they were user-defined or generated internally.


<a id="partition"></a>

### Partitioning of the training data (model training/testing)
</a>
TP and TN isoform sets need to be split into *training* and *test* data to correctly
generate the random forest model. Using the `--percent_training` or `-t` parameter,
users can specify the proportion of the data that is used for training. By default,
**80% of the isoforms will be used for training** (which is equivalent to supplying 
`-t 0.8`). The remaining 20% will be used to test the model.

The ML filter uses functions from the `caret` R package to build and test the random
forest classifier. Briefly, the model training/test includes the following internal 
steps:

1. Data partitioning via the `createDataPartition()` function.
2. Cross-validation (10x) using the `trainControl()` function.
3. Model training using the `train()` function.
4. Testing of the model using the `predict()` function from the `stats` package 
in base R.
5. Output model statistics using several other functions from the `caret` package 
(see the ML filter output section for details).

**Warning:** after training, the model will be stored as an `.RData` object in 
the specified output folder. If the filter is re-run on the same directory, SQANTI
will detect that a model already exists, skipping the training step and
**applying the extant model** to perform the Isoform/Artifact classification on the input data. 
This is useful on occasions where a user wants to apply an already-generated 
model to a different transcriptome; however, if you wish to train a new model, you will need to *delete/move/rename the previous model object*.



<a id="prob"></a>

### Probability threshold to define Isoforms/Artifacts
</a>
Ultimately, flagging transcripts as Isoforms or Artifacts is based on the 
random forest **probability** to correctly classify a transcript into either
group. These probabilities are included in the `POS_MLprob` and `NEG_MLprob` 
columns of the output `_MLresult_classification.txt` file. Intuitively:

- `POS_MLprob` is the probability of classifying the transcript as an Isoform.
- `NEG_MLprob` is the probability of classifying the transcript as an Artifact, and
is equivalent to `1 - POS_MLprob`.

The SQANTI ML filter will assign the Artifact and Isoform labels based on a 
**probability threshold parameter** supplied via the `--threshold` or `-j`
argument. By default, the filter is more stringent on the Isoform condition than
it is on the Artifact condition, that is: by default, **transcripts will be considered 
isoforms when `POS_MLprob` is larger than 0.7**. 

This results in some artifacts having `POS_MLprob > 0.5` (which can be 
counter-intuitive). Of course, in this scenario, the filter will exclude more 
transcripts than when using the usual `POS_MLprob > NEG_MLprob` approach. 
However, users may reproduce this lenient filter by lowering the threshold to 0.5.



<a id="cols"></a>

### Classification file columns excluded from the ML filter
</a>
Due to their lack of importance for artifact definition (or to prevent them from
having unwanted effects on the filtering), the following columns are removed 
from the classification table prior to running the SQANTI ML filter:

- Chromosome and strand info (`chrom`, `strand`)
- Reference gene/transcript info (`associated_gene`, `associated_transcript`)
- Structural information about the reference or the transcript (`ref_length`,
`ref_exons`, `ORF_length`, `CDS_length`, `CDS_start`, `CDS_end`, `CDS_genomic_start`,
`CDS_genomic_end`).
- Sequence information (`seq_A_downstream_TTS`, `ORF_seq`).
- SQANTI category (`structural_category`, `subcategory`).

Moreover, the following columns are excluded when TP and TN sets are built by
the ML filter, instead of supplied by the user:
- Redundant junction information (`all_canonical`) to prevent overfitting due to NNC non-canonical
or NNC transcripts being used as a TN set.
- Redundant distance to TSS/TTS information (all `dist_*` columns in the classification
file) to prevent bias towards reference-similar isoforms when RM are used as a TP set.


<a id="colrem"></a>

### Excluding additional data columns from the filter
</a>
The SQANTI ML filter incorporates the option to exclude columns from the 
*_classification.txt* file from model training to prevent filtering based on 
them. These should be provided as a single-column file including the column names
exactly as they are named in the classification file. A path to this file should
then be supplied via the `--remove_columns` or `-r` argument.

This is particularly useful to prevent **overfitting** during model training. 
In particular, there may be situations in which users wish to define TP and TN 
sets based on the values of specific SQ3 columns. We recommend that these are, 
in turn, excluded from the ML filter to prevent them from gaining too much 
importance during Isoform/Artifact classification. 

*Warning: we always encourage testing several ML filter configurations and TP/TN sets 
before making a final decision on which are the artifacts in your transcriptome.*


<a id="intraprim"></a>

### Intra-priming filter
</a>
The ML filter script in SQANTI3 includes a simple **threshold-based filter** to prevent
potential **intra-priming artifacts** from remaining in the transcriptome, which is
not specifically considered by the ML filter. Still, note that the `perc_A_downstream_TTS`
column, which is used to flag transcripts as potentially resulting from intra-priming,
is considered during random forest training.

To apply it, users may supply the **maximum percentage of adenines** ("A" nucleotides)
that is to be allowed in the genome sequence region immediately following the end of the
transcript (i.e. the TTS). This can be provided via the `--intrapriming` or `-i` 
argument, resulting in `perc_A_downstream_TTS > i` isoforms being flagged as 
intra-priming artifacts. By default, transcripts containing more than 60% A's at the genomic
3' end will be flagged.


<a id="retain"></a>

### Forcing the removal/retention of specific isoform groups
</a>
1. **Full-splice match (FSM) transcripts** can be excluded from the ML filter and
therefore retained by providing the `--force_fsm_in` or `-f` argument. By default,
FSM will be run through the filter along with the isoforms from the rest of the 
structural categories.

2. All **mono-exonic transcripts** can be automatically removed by providing the 
`--filter_mono_exonic` or `-e` parameter. Note that mono-exon transcripts will not
be evaluated by the ML filter (only by the intra-priming filter).


<a id="except"></a>

### Exceptions to running the SQANTI3 ML filter
</a>
The random forest classifier **cannot/will not** be run in any of the following
scenarios:

- All transcripts in the long read-defined transcriptome are classified as
mono-exon. The ML filter can only be applied to multi-exon transcripts.
- One (or both) of the user-defined TP and TN sets have less than 40 isoforms.
- One of the default TP and TN categories has less than 40 isoforms. 
This applies when no user-defined TP and TN sets are provided.

However, even in scenarios where no ML-based filtering can be applied, the
**intra-primming filter** will still be applied to all input transcripts.


<a id="mlout"></a>

## SQANTI ML filter output
</a>

The SQANTI ML filter output files are written to the path specified using the
`--dir` or `-d` argument, also appending the prefix provided via the `--output` 
or `-o` argument.

An example of the ML filter output can be found under the `example/MLfilter_output` folder, which is included in
the main SQANTI3 directory (see our [example dataset tutorial](https://github.com/ConesaLab/SQANTI3/wiki/Tutorial:-running-SQANTI3-on-an-example-dataset) for details).

The following output files are generated after running the filter:


<a id="input"></a>

### Input-related files
</a>
- `params.txt`: two-column file including the name of the argument and the
value that was used to run the filter.
- `TP_list.txt` and `TN_list.txt`: a single-column text file including the 
IDs of the isoforms used as TP and TN (generated only when no user-defined sets
were provided).

The supplied prefix will be appended to the filenames above.


<a id="model"></a>

### Model-related files
</a>
- `randomforest.RData`: R object containing the trained random forest model.
- `classifier_variagble-importance_table.txt`: two-column text file including the
names of the variables that were used for classification and a numeric value indicating
their importance in the random forest classifier.
- `intermediate_*_MLinput_table.txt`: training-ready classification table. SQANTI
ML filter performs a series of modifications of the classification table to 
enable running of the training/test steps. These include handling of `NAs`, value
formatting, column removal, etc. This file will only be output if the `--intermediate_files`
argument is supplied. Note that this file is **NOT intented for downstream anaylsis**.


<a id="test"></a>

### Test set files
</a>
All include the `testSet_` prefix:

- `stats.txt`: all statistics for the model testing.
- `ROC_curve.pdf`.
- `confusionMatrix.txt`: confusion matrix obtained after running the trained model
on the test set.
- `summary.txt`: summary statistics of the model testing, namely sensitivity,
specificity and AUROC.


<a id="reportml"></a>

### Filter report
</a>
The ML filter automatically returns a PDF report. This includes filter summary tables and plots as well as 
one diagnostic plot per variable used by the random forest classifier. This allows an evaluation of the behaviour of 
transcripts flagged as isoforms and artifacts in the light of each of the variables that were considered 
relevant to discriminate both.

Supplying the `--skip_report` flag will deactivate the report.


<a id="classifml"></a>

### Classification file
</a>
The `MLresult_classification.txt` constitutes a modified classification file 
including the results of the ML and intra-priming filters. The following
columns are added:

- `POS_MLprob` and `NEG_MLprob`, described above (see probability threshold section). 
- `ML_classifier` column, in which `ML_classifier == Positive` corresponds to transcripts
in which `POS_MLprob > t` (therefore flagged as Isoforms by the ML filter).
- `intra_priming` column, in which `intrapriming == TRUE` corresponds to transcripts
flagged as intra-priming atrifacts.
- `filter_result` column, including the combined result for the ML and intra-priming
filters. An isoform will be `filter_result == Isoform` if it passed both the
ML and intra-priming filters. 

In the specific case of mono-exonic transcripts, passing the intra-priming filter 
will be enough to be considered an isoform, as long as they are not removed using 
the `-e` argument. Similarly, all FSM transcripts will be flagged as true isoforms if
the `-f` option is provided.


<a id="inclusionml"></a>

### Inclusion list and filtering of SQANTI QC output files
</a>
Finally, `inclusion-list.txt` is a single-column file including the IDs of the transcripts 
classified as true isoforms (i.e. flagged as `Isoform` in the `filter_result` column).
This list will be used to provide filtered versions of supplied SQANTI QC output files. 
Use the `--isoAnnotGFF3`, `--isoforms`, `--gtf`, `--faa` and `--sam` arguments to provide 
each of the files you wish to obtain filtered.


