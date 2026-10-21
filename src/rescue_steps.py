import os
import json
import shutil
import pandas as pd

from src.wrapper_utils import sqanti_path
from src.module_logging import rescue_logger, message
from src.commands import (
    RSCRIPTPATH, run_command, PYTHONPATH, RESCUE_RANDOM_FOREST
)
from src.utilities.rescue.automatic_rescue import (
    get_lost_reference_id, generate_automatic_table
)
from src.utilities.rescue.rescue_helpers import (
    identify_rescue_candidates,
    get_rescue_gene_targets, get_rescue_reference_targets, read_classification
)

from src.utilities.rescue.candidate_mapping_helpers import (
    filter_transcriptome, process_sam_file, save_fasta, merge_target_records,
    mapping_fingerprint, is_mapping_reusable, save_mapping_fingerprint
)

from src.utilities.rescue.rescue_by_mapping import rescue_by_mapping

from src.rescue_output import (
    write_rescue_gtf, write_rescue_fasta
)

from src.parsers import read_counts_design
from src.utilities.rescue.evidence_check import (
    rule_column_matches, target_count_rules, count_samples,
    write_evidence_rules, apply_evidence_check
)
from src.utilities.rescue.sq_requant import load_counts

MINIMAP2_RESCUE_OPTS = "--secondary=yes -ax map-hifi"
# Rules removed from the JSON file before filtering the reference transcriptome,
# because they cannot be evaluated on it. A rule is removed when its column is one
# of these names or starts with one of them followed by "." or "_".
# - FL, FL.<sample>, prevalence, prevalence_<group>: long-read counts, which the
#   reference does not have.
# - subcategory: every reference FSM is reference_match, so the rule would test the
#   annotation itself instead of the evidence supporting each reference transcript.
RESCUE_IGNORED_RULES = ["FL", "prevalence", "subcategory"]

def run_automatic_rescue(classif_df,monoexons):
    message("Performing automatic rescue",rescue_logger)
    # Select the FSM  isoforms with more than one exon 
    rescue_classif = classif_df[
        (classif_df['structural_category'].isin(['full-splice_match'])) & 
        (classif_df['exons'] > 1)
    ]

    # Include monoexonic FSM if requested
    if monoexons != 'none':
        rescue_classif = pd.concat([
            rescue_classif,
            classif_df[
                (classif_df['structural_category'].isin(['full-splice_match'])) & 
                (classif_df['exons'] == 1)
            ]
        ])
 
    # Find the references that are lost and get the ones that are not represented by isoforms
    lost_ref = get_lost_reference_id(rescue_classif)
    if len(lost_ref) == 0:
       rescue_logger.info("No lost references found")
       rescue_logger.info("Automatic rescue is not needed")
       rescue_df = generate_automatic_table(pd.DataFrame({'associated_transcript': ['none']}),classif_df)
       return pd.DataFrame(columns=['isoform']), rescue_df
    rescue_logger.debug(f"Found {len(lost_ref)} lost references")
    # Lost references with an FSM, found in one pass over the table instead
    # of filtering it once per lost reference
    fsm_refs = rescue_classif.loc[
        rescue_classif['structural_category'] == 'full-splice_match',
        'associated_transcript'
    ].unique()
    has_fsm = pd.Series(lost_ref).isin(fsm_refs).to_numpy()
    rescue = pd.DataFrame({'isoform': lost_ref[has_fsm]})

    # Split into reference transcripts and ISM TODO: Eliminate this step?
    rescue_ref = rescue[rescue['isoform'].isin(rescue_classif['associated_transcript'])]
    rescue_logger.debug(f"Rescued {rescue_ref.shape[0]} transcripts")
    
    # Save the automatic rescue
    rescue_df = generate_automatic_table(rescue_ref,classif_df)
    return rescue_ref.iloc[:,0], rescue_df


def rescue_candidates(classif_df,monoexons,prefix):
    """
    Selection of rescue candidates from non-FSM artifacts.
    The ISM artifacts are selected if they are not associated with a FSM artifact already (they have already been rescued)
    """
    rescue_candidates = identify_rescue_candidates(classif_df, monoexons)
    
    # Write rescue candidates
    rescue_candidates.to_csv(f"{prefix}_rescue_candidates.tsv", 
                            sep="\t",
                            index=False)

    return rescue_candidates["isoform"].tolist()
    
    
def rescue_targets(classif_df,rescue_candidates,ref_gtf,prefix):
    # Get the genes with associated rescue candidates
    target_genes = get_rescue_gene_targets(classif_df,rescue_candidates)

    rescue_targets_lr = classif_df[
        (classif_df['associated_gene'].isin(target_genes)) &
        (classif_df['filter_result'] == 'Isoform')
    ]['isoform']

    rescue_targets_ref = get_rescue_reference_targets(ref_gtf,target_genes)
    rescue_logger.debug(f"Found {rescue_targets_ref.shape[0]} reference targets")
    rescue_logger.debug(f"Found {rescue_targets_lr.shape[0]} long-read targets")
    # Merge both groups and remove duplicates
    rescue_targets = pd.concat([rescue_targets_lr, rescue_targets_ref]).drop_duplicates().reset_index(drop=True)  
    rescue_logger.debug(f"Rescue targets: {rescue_targets.shape[0]}")
    
    rescue_targets.to_csv(f"{prefix}_rescue_targets.tsv",
                        sep="\t",
                        index=False)
    return rescue_targets.tolist()

## Run mapping of rescue candidates (artifacts) to targets
def run_candidate_mapping(ref_trans_fasta,targets_list,candidates_list,
                          corrected_isoforms, out_dir, out_prefix, cpus):
    prefix = f"{out_dir}/{out_prefix}"
    #### PREPARATION OF FILES FOR MINIMAP2 ####
    message("Preparation of files for artifact mapping:", rescue_logger)
    targets_fasta = f"{prefix}_rescue_targets.fasta"
    candidates_fasta = f"{prefix}_rescue_candidates.fasta"

    ## Filter reference transcriptome FASTA to only include target ref transcripts
    rescue_logger.info("Filtering reference transcriptome FASTA to only rescue targets.")
    ref_targets = filter_transcriptome(ref_trans_fasta,targets_list)

    ## Filter SQ3 transcriptome FASTA to only include target LR transcripts
    rescue_logger.info("Filtering supplied long read transcriptome FASTA (--isoforms) to only include rescue targets...")
    LR_targets = filter_transcriptome(corrected_isoforms,targets_list)

    ## join both FASTA files and remove duplicated IDs (some tools like Bambu have the same ids for reference and long read transcripts,
    ## and duplicated IDs within any of the FASTA files would end up as duplicated @SQ lines in the minimap2 SAM header)
    all_targets, dropped_ids = merge_target_records(LR_targets, ref_targets)
    if dropped_ids:
        rescue_logger.warning(f"Removed {len(dropped_ids)} rescue targets with duplicated IDs (long read transcripts are kept over reference ones).")
        rescue_logger.debug(f"Duplicated target IDs: {', '.join(sorted(set(dropped_ids)))}")

    save_fasta(all_targets,targets_fasta)

    ## Filter SQ3 FASTA to include rescue candidates
    rescue_logger.info("Creating rescue candidate FASTA from supplied long read transcriptome fasta (--isoforms)...")

    # make file names
    candidate_filt = filter_transcriptome(corrected_isoforms,candidates_list)
    if len(candidate_filt) == 0:
        rescue_logger.warning("No rescue candidates found in the supplied long read transcriptome FASTA (--corrected_isoforms_fasta).")
        rescue_logger.warning("Are you sure the file is the correct one? It should be the corrected output from SQANTI3 qc, not filter")
    save_fasta(candidate_filt,candidates_fasta)
    
    #### MAPPING ARTIFACTS (CANDIDATES) WITH MINIMAP2 ####
    message("Artifact mapping (candidates vs targets)",rescue_logger)    
    # Mapping
    rescue_logger.info("Mapping rescue candidates to rescue targets with minimap2...")
    # make file names
    sam_file = f"{prefix}_mapped_rescue.sam"
    # The mapping is only reused if it was generated from the same targets, candidates and options
    fingerprint = mapping_fingerprint([targets_fasta, candidates_fasta], MINIMAP2_RESCUE_OPTS)
    if is_mapping_reusable(sam_file, fingerprint):
        rescue_logger.info("Mapping file already exists for the same rescue targets and candidates, skipping mapping step.")
    else:
        if os.path.isfile(sam_file):
            rescue_logger.info("Existing mapping file was generated from different inputs, mapping again.")
        # make command (written to a temporary file so an interrupted run never leaves a truncated SAM behind)
        tmp_sam_file = f"{sam_file}.tmp"
        minimap_cmd = f"minimap2 {MINIMAP2_RESCUE_OPTS} -t {cpus} {targets_fasta} {candidates_fasta} > {tmp_sam_file}"
        # run
        logFile=f"{out_dir}/logs/rescue/minimap2.log"
        run_command(minimap_cmd,rescue_logger,logFile,"Mapping rescue candidates to targets")
        os.replace(tmp_sam_file, sam_file)
        save_mapping_fingerprint(sam_file, fingerprint)
    # Filter mapping results (select SAM columns)
    rescue_logger.info("Building candidate-target table of mapping hits...")
    hits_df = process_sam_file(sam_file)
    # Save as TSV
    hits_file = f"{out_dir}/{out_prefix}_rescue_mapping_hits.tsv"
    hits_df.to_csv(hits_file, sep="\t", index=False, header=True)

    rescue_logger.info(f"Mapping hit table was saved to {hits_file}") 
    rescue_logger.debug("Candidate-target mapping process has been executed successfully.")
    return hits_df

def is_ignored_rule(column, ignored=RESCUE_IGNORED_RULES):
    """Check whether a rule column cannot be evaluated on the reference transcriptome.

    Args:
        column (str): column name used as key in the rules JSON file
        ignored (list): names from RESCUE_IGNORED_RULES (see rule_column_matches)

    Returns:
        bool: True if the rule has to be removed for the reference
    """
    return rule_column_matches(column, ignored)


def write_reference_rules(json_filter, out_file):
    """Write a copy of the rules without the requisites ignored for the reference.

    The structure of the file is kept: every structural category and every rule
    stays, only the ignored requisites are removed. A rule left empty accepts
    every reference transcript of its structural category.

    Args:
        json_filter (str): rules JSON file used to filter the long-read transcriptome
        out_file (str): path of the JSON file to write

    Returns:
        list: sorted names of the removed requisites
    """
    with open(json_filter, 'r') as f:
        rules = json.load(f)

    removed = set()
    reference_rules = {}
    for sc, rule_sets in rules.items():
        reference_rules[sc] = []
        for rule_set in rule_sets:
            kept = {}
            for column, value in rule_set.items():
                if is_ignored_rule(column):
                    removed.add(column)
                else:
                    kept[column] = value
            reference_rules[sc].append(kept)

    with open(out_file, 'w') as f:
        json.dump(reference_rules, f, indent=4)

    return sorted(removed)


## Run rescue steps specific to rules filter
def run_rules_rescue(filter_classification, reference_classification, hits_df, 
                     rescue_df, automatic_inclusion_list, out_dir, json_filter):
    ## Run rules filter on reference transcriptome
    message("Rules rescue selected!",rescue_logger)
    rescue_logger.info("Applying provided rules (--json_filter) to reference transcriptome classification file.")
    ref_out = "reference"
    ref_dir = f"{out_dir}/reference_rules_filter"
    os.makedirs(ref_dir, exist_ok=True)
    # Rules that need long-read evidence cannot be evaluated on the reference
    ref_json = f"{ref_dir}/reference_rules.json"
    removed = write_reference_rules(json_filter, ref_json)
    if removed:
        rescue_logger.info(f"Requisites not applied to the reference transcriptome: {removed}")
    FILTER_PATH = sqanti_path("sqanti3_filter.py")
    # Actual command
    refRules_cmd = f"{PYTHONPATH} {FILTER_PATH} rules --sqanti_class {reference_classification} -j {ref_json} -o {ref_out} -d {ref_dir} --skip_report"
    logFile=f"{out_dir}/logs/refRules.log"
    run_command(refRules_cmd,rescue_logger,logFile,description="Run rules filter on reference transcriptome")

    ## run rescue-by-mapping
    rescue_logger.info("Running rescue-by-mapping for rules filter.")
    # Filenames
    ref_rules = f"{out_dir}/reference_rules_filter/reference_RulesFilter_classification.txt"
    inclusion_list, rescue_df = rescue_by_mapping(hits_df,ref_rules,filter_classification, automatic_inclusion_list,
                                                  rescue_df,"rules")
    return inclusion_list, rescue_df

## Run rescue steps specific to the ML filter
def run_ML_rescue(filter_classification, reference_classification, hits_df, rescue_df,
                  automatic_inclusion_list, out_dir,out_prefix, random_forest, thr):
    prefix = f"{out_dir}/{out_prefix}"
    ## run pre-trained ML classifier on reference transcriptome
    message("ML rescue selected!",rescue_logger)
    rescue_logger.info("Running pre-trained random forest on reference transcriptome classification file.")
    
    # define Rscript command with run_randomforest_on_reference.R args
    refML_cmd = f"{RSCRIPTPATH} {RESCUE_RANDOM_FOREST} -c {reference_classification} -o {out_prefix} -d {out_dir} -r {random_forest}"
    logFile=f"{out_dir}/logs/refML.log"
    run_command(refML_cmd,rescue_logger,logFile,description="Run random forest on reference transcriptome")
    
    ## run rescue-by-mapping
    rescue_logger.info("Running rescue-by-mapping for ML filter.")

    # input file name
    ref_isoform_predict = f"{prefix}_reference_isoform_predict.tsv"

    inclusion_list, rescue_df = rescue_by_mapping(hits_df,ref_isoform_predict,filter_classification, 
                                                  automatic_inclusion_list, rescue_df,"ml",thr)
    
    return inclusion_list, rescue_df

## Evidence check of reference targets (count requisites of the rules)
def run_fallback_mapping(classif_df, rescue_df, ref_trans_fasta, ref_gtf,
                         corrected_isoforms, out_dir, out_prefix, cpus):
    """Map the FSM artifacts of automatic rescue to the rescue targets.

    Automatic rescue does not map its artifacts. Their hits are only used to
    reassign them to a long-read isoform if their reference transcript fails
    the evidence check, so they are kept apart from the rescue-by-mapping hits.
    """
    fsm_artifacts = rescue_df.loc[rescue_df["rescue_mode"] == "automatic", "artifact"].unique().tolist()
    if not fsm_artifacts:
        return pd.DataFrame(columns=["rescue_candidate", "mapping_hit", "alignment_type", "alignment_score"])

    # The hit table is not reused as is: run_candidate_mapping() only reuses the SAM
    # file if it was made from the same targets and artifacts (mapping fingerprint)
    fb_prefix = f"{out_prefix}_fallback"
    message("Fallback mapping of automatic rescue artifacts", rescue_logger)
    targets = rescue_targets(classif_df, fsm_artifacts, ref_gtf, f"{out_dir}/{fb_prefix}")
    return run_candidate_mapping(ref_trans_fasta, targets, fsm_artifacts,
                                 corrected_isoforms, out_dir, fb_prefix, cpus)


def run_evidence_check(classif_df, rescue_df, inclusion_list, hits_df, json_filter,
                       counts_design, counts_file, out_dir, min_expression):
    """Apply the count requisites of the rules to the rescued reference targets.

    The requisites are written to {out_dir}/evidence_check_rules.json and
    evaluated by the rules filter on the counts that requantification would give
    each target from the --counts file, with the expression threshold of QC
    (min_expression).
    """
    message("Evidence check of rescued reference transcripts", rescue_logger)
    with open(json_filter) as f:
        rule_sets = target_count_rules(json.load(f))
    counts_df = load_counts(counts_file, classif_df)
    # The design and the rules were checked against --counts by rescue_args_validation;
    # with a single sample the design is ignored
    design = {}
    if counts_design is not None and len(count_samples(counts_df)) > 1:
        design = read_counts_design(counts_design, rescue_logger)
    rules_dict = write_evidence_rules(rule_sets, f"{out_dir}/evidence_check_rules.json")
    return apply_evidence_check(rescue_df, inclusion_list, classif_df, counts_df, hits_df, rules_dict,
                                design, min_expression)


def concatenate_gtf_files(input_files, output_file):
    """
    Concatenate multiple GTF files into a single GTF file.
    Args:
        input_files (list): List of input GTF file paths.
        output_file (str): Path to the output GTF file.
    """
    with open(output_file, 'w') as outfile:
        for fname in input_files:
            with open(fname) as infile:
                shutil.copyfileobj(infile, outfile)

def save_rescue_results(out_dir,out_prefix, rescued_transcripts, rescue_df, refGTF,
                        filtered_isoforms_gtf,corrected_isoforms_fasta,
                        rescued_class,ref_class):
    prefix = f"{out_dir}/{out_prefix}"
    
    ## Save inclusion list
    inclusion_file = f"{prefix}_rescue_inclusion_list.tsv"
    rescued_transcripts.to_frame(name='isoform').to_csv(inclusion_file,
                                                    sep="\t",
                                                    index=False)
    rescue_logger.info(f"Final inclusion list written to file: {inclusion_file}")

    ## Save rescue table
    rescue_table_file = f"{prefix}_rescue_table.tsv"
    rescue_df.to_csv(rescue_table_file, sep="\t", index=False)
    rescue_logger.info(f"Final rescue table written to file: {rescue_table_file}")
    
    ## Create new GTF including rescued transcripts #
    output_gtf = write_rescue_gtf(filtered_isoforms_gtf, refGTF, rescued_transcripts.to_list(), prefix)
    rescue_logger.info(f"Final output GTF written to file:  {output_gtf}")
    
    ## Create new FASTA including rescued transcripts #
    good_transcripts = rescued_class[rescued_class['filter_result'] == 'Isoform']['isoform']
    ref_fasta_file = os.path.join(out_dir, 
                                  os.path.basename(refGTF).replace('.gtf', '.isoforms.fasta'))
    write_rescue_fasta(corrected_isoforms_fasta,ref_fasta_file, good_transcripts, rescued_transcripts, prefix)
    rescue_logger.info(f"Rescued FASTA written to file: {prefix}_rescued.fasta")

    # Save new classification
    rescued_class = rescued_class[rescued_class['isoform'].isin(good_transcripts)]
    if ref_class is None:
        rescue_logger.warning("No reference classification provided.")
        rescue_logger.warning("Rescued classification will only include the user-defined true isoforms.")
    else:
        rClass = read_classification(ref_class) 
        tClass = rescued_class
        rescued_class = pd.concat([tClass, rClass[rClass['isoform'].isin(rescued_transcripts)]])

    rescued_class.to_csv(f"{prefix}_rescued_classification.txt", sep="\t", index=False)
    rescue_logger.info(f"Rescued classification written to file: {prefix}_rescued_classification.txt")
    return rescued_class
