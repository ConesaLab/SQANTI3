import os
import shutil
import pandas as pd

from src.wrapper_utils import (sqanti_path)
from src.module_logging import rescue_logger, message
from src.commands import (
    RSCRIPTPATH, run_command, PYTHONPATH, RESCUE_RANDOM_FOREST
)
from src.utilities.rescue.automatic_rescue import (
    read_classification, rescue_fsm_monoexons,
    get_lost_reference_id, 
    rescue_lost_reference, save_automatic_rescue
)
from src.utilities.rescue.rescue_helpers import (
    get_rescue_gene_targets,get_rescue_reference_targets
)

from src.utilities.rescue.candidate_mapping_helpers import (
    filter_transcriptome,
    prepare_fasta_transcriptome,
    process_sam_file,
    save_fasta
)

from src.utilities.rescue.rescue_by_mapping import rescue_by_mapping


def run_automatic_rescue(classification_file,monoexons,prefix):
    # Load classification
    message("Reading filter classification file",rescue_logger)
    classif_df = read_classification(classification_file)

    message("Performing automatic rescue",rescue_logger)
    # Select the FSM and ISM isoforms with more than one exon 
    rescue_classif = classif_df[
        (classif_df['structural_category'].isin(['full-splice_match'])) & 
        (classif_df['exons'] > 1)
    ]
 
    # Find the references that are lost and get the ones that are not represented by isoforms
    lost_ref = get_lost_reference_id(rescue_classif)
    if len(lost_ref) == 0:
       rescue_logger.info("No lost references found")
       rescue_logger.info("Automatic rescue is not needed")
       save_automatic_rescue(pd.DataFrame({'associated_transcript': ['none']}),classif_df,prefix)
       return
    rescue_logger.debug(f"Found {len(lost_ref)} lost references")
    rescue = pd.DataFrame()
    for ref_id in lost_ref:
        rescue = pd.concat([rescue,rescue_lost_reference(ref_id, rescue_classif)])

    # Split into reference transcripts and ISM TODO: Eliminate this step?
    rescue_ref = rescue[rescue['isoform'].isin(rescue_classif['associated_transcript'])]
    # Adding monoexons
    if monoexons in ['all','fsm']:
        rescue_fsm_me = rescue_fsm_monoexons(classif_df)
        rescue_auto = pd.concat([rescue_ref,rescue_fsm_me])
    else:
        rescue_auto = rescue_ref
    rescue_logger.debug(f"Rescued {rescue_auto.shape[0]} transcripts")
    
    # Save the automatic rescue
    save_automatic_rescue(rescue_auto,classif_df,prefix)
    return

def rescue_candidates(classification_file,monoexons,prefix):
    """
    Selection of rescue candidates from non-FSM artifacts.
    The ISM artifacts are selected if they are not associated with a FSM artifact already (they have already been rescued)
    """
    # Load classification
    classif_df = read_classification(classification_file)
    
    # Identify transcripts that have a full-splice_match Artifact
    transcripts_with_fsm_artifact = set(
        classif_df[
            (classif_df['structural_category'] == 'full-splice_match') &
            (classif_df['filter_result'] == 'Artifact')
        ]['associated_transcript']
    )

    # Get initial set of candidates
    rescue_candidates = classif_df[
        (classif_df['structural_category'].isin(['incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog'])) &
        (classif_df['filter_result'] == 'Artifact')
    ]

    # Exclude incomplete-splice_match with a conflicting full-splice_match Artifact
    rescue_candidates = rescue_candidates[
        ~(
            (rescue_candidates['structural_category'] == 'incomplete-splice_match') &
            (rescue_candidates['associated_transcript'].isin(transcripts_with_fsm_artifact))
        )
    ]

    if monoexons != 'all':
        rescue_novel = rescue_novel[rescue_novel['exons'] > 1]
    
    # Write rescue candidates
    rescue_candidates.to_csv(f"{prefix}_rescue_candidates.tsv", 
                            sep="\t",
                            index=False)

    return rescue_candidates["isoform"].tolist()
    
def rescue_targets(classification_file,rescue_candidates,ref_gtf,prefix):
    # Load classification
    classif_df = read_classification(classification_file)
    
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
                          corrected_isoforms, out_dir, out_prefix):
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

    ## join both FASTA files
    all_targets = ref_targets + LR_targets
    save_fasta(all_targets,targets_fasta)

    ## Filter SQ3 FASTA to include rescue candidates
    rescue_logger.info("Creating rescue candidate FASTA from supplied long read transcriptome fasta (--isoforms)...")

    # make file names
    candidate_filt = filter_transcriptome(corrected_isoforms,candidates_list)
    save_fasta(candidate_filt,candidates_fasta)
    
    #### MAPPING ARTIFACTS (CANDIDATES) WITH MINIMAP2 ####
    # TODO: eliminate file logic in the process
    message("Artifact mapping (candidates vs targets)",rescue_logger)
    
    # Mapping
    rescue_logger.info("Mapping rescue candidates to rescue targets with minimap2...")

    # make file names
    sam_file = f"{prefix}_mapped_rescue.sam"
    if os.path.isfile(sam_file):
        rescue_logger.info("Mapping file already exists, skipping mapping step.")
    else:
        # make command
        minimap_cmd = f"minimap2 --secondary=yes -ax map-hifi {targets_fasta} {candidates_fasta} > {sam_file}"

        # run
        logFile=f"{out_dir}/logs/rescue/minimap2.log"
        run_command(minimap_cmd,rescue_logger,logFile,"Mapping rescue candidates to targets")

    # Filter mapping results (select SAM columns)
    rescue_logger.info("Building candidate-target table of mapping hits...")

    process_sam_file(sam_file,out_dir,out_prefix)

    rescue_logger.debug("Candidate-target mapping process has been executed successfully.")


## Run rescue steps specific to rules filter
def run_rules_rescue(filter_classification, reference_classification,
                     out_dir, out_prefix, json_filter):
    prefix = f"{out_dir}/{out_prefix}"
    ## Run rules filter on reference transcriptome
    message("Rules rescue selected!",rescue_logger)
    rescue_logger.info("Applying provided rules (--json_filter) to reference transcriptome classification file.")

    # create reference out prefix and dir
    ref_out = "reference"
    ref_dir = f"{out_dir}/reference_rules_filter"
    FILTER_PATH = sqanti_path("sqanti3_filter.py")
    # define command
    refRules_cmd = f"{PYTHONPATH} {FILTER_PATH} rules --sqanti_class {reference_classification} -j {json_filter} -o {ref_out} -d {ref_dir} --skip_report"
    logFile=f"{out_dir}/logs/refRules.log"
    run_command(refRules_cmd,rescue_logger,logFile,description="Run rules filter on reference transcriptome")

    ## run rescue-by-mapping
    rescue_logger.info("Running rescue-by-mapping for rules filter.")
    # Filenames
    mapping_hits = f"{prefix}_rescue_mapping_hits.tsv"
    ref_rules = f"{out_dir}/reference_rules_filter/reference_RulesFilter_result_classification.txt"
    rescue_by_mapping(mapping_hits,ref_rules,filter_classification,
                      f"{prefix}_automatic_rescue_table.tsv",prefix,"rules")


## Run rescue steps specific to the ML filter
def run_ML_rescue(filter_classification, reference_classification,
                  out_dir,out_prefix, random_forest, thr):
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
    mapping_hits = f"{prefix}_rescue_mapping_hits.tsv"
    ref_isoform_predict = f"{prefix}_reference_isoform_predict.tsv"

    rescue_by_mapping(mapping_hits,ref_isoform_predict,filter_classification,
                        f"{prefix}_automatic_rescue_table.tsv",prefix,"ml",thr)

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

def save_rescue_results(out_dir,out_prefix, mode, refGTF,
                        filtered_isoforms_gtf,corrected_isoforms_fasta,
                        filter_class,ref_class):
    prefix = f"{out_dir}/{out_prefix}"
    # create file names
    tmp_gtf = f"{out_dir}/rescued_only_tmp.gtf"
    output_gtf = f"{prefix}_rescued.gtf"

    # Select the propper inclusion list
    if mode == "full":
        rescued_list = f"{prefix}_full_inclusion_list.tsv"
    else:
        rescued_list = f"{prefix}_automatic_inclusion_list.tsv"
    # Read the rescued transcripts from the inclusion list
    rescued_transcripts = set()
    with open(rescued_list, 'r') as f:
      for line in f:
        rescued_transcripts.add(line.strip())
    f.close()    
    output_gtf = write_rescue_gtf(filtered_isoforms_gtf, refGTF, rescued_transcripts, prefix)
    rescue_logger.info(f"Final output GTF written to file:  {prefix}_rescued.gtf")
    
    ## Create new FASTA including rescued transcripts #
    good_transcripts = get_good_transcripts(filter_class)
    ref_fasta_file = os.path.join(out_dir, 
                                  os.path.basename(refGTF).replace('.gtf', '.fasta'))
    write_rescue_fasta(corrected_isoforms_fasta,ref_fasta_file, good_transcripts, rescued_transcripts, prefix)
    rescue_logger.info(f"Rescued FASTA written to file: {prefix}_rescued.fasta")

    # Save new classification
    if ref_class is None:
        rescue_logger.warning("No reference classification provided.")
        rescue_logger.warning("Rescued classification will only include the user-defined isoforms.")
        rescued_class = read_classification(filter_class)
    else:
        rClass = read_classification(ref_class) 
        tClass = read_classification(filter_class)
        rescued_class = pd.concat([tClass[tClass['isoform'].isin(good_transcripts)],
                                    rClass[rClass['isoform'].isin(rescued_transcripts)]])
    rescued_class.to_csv(f"{prefix}_rescued_classification.tsv", sep="\t", index=False)
    rescue_logger.info(f"Rescued classification written to file: {prefix}_rescued_classification.tsv")
    return(rescued_list,output_gtf)
