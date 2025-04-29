import os
import sys

import pandas as pd

from src.utilities.rescue.rescue_by_mapping_rules import rescue_rules
from src.wrapper_utils import (sqanti_path)
from src.module_logging import rescue_logger, message
from src.commands import (
    RSCRIPTPATH, RESCUE_AUTO_PATH, utilitiesPath, run_command ,
    PYTHONPATH, RSCRIPT_RESCUE_RULES, RSCRIPT_RESCUE_ML, RESCUE_RANDOM_FOREST
)
from src.utilities.rescue.automatic_rescue import (
    read_classification, rescue_fsm_monoexons, add_ism_monoexons,
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

def run_automatic_rescue(classification_file,monoexons,mode,prefix):
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
       return
    rescue_logger.debug(f"Found {len(lost_ref)} lost references")
    rescue = pd.DataFrame()
    for ref_id in lost_ref:
        rescue_df = rescue_lost_reference(ref_id, rescue_classif)
        rescue = pd.concat([rescue,rescue_df])

    # Split into reference transcripts and ISM
    rescue_ref = rescue[rescue['isoform'].isin(rescue_classif['associated_transcript'])]
    # Adding monoexons
    if monoexons in ['all','fsm']:
        rescue_fsm_me = rescue_fsm_monoexons(classif_df)
        rescue_auto = pd.concat([rescue_ref,rescue_fsm_me])
    else:
        rescue_auto = rescue_ref
    rescue_logger.debug(f"Rescued {rescue_auto.shape[0]} transcripts")

    # Save the automatic rescue
    save_automatic_rescue(rescue_auto,rescue_classif,mode,prefix)

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
def run_candidate_mapping(args,targets_list,candidates_list):
    
    #### PREPARATION OF FILES FOR MINIMAP2 ####
    message("Preparation of files for artifact mapping:", rescue_logger)
    targets_fasta = f"{args.dir}/{args.output}_rescue_targets.fasta"
    candidates_fasta = f"{args.dir}/{args.output}_rescue_candidates.fasta"

    ## Convert reference transcriptome GTF to FASTA
    ref_trans_fasta = prepare_fasta_transcriptome(args.refGTF,args.refFasta,args.dir)

    ## Filter reference transcriptome FASTA to only include target ref transcripts
    rescue_logger.info("Filtering reference transcriptome FASTA to only rescue targets.")

    # make file names
    ref_targets = filter_transcriptome(ref_trans_fasta,targets_list)

    ## Filter SQ3 transcriptome FASTA to only include target LR transcripts
    rescue_logger.info("Filtering supplied long read transcriptome FASTA (--isoforms) to only include rescue targets...")

    # make file names
    LR_targets = filter_transcriptome(args.rescue_isoforms,targets_list)

    ## join both FASTA files
    all_targets = ref_targets + LR_targets
    save_fasta(all_targets,targets_fasta)

    ## Filter SQ3 FASTA to include rescue candidates
    rescue_logger.info("Creating rescue candidate FASTA from supplied long read transcriptome fasta (--isoforms)...")

    # make file names
    candidate_filt = filter_transcriptome(args.rescue_isoforms,candidates_list)
    save_fasta(candidate_filt,candidates_fasta)
    
    #### MAPPING ARTIFACTS (CANDIDATES) WITH MINIMAP2 ####
    # TODO: eliminate file logic in the process
    message("Artifact mapping (candidates vs targets)",rescue_logger)
    
    # Mapping
    rescue_logger.info("Mapping rescue candidates to rescue targets with minimap2...")

    # make file names
    sam_file = f"{args.dir}/{args.output}_mapped_rescue.sam"

    # make command
    minimap_cmd = f"minimap2 --secondary=yes -ax map-hifi {targets_fasta} {candidates_fasta} > {sam_file}"

    # run
    logFile=f"{args.dir}/logs/rescue/minimap2.log"
    run_command(minimap_cmd,rescue_logger,logFile,"Mapping rescue candidates to targets")

    if os.path.isfile(sam_file):
        rescue_logger.info(f"Minimap2 results were saved to {sam_file}")
        rescue_logger.debug("minimap2 command used:")
        rescue_logger.debug(minimap_cmd)

    # Filter mapping results (select SAM columns)
    rescue_logger.info("Building candidate-target table of mapping hits...")

    process_sam_file(sam_file,args.dir,args.output)

    rescue_logger.debug("Candidate-target mapping process has been executed successfully.")


## Run rescue steps specific to rules filter
def run_rules_rescue(args):

    ## Run rules filter on reference transcriptome

    message("Rules rescue selected",rescue_logger)
    rescue_logger.info("Applying provided rules (--json_filter) to reference transcriptome classification file.")

    # create reference out prefix and dir
    ref_out = "reference"
    ref_dir = f"{args.dir}/reference_rules_filter"

    FILTER_PATH = sqanti_path("sqanti3_filter.py")
    # define command
    refRules_cmd = f"{PYTHONPATH} {FILTER_PATH} rules --sqanti_class {args.refClassif} -j {args.json_filter} -o {ref_out} -d {ref_dir} --skip_report"

    # print command
    run_command(refRules_cmd,rescue_logger,"log/rescue/refRules.log",description="Run rules filter on reference transcriptome")
        # make file names
    ref_rules = f"{args.dir}/reference_rules_filter/reference_RulesFilter_result_classification.txt"

    if os.path.isfile(ref_rules):
        ## run rescue-by-mapping
        rescue_logger.info("Running rescue-by-mapping for rules filter.")

        # input file name
        mapping_hits = f"{args.dir}/{args.output}_rescue_mapping_hits.tsv"

        # define Rscript command with rescue_by_mapping_rules.R args
        rescue_cmd = f"{RSCRIPTPATH} {RSCRIPT_RESCUE_RULES} -c {args.filter_class} \
        -o {args.output} -d {args.dir} -u {utilitiesPath} -m {mapping_hits} -r {ref_rules}"


        # expected output name
        rescued_file = f"{args.dir}/{args.output}_rescue_inclusion-list.tsv"
        automatic_rescue_file = f"{args.dir}/{args.output}_automatic_rescue_table.tsv"
        run_command(rescue_cmd,rescue_logger,"log/rescue/rescue.log",description="Run rescue by mapping")
        # TODO: Find a way to run this part in python 
        # print(mapping_hits, ref_rules, args.filter_class, automatic_rescue_file, f"{args.dir}/{args.output}")
        # rescue_rules(mapping_hits, ref_rules, args.filter_class, automatic_rescue_file, f"{args.dir}/{args.output}")
        if os.path.isfile(rescued_file):
            # load output list of rescued transcripts
            rescued_df = pd.read_table(rescued_file, header = None, \
            names = ["transcript"])
            rescued_list = list(rescued_df["transcript"])

            # return rescued transcript list
            return(rescued_list)

        else:
            rescue_logger.error("ERROR: rescue inclusion list not created -file not found!")
            sys.exit(1)

    else:
        rescue_logger.error("ERROR: reference filter classification not found!")
        sys.exit(1)


## Run rescue steps specific to the ML filter
def run_ML_rescue(args):

  ## run pre-trained ML classifier on reference transcriptome
  message("ML rescue selected!",rescue_logger)
  rescue_logger.info("Running pre-trained random forest on reference transcriptome classification file.")

  # define Rscript command with run_randomforest_on_reference.R args
  refML_cmd = f"{RSCRIPTPATH} {RESCUE_RANDOM_FOREST} -c {args.refClassif} -o {args.output} -d {args.dir} -r {args.random_forest}"

  # print command
  rescue_logger.debug(refML_cmd)

  # run R script via terminal
  run_command(refML_cmd,rescue_logger,"log/rescue/refML.log",description="Run random forest on reference transcriptome")
  # make expected output file name
  ref_isoform_predict = f"{args.dir}/{args.output}_reference_isoform_predict.tsv"

  if os.path.isfile(ref_isoform_predict):

    ## run rescue-by-mapping
    rescue_logger.info("Running rescue-by-mapping for ML filter.")

    # input file name
    mapping_hits = f"{args.dir}/{args.output}_rescue_mapping_hits.tsv"

    # define Rscsript command with rescue_by_mapping_ML.R args
    rescue_cmd = f"{RSCRIPTPATH} {RSCRIPT_RESCUE_ML} -c {args.filter_class} -o {args.output} -d {args.dir} -u {utilitiesPath} -m {mapping_hits} -r {ref_isoform_predict} -j {args.threshold}"

    # expected output name
    rescued_file = f"{args.dir}/{args.output }_rescue_inclusion-list.tsv"

    # run R script via terminal
    run_command(rescue_cmd,rescue_logger,"log/rescue/rescue.log",description="Run rescue by mapping")

    if os.path.isfile(rescued_file):
      # load output list of rescued transcripts
      rescued_df = pd.read_table(rescued_file, header = None, \
      names = ["transcript"])
      rescued_list = list(rescued_df["transcript"])

      # return rescued transcript list
      return(rescued_list)

    else:
      rescue_logger.error("Rescue inclusion list not created -file not found!")
      sys.exit(1)

  else:
    rescue_logger.error("Reference isoform predictions not found!")
    sys.exit(1)