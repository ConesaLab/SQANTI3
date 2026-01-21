#!/usr/bin/env python3
__author__  = "angeles.arzalluz@gmail.com"

###################################################
##########     SQANTI3 RESCUE WRAPPER    ##########
###################################################

#### PREPARATION ####

## Module import
import os

import pandas as pd

from src.parsers import parse_counts
from src.rescue_argparse import rescue_argparse
from src.module_logging import rescue_logger, message, update_logger
from src.logging_config import rescue_art, art_logger

from src.argparse_utils import rescue_args_validation
from src.config import __version__
from src.rescue_steps import (
  run_automatic_rescue,
  rescue_candidates, rescue_targets,
  run_candidate_mapping, run_rules_rescue, run_ML_rescue,
  save_rescue_results
)
from src.utilities.rescue.candidate_mapping_helpers import prepare_fasta_transcriptome
from src.utilities.rescue.rescue_helpers import read_classification
from src.utilities.rescue.sq_requant import parse_files, run_requant, to_tpm

def main():
  art_logger.info(rescue_art())
  args = rescue_argparse().parse_args()
  update_logger(rescue_logger,args.dir,"rescue",args.log_level)
  # Create the log directory if it does not exist
  os.makedirs(f"{args.dir}/logs", exist_ok=True)
  # Check if the logs directory exists, if not create it
  rescue_args_validation(args)
  rescue_logger.info(f"Running SQANTI3 rescue pipeline version {__version__}")

  #### RUN AUTOMATIC RESCUE ####
  # this part is run for both rules and ML and if all arg tests passed
  message(f"Initializing SQANTI3 rescue pipeline in {args.mode} mode",rescue_logger)
  prefix = f"{args.dir}/{args.output}"
  # Load classification
  message("Reading filter classification file",rescue_logger)
  class_df = read_classification(args.filter_class)
  ## AUTOMATIC RESCUE ##
  inclusion_list, rescue_df = run_automatic_rescue(class_df,args.rescue_mono_exonic,prefix)
  message("Automatic rescue completed",rescue_logger)
  ## Convert reference transcriptome GTF to FASTA
  ref_trans_fasta = prepare_fasta_transcriptome(args.refGTF,args.refFasta,args.dir)

  ### RUN FULL RESCUE (IF REQUESTED) ###
  if args.mode == "full":
    candidates = rescue_candidates(class_df,args.rescue_mono_exonic,
                                   prefix)
    rescue_logger.debug(f"Rescue candidates: {len(candidates)}")

    targets = rescue_targets(class_df,candidates,
                             args.refGTF,prefix)
    rescue_logger.debug(f"Rescue targets: {len(targets)}")    
    
    #### RUN MAPPING
    # when in full mode, rescue maps candidates not included in the
    # automatic rescue (ISM, NIC, NNC) to long-read and reference
    # isoforms passing the filter (targets)

    if os.path.isfile(f"{prefix}_rescue_mapping_hits.tsv"):
      rescue_logger.info("Mapping hits already exist, skipping mapping step.")
      hits_df = pd.read_csv(f"{prefix}_rescue_mapping_hits.tsv",sep="\t")
    else:
      hits_df = run_candidate_mapping(ref_trans_fasta,targets,candidates,
                           args.corrected_isoforms_fasta,args.dir,args.output)

    #### RUN ML FILTER RESCUE ####
    # this part combines reference ML filter run with mapping results
    # and is therefore run only for ML filter
    if args.strategy == "ml":
      message("Rescue-by-mapping for ML filter",rescue_logger)
      # run ML-specific steps of rescue
      inclusion_list, rescue_df = run_ML_rescue(class_df, args.refClassif, hits_df, rescue_df,
                                                inclusion_list, args.dir, args.output, 
                                                args.random_forest, args.threshold)

    #### RUN RULES FILTER RESCUE ####
    # this part runs SQ3 rules filter for the reference transcriptome
    # and combines the results with the mapping hits obtained in the previous step
    if args.subcommand == "rules":
      message("Rescue-by-mapping for rules filter", rescue_logger)
      # run rules-specific steps of rescue
      inclusion_list, rescue_df = run_rules_rescue(class_df, args.refClassif, hits_df, rescue_df,
                                                   inclusion_list,args.dir, args.json_filter)

  #### WRITE FINAL OUTPUTS OF RESCUE ####
  # Create new GTF including rescued transcripts #
  if args.rescue_gtf is None:
    rescue_logger.warning("No filtered GTF provided.")
    rescue_logger.warning("Rescue will be performed but no GTF will be generated.")
  else:
    message("Generating rescued GTF.",rescue_logger)
    rescue_class = save_rescue_results(args.dir, args.output, inclusion_list, rescue_df,
                                       args.refGTF, args.filtered_isoforms_gtf,args.corrected_isoforms_fasta,
                                       class_df,args.refClassif)

  ## END ##
  message("Rescue finished successfully!",rescue_logger)

  if args.requant:  
    message("Running requantification.",rescue_logger)
    #TODO: Make this take the variables from python directly
    counts_df = parse_counts(args.counts)
    rescue_logger.info("Counts file parsed.")
    requant_df = run_requant(counts_df, rescue_df, class_df, prefix)
    rescue_logger.info("Requantification of counts completed.")
    rescue_logger.info(f"New count table saved to {prefix}_reassigned_counts.tsv")
    # Doing this, we loose the counts assigned to multi_transcript and artifacts (they have no length, so TPM cannot be calculated)
    to_tpm(requant_df,rescue_class, prefix)
    rescue_logger.info("Requantification finished!")

## Run main()
if __name__ == "__main__":
    main()
