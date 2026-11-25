#!/usr/bin/env python
__author__  = "angeles.arzalluz@gmail.com"

###################################################
##########     SQANTI3 RESCUE WRAPPER    ##########
###################################################

#### PREPARATION ####

## Module import
import os

import pandas as pd

from src.rescue_argparse import rescue_argparse
from src.module_logging import rescue_logger, message, update_logger
from src.logging_config import rescue_art, art_logger

from src.argparse_utils import rescue_args_validation
from src.config import __version__
from src.rescue_steps import (
  run_automatic_rescue,
  rescue_candidates, rescue_targets,
  run_candidate_mapping, run_rules_rescue, run_ML_rescue,
  save_rescue_results,
  run_fallback_mapping, run_evidence_check
)
from src.utilities.rescue.candidate_mapping_helpers import prepare_fasta_transcriptome
from src.utilities.rescue.rescue_helpers import read_classification
from src.utilities.rescue.sq_requant import requantification_pipeline
from src.utilities.rescue.evidence_check import rows_for_requant
from src.write_parameters import write_rescue_parameters

def main():
  art_logger.info(rescue_art())
  args = rescue_argparse().parse_args()
  update_logger(rescue_logger,args.dir,"rescue",args.log_level)
  # Create the log directory if it does not exist
  os.makedirs(f"{args.dir}/logs", exist_ok=True)
  # Check if the logs directory exists, if not create it
  evidence_check = rescue_args_validation(args)
  write_rescue_parameters(args)
  rescue_logger.info(f"Running SQANTI3 rescue pipeline version {__version__}")

  #### RUN AUTOMATIC RESCUE ####
  # this part is run for both rules and ML and if all arg tests passed
  message(f"Initializing SQANTI3 rescue pipeline in {args.mode} mode",rescue_logger)
  prefix = f"{args.dir}/{args.output}"
  # Load classification
  message("Reading filter classification file",rescue_logger)
  class_df = read_classification(args.filter_class)
  ## AUTOMATIC RESCUE ##
  inclusion_list, rescue_df = run_automatic_rescue(class_df,args.rescue_mono_exonic)
  message("Automatic rescue completed",rescue_logger)
  ## Convert reference transcriptome GTF to FASTA
  ref_trans_fasta = prepare_fasta_transcriptome(args.refGTF,args.refFasta,args.dir)

  hits_df = pd.DataFrame(columns=["rescue_candidate", "mapping_hit", "alignment_type", "alignment_score"])
  ### RUN FULL RESCUE (IF REQUESTED) ###
  candidates = None
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

    # The hit table is not reused as is: run_candidate_mapping() only reuses the SAM
    # file if it was made from the same targets and candidates (mapping fingerprint)
    hits_df = run_candidate_mapping(ref_trans_fasta,targets,candidates,
                         args.corrected_isoforms_fasta,args.dir,args.output,
                         args.cpus)

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
    elif args.strategy == "rules":
      message("Rescue-by-mapping for rules filter", rescue_logger)
      # run rules-specific steps of rescue
      inclusion_list, rescue_df = run_rules_rescue(class_df, args.refClassif, hits_df, rescue_df,
                                                   inclusion_list,args.dir, args.json_filter)

  #### EVIDENCE CHECK OF REFERENCE TARGETS ####
  # The requisites on long-read counts, removed to filter the reference, are applied
  # to the counts each reintroduced reference transcript receives from its artifacts
  if evidence_check:
    # FSM artifacts of automatic rescue are not mapped: if their reference fails, their
    # counts go to the gene residual. With --map_automatic_fsm they are mapped, so that
    # they can fall back to a long-read isoform (in any rescue mode)
    fallback_hits = hits_df
    if args.map_automatic_fsm:
      fsm_hits = run_fallback_mapping(class_df, rescue_df, ref_trans_fasta, args.refGTF,
                                      args.corrected_isoforms_fasta, args.dir, args.output,
                                      args.cpus)
      fallback_hits = pd.concat([hits_df, fsm_hits], ignore_index=True)
    inclusion_list, rescue_df = run_evidence_check(class_df, rescue_df, inclusion_list, fallback_hits,
                                                   args.json_filter, args.counts_design, args.counts, args.dir,
                                                   args.min_expression)

  #### WRITE FINAL OUTPUTS OF RESCUE ####
  # Create new GTF including rescued transcripts #
  if args.filtered_isoforms_gtf is None:
    rescue_class = None
    rescue_logger.warning("No filtered GTF provided.")
    rescue_logger.warning("Rescue will be performed but no GTF will be generated.")
  else:
    message("Generating rescued GTF.",rescue_logger)
    rescue_class = save_rescue_results(args.dir, args.output, inclusion_list, rescue_df,
                                       args.refGTF, args.filtered_isoforms_gtf,args.corrected_isoforms_fasta,
                                       class_df,args.refClassif)

  ## END ##
 

  if not args.skip_requant:
    if rescue_class is None:
      rescue_logger.warning(
        "Requantification skipped: it requires the rescued classification, which is only "
        "generated when --filtered_isoforms_gtf is provided."
      )
    else:
      message("Running requantification.", rescue_logger)
      requantification_pipeline(args.dir, args.output, args.counts, rows_for_requant(rescue_df),
                                class_df, rescue_class, args.counts_design, args.min_expression)
  else:
    rescue_logger.info("Requantification skipped (--skip_requant). The count matrix will not include rescued transcripts.")

  message("Rescue finished successfully!", rescue_logger)

## Run main()
if __name__ == "__main__":
    main()
