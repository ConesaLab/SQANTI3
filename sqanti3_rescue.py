#!/usr/bin/env python3
__author__  = "angeles.arzalluz@gmail.com"

###################################################
##########     SQANTI3 RESCUE WRAPPER    ##########
###################################################

#### PREPARATION ####

## Module import
import os
import sys

from src.rescue_argparse import rescue_argparse
from src.module_logging import rescue_logger, message, update_logger
from src.logging_config import rescue_art, art_logger

from src.argparse_utils import rescue_args_validation
from src.config import __version__
from src.rescue_output import write_rescue_fasta, write_rescue_gtf
from src.rescue_steps import (
  concatenate_gtf_files,
  run_automatic_rescue,
  rescue_candidates, rescue_targets,
  run_candidate_mapping, run_rules_rescue, run_ML_rescue
)
from src.utilities.rescue import sq_requant
from src.utilities.rescue.rescue_helpers import get_good_transcripts


def main():
  art_logger.info(rescue_art())
  args = rescue_argparse().parse_args()
  update_logger(rescue_logger,args.dir,args.log_level)
  # Check if the logs directory exists, if not create it
  rescue_args_validation(args)
  rescue_logger.info(f"Running SQANTI3 rescue pipeline version {__version__}")

  #### RUN AUTOMATIC RESCUE ####
  # this part is run for both rules and ML and if all arg tests passed
  message(f"Initializing SQANTI3 rescue pipeline in {args.mode} mode",rescue_logger)
  prefix = f"{args.dir}/{args.output}"
  # TODO: Pre-read the filter_class into a pd dataframe and pass it to the functions (once we have pythonzed everything)
  run_automatic_rescue(args.filter_class,args.rescue_mono_exonic,
                              args.mode,prefix)
  message("Automatic rescue completed",rescue_logger)

  ### RUN FULL RESCUE (IF REQUESTED) ###
  if args.mode == "full":
    candidates = rescue_candidates(args.filter_class,args.rescue_mono_exonic,
                                   prefix)
    rescue_logger.debug(f"Rescue candidates: {len(candidates)}")

    targets = rescue_targets(args.filter_class,candidates,
                             args.refGTF,prefix)
    rescue_logger.debug(f"Rescue targets: {len(targets)}")    
    
    #### RUN MAPPING
    # when in full mode, rescue maps candidates not included in the
    # automatic rescue (ISM, NIC, NNC) to long-read and reference
    # isoforms passing the filter (targets)

    run_candidate_mapping(args.refGTF,args.refFasta,targets,candidates,
                          args.corrected_isoforms_fasta,args.dir,args.output)


    #### RUN ML FILTER RESCUE ####
    # this part combines reference ML filter run with mapping results
    # and is therefore run only for ML filter

    if args.strategy == "ml":

      message("Rescue-by-mapping for ML filter",rescue_logger)
      # run ML-specific steps of rescue
      rescued = run_ML_rescue(args.filter_class, args.refClassif,
                              args.dir, args.output, args.random_forest, args.threshold)


    #### RUN RULES FILTER RESCUE ####
    # this part runs SQ3 rules filter for the reference transcriptome
    # and combines the results with the mapping hits obtained in the previous step
    if args.strategy == "rules":
      message("Rescue-by-mapping for rules filter", rescue_logger)
      # run rules-specific steps of rescue
      rescued = run_rules_rescue(args.filter_class, args.refClassif,
                                  args.dir, args.output, args.json_filter)


    # Finish print if output exists (same for rules and ML) ####
    inclusion_list = f"{prefix}_full_inclusion_list.tsv"
    if os.path.isfile(inclusion_list):
      message(f"Rescue {args.strategy} finished successfully!",rescue_logger)
      rescue_logger.info(f"Final rescued transcript list written to file: {inclusion_list}")
    else:
      rescue_logger.error(f"Something went wrong, inclusion list not found: {inclusion_list}")
      sys.exit(1)
  ### End of condition (mode == "full")

  #### WRITE FINAL OUTPUTS OF RESCUE ####
  # Create new GTF including rescued transcripts #
  if args.filtered_isoforms_gtf is None:
    rescue_logger.warning("No filtered GTF provided.")
    rescue_logger.warning("Rescue will be performed but no GTF will be generated.")
  else:
    message("Generating rescued GTF.",rescue_logger)

    # Select the propper inclusion list
    if args.mode == "full":
        rescued_list = f"{prefix}_full_inclusion_list.tsv"
    else:
        rescued_list = f"{prefix}_automatic_inclusion_list.tsv"
    # Read the rescued transcripts from the inclusion list
    rescued_transcripts = set()
    with open(rescued_list, 'r') as f:
      for line in f:
        rescued_transcripts.add(line.strip())
    f.close()    
    write_rescue_gtf(args.filtered_isoforms_gtf, args.refGTF, rescued_transcripts, prefix)
    rescue_logger.info(f"Final output GTF written to file:  {prefix}_rescued.gtf")
    
    ## Create new FASTA including rescued transcripts #
    good_transcripts = get_good_transcripts(args.filter_class)
    ref_fasta_file = os.path.join(args.dir, 
                                  os.path.basename(args.refGTF).replace('.gtf', '.fasta'))
    write_rescue_fasta(args.corrected_isoforms_fasta,ref_fasta_file, good_transcripts, rescued_transcripts, prefix)
    rescue_logger.info(f"Rescued FASTA written to file: {prefix}_rescued.fasta")
  ## END ##
  message("Rescue finished successfully!",rescue_logger)

  if args.requant:
    if not args.counts:
      rescue_logger.error("Counts file is required for the requantification module.")
    message("\nRunning requantification.",rescue_logger)

    rescue_gtf, inclusion_list, counts, rescued = sq_requant.parse_files(args)
    sq_requant.run_requant(rescue_gtf, inclusion_list, counts, rescued, args.dir, args.output)
    sq_requant.to_tpm(rescue_gtf, args.dir, args.output)
    rescue_logger.info("\nRequantification finished!\n")

## Run main()
if __name__ == "__main__":
    main()
