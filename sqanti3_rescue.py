#!/usr/bin/env python3
__author__  = "angeles.arzalluz@gmail.com"

###################################################
##########     SQANTI3 RESCUE WRAPPER    ##########
###################################################

#### PREPARATION ####

## Module import
import os, sys, shutil
import pandas as pd

from src.rescue_argparse import rescue_argparse
from src.module_logging import rescue_logger, message
from src.logging_config import rescue_art, art_logger
from src.argparse_utils import rescue_args_validation
from src.commands import run_command, utilitiesPath
from src.config import __version__
from src.rescue_steps import (
  run_automatic_rescue_py,
  rescue_candidates, rescue_targets,
  run_candidate_mapping, run_rules_rescue, run_ML_rescue
)

## Set general path variables
Rscript_path = shutil.which('Rscript')
gffread_path = shutil.which('gffread')
python_path = shutil.which('python')

## Set path variables to call R scripts

## Set path variables to call SQ3 scripts
filter_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "sqanti3_filter.py")

#### MAIN ####
## Define main()
def main():
  art_logger.info(rescue_art())
  args = rescue_argparse().parse_args()
  rescue_args_validation(args)
  rescue_logger.info(f"Running SQANTI3 rescue pipeline version {__version__}")
  #### RUN AUTOMATIC RESCUE ####
  # this part is run for both rules and ML and if all arg tests passed
  message(f"Initializing SQANTI3 rescue pipeline in {args.mode} mode",rescue_logger)
  prefix = f"{args.dir}/{args.output}"
  #run_automatic_rescue(args)
  rescue_ism = run_automatic_rescue_py(args.filter_class,args.rescue_mono_exonic,
                              args.mode,prefix)

  ### RUN FULL RESCUE (IF REQUESTED) ###
  if args.mode == "full":

    candidates = rescue_candidates(args.filter_class,args.rescue_mono_exonic,
                                   rescue_ism,prefix)
    rescue_logger.debug(f"Rescue candidates: {len(candidates)}")

    targets = rescue_targets(args.filter_class,candidates,
                             args.refGTF,prefix)
    rescue_logger.debug(f"Rescue targets: {len(targets)}")    
    
    #### RUN MAPPING
    # when in full mode, rescue maps candidates not included in the
    # automatic rescue (ISM, NIC, NNC) to long-read and reference
    # isoforms passing the filter (targets)

    run_candidate_mapping(args,targets,candidates)
    

    #### RUN ML FILTER RESCUE ####
    # this part combines reference ML filter run with mapping results
    # and is therefore run only for ML filter

    if args.subcommand == "ml":

      message("Rescue-by-mapping for ML filter",rescue_logger)

      # run ML-specific steps of rescue
      rescued = run_ML_rescue(args)


    #### RUN RULES FILTER RESCUE ####
    # this part runs SQ3 rules filter for the reference transcriptome
    # and combines the results with the mapping hits obtained in the previous step

    if args.subcommand == "rules":

      rescue_logger.info("**** RESCUE-BY-MAPPING FOR RULES FILTER")

      # run rules-specific steps of rescue
      rescued = run_rules_rescue(args)


    # Finish print if output exists (same for rules and ML) ####
    inclusion_list = f"{args.dir}/{args.output}_rescue_inclusion-list.tsv"

    if os.path.isfile(inclusion_list):
      rescue_logger.info(f"Final rescued transcript list witten to file: {inclusion_list}")

  ### End of condition (mode == "full")



  #### WRITE FINAL OUPTUTS OF RESCUE ####
  # Create new GTF including rescued transcripts #

  rescue_logger.info("Adding rescued transcripts to provided SQ3 filtered GTF.")

  # create file names
  tmp_gtf = f"{args.dir}/rescued_only_tmp.gtf"
  output_gtf = f"{args.dir}/{args.output}_rescued.gtf"

  # condition: inclusion list file only produced for mode = full
  # in mode = automatic, it is replaced by automatic rescue list
  if args.mode == "full":
      rescued_list = f"{args.dir}/{args.output}_rescue_inclusion-list.tsv"
  else:
      rescued_list = f"{args.dir}/{args.output}_automatic_rescued_list.tsv"

  # filter reference GTF to create tmp_gtf
  gtf_cmd = f"gffread --ids {rescued_list} -T -o {tmp_gtf} {args.refGTF}"

  run_command(gtf_cmd,rescue_logger,"log/rescue/gtf.log",description="Filter reference GTF to create tmp GTF")
  # concatenate with filtered GTF
  cat_cmd = f"cat {args.rescue_gtf} {tmp_gtf} > {output_gtf}"
  run_command(cat_cmd,rescue_logger,"log/rescue/cat.log",description="Concatenate filtered GTF with tmp GTF")
  
  rescue_logger.info(f"Added rescued reference transcripts to provided GTF ({args.rescue_gtf} )")
  rescue_logger.info(f"Final output GTF written to file:  {output_gtf} ")

  # remove tmp_gtf
  os.remove(tmp_gtf)

  ## END ##
  rescue_logger.info(f"Rescue finished successfully!")



## Run main()
if __name__ == "__main__":
    main()
