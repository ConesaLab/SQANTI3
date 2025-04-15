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
  run_candidate_mapping
)

## Set general path variables
Rscript_path = shutil.which('Rscript')
gffread_path = shutil.which('gffread')
python_path = shutil.which('python')

## Set path variables to call R scripts
run_randomforest_path = "rescue/run_randomforest_on_reference.R"
rescue_by_mapping_ML_path = "rescue/rescue_by_mapping_ML.R"
rescue_by_mapping_rules_path = "rescue/rescue_by_mapping_rules.R"

## Set path variables to call SQ3 scripts
filter_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "sqanti3_filter.py")

# ## Check that Rscript is working
# if os.system(Rscript_path + " --version") != 0:
#     rescue_logger.error("Rscript executable not found. Abort!")
#     sys.exit(1)

# ## Check that gffread is working
# if os.system(gffread_path + " --version") != 0:
#   rescue_logger.error("Cannot find gffread executable. Abort!")
#   sys.exit(1)


#### DEFINE FUNCTIONS ####

## Run rescue steps specific to the ML filter
def run_ML_rescue(args):

  ## run pre-trained ML classifier on reference transcriptome

  rescue_logger.info("ML rescue selected!")
  rescue_logger.info("Running pre-trained random forest on reference transcriptome classification file.")

  # define Rscript command with run_randomforest_on_reference.R args
  refML_cmd = f"{Rscript_path} {utilitiesPath}/{run_randomforest_path} -c {args.refClassif} -o {args.output} -d {args.dir} -r {args.randomforest}"
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
    rescue_cmd = f"{Rscript_path} {utilitiesPath}/{rescue_by_mapping_ML_path} -c {args.filter_class} -o {args.output} -d {args.dir} -u {utilitiesPath} -m {mapping_hits} -r {ref_isoform_predict} -j {args.threshold}"


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



## Run rescue steps specific to rules filter
def run_rules_rescue(args):

  ## Run rules filter on reference transcriptome

  rescue_logger.info("**** Rules rescue selected!")
  rescue_logger.info("Applying provided rules (--json_filter) to reference transcriptome classification file.")

  # create reference out prefix and dir
  ref_out = "reference"
  ref_dir = f"{args.dir}/reference_rules_filter"

  # define command
  refRules_cmd = f"{python_path} {filter_path} rules {args.refClassif} -j {args.json_filter} -o {ref_out} -d {ref_dir}"


  # print command
  run_command(refRules_cmd,rescue_logger,"log/rescue/refRules.log",description="Run rules filter on reference transcriptome")
    # make file names
  ref_rules = f"{args.dir}/reference_rules_filter/reference_RulesFilter_result_classification.txt"

  if os.path.isfile(ref_rules):
    ## run rescue-by-mapping
    rescue_logger.info("Running rescue-by-mapping for rules filter.")

    # input file name
    mapping_hits = f"{args.dir}/{args.output}_rescue_mapping_hits.tsv"

    # define Rscsript command with rescue_by_mapping_ML.R args
    rescue_cmd = f"{Rscript_path} {utilitiesPath}/{rescue_by_mapping_rules_path} -c {args.filter_class} \
      -o {args.output} -d {args.dir} -u {utilitiesPath} -m {mapping_hits} -r {ref_rules}"


    # expected output name
    rescued_file = f"{args.dir}/{args.output}_rescue_inclusion-list.tsv"
    run_command(rescue_cmd,rescue_logger,"log/rescue/rescue.log",description="Run rescue by mapping")
  
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
    rescue_logger.debug(f"Rescue candidates: {candidates.shape[0]}")

    targets = rescue_targets(args.filter_class,candidates,
                             args.refGTF,prefix)
    rescue_logger.debug(f"Rescue targets: {targets.shape[0]}")    
    
    #### RUN MAPPING
    # when in full mode, rescue maps candidates not included in the
    # automatic rescue (ISM, NIC, NNC) to long-read and reference
    # isoforms passing the filter (targets)

    run_candidate_mapping(args,targets.tolist(),candidates.tolist())
    

    #### RUN ML FILTER RESCUE ####
    # this part combines reference ML filter run with mapping results
    # and is therefore run only for ML filter

    if args.subcommand == "ml":

      rescue_logger.info("**** RESCUE-BY-MAPPING FOR ML FILTER")

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
