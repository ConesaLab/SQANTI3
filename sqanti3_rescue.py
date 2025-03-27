#!/usr/bin/env python3
__author__  = "angeles.arzalluz@gmail.com"
import shutil
from src.commands import run_command
from src.config import __version__
###################################################
##########     SQANTI3 RESCUE WRAPPER    ##########
###################################################

#### PREPARATION ####

## Module import
import os, sys, subprocess
import pandas as pd

from src.rescue_argparse import rescue_argparse
from src.logging_config import rescue_art, art_logger
from src.module_logging import rescue_logger
## Set general path variables
Rscript_path = shutil.which('Rscript')
gffread_path = shutil.which('gffread')
python_path = shutil.which('python')
utilities_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "src/utilities")

## Set path variables to call R scripts
automatic_rescue_path = "rescue/automatic_rescue.R"
run_randomforest_path = "rescue/run_randomforest_on_reference.R"
rescue_by_mapping_ML_path = "rescue/rescue_by_mapping_ML.R"
rescue_by_mapping_rules_path = "rescue/rescue_by_mapping_rules.R"

## Set path variables to call SQ3 scripts
filter_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "sqanti3_filter.py")

## Check that Rscript is working
if os.system(Rscript_path + " --version") != 0:
    rescue_logger.error("Rscript executable not found. Abort!")
    sys.exit(1)

## Check that gffread is working
if os.system(gffread_path + " --version") != 0:
  rescue_logger.error("Cannot find gffread executable. Abort!")
  sys.exit(1)


#### DEFINE FUNCTIONS ####

## Run automatic rescue
def run_automatic_rescue(args):

  ## prepare to run script

  # define Rscript command with automatic_rescue.R args
  auto_cmd = Rscript_path + " {u}/{s} -c {c} -o {o} -d {d} -u {u} \
  -g {g} -e {e} -m {m}".format(u = utilities_path, s = automatic_rescue_path, \
  c = args.filter_class, o = args.output, d = args.dir, \
  g = args.refGTF, e = args.rescue_mono_exonic, m = args.mode)

  # print command
  rescue_logger.info("Automatic rescue run via the following command:")
  rescue_logger.info(auto_cmd)

  ## run automatic rescue script via terminal
  if subprocess.check_call(auto_cmd, shell = True) != 0:
    rescue_logger.error("ERROR running automatic rescue: {auto_cmd}")
    sys.exit(1)

  ## load output: transcripts rescued as a result of automatic rescue

  # make file name
  automatic_rescued_list = f"{args.dir}/{args.output}_automatic_rescued_list.tsv"

  # set object containing rescued list from the output file
  auto_rescue = set(line.strip() for line in open(automatic_rescued_list))

  ## return automatic rescue outputs
  return(auto_rescue)


## Run mapping of rescue candidates (artifacts) to targets
def run_candidate_mapping(args):

  #### PREPARATION OF FILES FOR MINIMAP2 ####
  rescue_logger.info("**** Preparation of files for artifact mapping:")
  ## Convert reference transcriptome GTF to FASTA

  rescue_logger.info("Creating reference transcriptome FASTA from provided GTF (--refGTF).")

  # make FASTA file name
  pre, _ = os.path.splitext(os.path.basename(args.refGTF))
  ref_trans_Fasta = f" {args.dir}/{pre}.fasta"

  # build gffread command
  ref_cmd = f"gffread -w {ref_trans_Fasta} -g {args.refFasta} {args.refGTF}"

  # run gffread
  try:
      subprocess.check_call(ref_cmd, shell=True)
      if os.path.isfile(ref_trans_Fasta):
          rescue_logger.info(f"Reference transcriptome FASTA was saved to {ref_trans_Fasta}")
          rescue_logger.info("gffread command used:")
          rescue_logger.info(ref_cmd)
      else:
          rescue_logger.error("Reference transcriptome FASTA was not created - file not found!")
          sys.exit(1)
  except subprocess.CalledProcessError:
      rescue_logger.error(f"Error converting reference transcriptome GTF to FASTA. Command used: {ref_cmd}")
      sys.exit(1)



  ## Filter reference transcriptome FASTA to only include target ref transcripts

  rescue_logger.info("Filtering reference transcriptome FASTA to only rescue targets.")

  # make file names
  target_file = f"{args.dir}/{args.output}_rescue_targets.tsv"
  ref_target_fasta = f"{args.dir}/{args.output}_rescue_targets.ref.fasta"

  # make command
  fasta_cmd = f"seqtk subseq {ref_trans_Fasta} {target_file} > {ref_target_fasta}"

  # run
  try:
      subprocess.check_call(fasta_cmd, shell=True)
      if os.path.isfile(ref_target_fasta):
          rescue_logger.info(f"Target reference transcript sequences were saved to {ref_target_fasta}")
          rescue_logger.info("seqtk command used:")
          rescue_logger.info(fasta_cmd)
      else:
          rescue_logger.error("Target reference transcript FASTA was not created - file not found!")
          sys.exit(1)
  except subprocess.CalledProcessError:
      rescue_logger.error(f"Error retrieving target reference transcripts from FASTA. Command used: {fasta_cmd}")
      sys.exit(1)


  ## Filter SQ3	transcriptome FASTA to only include target LR transcripts

  rescue_logger.info("Filtering supplied long read transcriptome FASTA (--isoforms) to only include rescue targets...")

  # make file names
  LR_target_fasta = f"{args.dir}/{args.output}_rescue_targets.LR.fasta"

  # make command
  fasta_cmd = f"seqtk subseq {args.rescue_isoforms} {target_file} > {LR_target_fasta}"

  # run
  try:
      subprocess.check_call(fasta_cmd, shell=True)
      if os.path.isfile(LR_target_fasta):
          rescue_logger.info(f"Target long read transcript sequences were saved to {LR_target_fasta}")
          rescue_logger.info("seqtk command used:")
          rescue_logger.info(fasta_cmd)
      else:
          rescue_logger.error("Target long read transcript FASTA was not created - file not found!")
          sys.exit(1)
  except subprocess.CalledProcessError:
      rescue_logger.error(f"Error retrieving target long-read transcripts from FASTA. Command used: {fasta_cmd}")
      sys.exit(1)


  ## join both FASTA files
  rescue_logger.info("Joining reference and LR rescue target FASTA files...")

  target_fasta = f"{args.dir}/{args.output}_rescue_targets.fasta"
  cat_cmd = f"cat {ref_target_fasta} {LR_target_fasta} > {target_fasta}"

  try:
      subprocess.check_call(cat_cmd, shell=True)
      if os.path.isfile(target_fasta):
          rescue_logger.info(f"Rescue target FASTA was saved to {target_fasta}")
          rescue_logger.info("Command used:")
          rescue_logger.info(cat_cmd)

          # Remove intermediate target FASTA files (LR and ref)
          rescue_logger.info("Removing intermediate target FASTA files...")
          rm_cmd = f"rm {ref_target_fasta} {LR_target_fasta}"
          try:
              subprocess.call(rm_cmd, shell=True)
          except subprocess.CalledProcessError:
              rescue_logger.warning(f"Failed to remove intermediate files. Command used: {rm_cmd}")
      else:
          rescue_logger.error("Target FASTA was not created - file not found!")
          sys.exit(1)
  except subprocess.CalledProcessError:
      rescue_logger.error(f"Error joining target long-read and reference FASTA files. Command used: {cat_cmd}")
      sys.exit(1)


  ## Filter SQ3 FASTA to include rescue candidates
  rescue_logger.info("Creating rescue candidate FASTA from supplied long read transcriptome fasta (--isoforms)...")

  # make file names
  candidate_file = f"{args.dir}/{args.output}_rescue_candidates.tsv"
  candidate_fasta = f"{args.dir}/{args.output}_rescue_candidates.fasta"

  # make command
  fasta_cmd = f"seqtk subseq {args.rescue_isoforms} {candidate_file} > {candidate_fasta}"

  # run
  try:
      subprocess.check_call(fasta_cmd, shell=True)
      if os.path.isfile(candidate_fasta):
          rescue_logger.info(f"Rescue candidate FASTA was saved to {candidate_fasta}")
          rescue_logger.info("seqtk command used:")
          rescue_logger.info(fasta_cmd)
      else:
          rescue_logger.error("Candidate FASTA was not created - file not found!")
          sys.exit(1)
  except subprocess.CalledProcessError:
      rescue_logger.error(f"Error retrieving rescue candidate sequences from FASTA. Command used: {fasta_cmd}")
      sys.exit(1)


  #### MAPPING ARTIFACTS (CANDIDATES) WITH MINIMAP2 ####
  # TODO: eliminate file logic in the process
  rescue_logger.info("**** Artifact mapping (candidates vs targets)")
  
  # Mapping
  rescue_logger.info("Mapping rescue candidates to rescue targets with minimap2...")

  # make file names
  sam_file = f"{args.dir}/{args.output}_mapped_rescue.sam"

  # make command
  minimap_cmd = f"minimap2 --secondary=yes -ax map-hifi {target_fasta} {candidate_fasta} > {sam_file}"

  # run
  logFile=f"{args.dir}/logs/rescue/minimap2.log"
  run_command(minimap_cmd,rescue_logger,logFile,"Maping rescue candidates to targets")

  if os.path.isfile(sam_file):
    rescue_logger.info(f"Minimap2 results were saved to {sam_file}")
    rescue_logger.info("minimap2 command used:")
    rescue_logger.info(minimap_cmd)

    # Filter mapping results (select SAM columns)
    rescue_logger.info("Building candidate-target table of mapping hits...")

    # remove header from SAM
    sam_tmp_file = f"{args.dir}/{args.output}_mapped_rescue_noheader.sam"
    sam_cmd = f"grep -v '@' {sam_file} > {sam_tmp_file}"

    run_command(sam_cmd,rescue_logger,"log/rescue/sam_noheader.log",description="Remove header from SAM file")
    if os.path.isfile(sam_tmp_file):
      # get cols with candidate-target pairs + alignment type
      hits_file = f"{args.dir}/{args.output}_rescue_mapping_hits.tsv"
      hits_cmd = f"cut -f1-3 {sam_tmp_file} > {hits_file}"

      run_command(hits_cmd,rescue_logger,"log/rescue/hits.log",description="Extract candidate-target pairs from SAM file")

      if os.path.isfile(hits_file):
        rescue_logger.info(f"Mapping hit table was saved to {hits_file}")

        # delete altered SAM file
        os.remove(rm_cmd)


## Run rescue steps specific to the ML filter
def run_ML_rescue(args):

  ## run pre-trained ML classifier on reference transcriptome

  rescue_logger.info("ML rescue selected!")
  rescue_logger.info("Running pre-trained random forest on reference transcriptome classification file.")

  # define Rscript command with run_randomforest_on_reference.R args
  refML_cmd = f"{Rscript_path} {utilities_path}/{run_randomforest_path} -c {args.refClassif} -o {args.output} -d {args.dir} -r {args.randomforest}"
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
    rescue_cmd = f"{Rscript_path} {utilities_path}/{rescue_by_mapping_ML_path} -c {args.filter_class} -o {args.output} -d {args.dir} -u {utilities_path} -m {mapping_hits} -r {ref_isoform_predict} -j {args.threshold}"


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
  rescue_logger.info("Applying provided rules (--json) to reference transcriptome classification file.")

  # create reference out prefix and dir
  ref_out = "reference"
  ref_dir = f"{args.dir}/reference_rules_filter"

  # define command
  refRules_cmd = f"{python_path} {filter_path} rules {args.refClassif} -j {args.json_filter} -o {ref_out} -d {ref_dir}"


  # print command
  rescue_logger.debug(refRules_cmd)
  run_command(refRules_cmd,rescue_logger,"log/rescue/refRules.log",description="Run rules filter on reference transcriptome")
    # make file names
  ref_rules = f"{args.dir}/reference_rules_filter/reference_RulesFilter_result_classification.txt"

  if os.path.isfile(ref_rules):
    ## run rescue-by-mapping
    rescue_logger.info("Running rescue-by-mapping for rules filter.")

    # input file name
    mapping_hits = f"{args.dir}/{args.output}_rescue_mapping_hits.tsv"

    # define Rscsript command with rescue_by_mapping_ML.R args
    rescue_cmd = f"{Rscript_path} {utilities_path}/{rescue_by_mapping_rules_path} -c {args.filter_class} \
      -o {args.output} -d {args.dir} -u {utilities_path} -m {mapping_hits} -r {ref_rules}"


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

  ## Check that common arguments are valid
  args.filter_class = os.path.abspath(args.filter_class)
  if not os.path.isfile(args.filter_class):
      rescue_logger.error(f"{args.filter_class} doesn't exist. Abort!")
      sys.exit(1)

  if not os.path.isfile(args.rescue_isoforms):
      rescue_logger.error(f"{args.rescue_isoforms} doesn't exist. Abort!")
      sys.exit(1)

  if not os.path.isfile(args.rescue_gtf):
      rescue_logger.error(f"{args.rescue_gtf} doesn't exist. Abort!")
      sys.exit(1)

  if not os.path.isfile(args.refGTF):
      rescue_logger.error(f"{args.refGTF} doesn't exist. Abort!")
      sys.exit(1)

  if not os.path.isfile(args.refFasta):
      rescue_logger.error(f"{args.refFasta} doesn't exist. Abort!")
      sys.exit(1)

  # TODO: Condition to run SQANTI_QC on the reference transcriptome
  if not os.path.isfile(args.refClassif):
      rescue_logger.error(f"{args.refClassif} doesn't exist. Abort!")
      sys.exit(1)

  ## Define output dir and output name in case it was not defined
  if args.dir is None:
      args.dir = os.path.dirname(args.filter_class)
      rescue_logger.warning(f"Output directory not defined. All the outputs will be stored at {args.dir} directory")
  else:
      if not os.path.exists(args.dir):
          os.makedirs(args.dir)
          rescue_logger.info(f"Created output directory: {args.dir}")

  ## Check that ML-specific args are valid
  if args.subcommand == "ml":
      if not os.path.isfile(args.randomforest):
          rescue_logger.error(f"{args.randomforest} doesn't exist. Abort!")
          sys.exit(1)

      if args.threshold < 0 or args.threshold > 1.:
          rescue_logger.error(f"--threshold must be between 0-1, value of {args.threshold} was supplied. Abort!")
          sys.exit(1)

  ## Check that rules-specific args are valid
  if args.subcommand == "rules":
      if not os.path.isfile(args.json_filter):
          rescue_logger.error(f"{args.json_filter} doesn't exist. Abort!")
          sys.exit(1)



  #### RUN AUTOMATIC RESCUE ####
  # this part is run for both rules and ML and if all arg tests passed

  auto_result = run_automatic_rescue(args)


  ### RUN FULL RESCUE (IF REQUESTED) ###
  if args.mode == "full":

    #### RUN MAPPING
    # when in full mode, rescue maps candidates not included in the
    # automatic rescue (ISM, NIC, NNC) to long-read and reference
    # isoforms passing the filter (targets)

    run_candidate_mapping(args)


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
