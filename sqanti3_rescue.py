#!/usr/bin/env python3
__author__  = "angeles.arzalluz@gmail.com"
__version__ = "5.1"

###################################################
##########     SQANTI3 RESCUE WRAPPER    ##########
###################################################

#### PREPARATION ####

## Module import
import os, sys, argparse, subprocess
import distutils.spawn

## Set general path variables
Rscript_path = distutils.spawn.find_executable('Rscript')
gffread_path = distutiles.spawn.find_executable('gffread')
utilities_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")

## Set path variables to call R scripts
automatic_rescue_path = "rescue/automatic_rescue.R"

## Check that Rscript is working
if os.system(Rscript_path + " --version") != 0:
    print("Rscript executable not found. Abort!", file = sys.stderr)
    sys.exit(-1)
    
## Check that gffread is working
if os.system(gffread_path + " --version") != 0:
  print("Cannot find gffread executable. Abort!", file = sys.stderr)
    
    
#### DEFINE FUNCTIONS ####

## Run automatic rescue
def run_automatic_rescue(args):
  
  ## prepare to run script
  
  # define Rscript command with automatic_rescue.R args
  auto_cmd = Rscript_path + " {u}/{s} -c {c} -o {o} -d {d} -u {u} \
  -g {g} -e {e}".format(u = utilities_path, s = automatic_rescue_path, \
  c = args.sqanti_MLclassif, o = args.output, d = args.dir, \
  g = args.refGTF, e = args.rescue_mono_exonic)
  
  # print command
  print(auto_cmd + "\n")
  
  ## run automatic rescue script via terminal
  subprocess.call(auto_cmd, shell = True)
  
  ## load outputs
  # transcripts rescued as a result of automatic rescue
  automatic_rescued_list = args.dir + "/" + args.output + "_automatic_rescued_list.tsv"
  # set object containing rescued list from the output file
  auto_rescue = set(line.strip() for line in open(automatic_rescued_list))
  
  # rescue candidate list
  rescue_candidate_list = args.dir + "/" + args.output + "_rescue_candidates.tsv"
  rescue_candidates = set(line.strip() for line in open(rescue_candidate_list))
  
  # rescue target list
  rescue_target_list = args.dir + "/" + args.output + "_rescue_targets.tsv"
  rescue_targets = set(line.strip() for line in open(rescue_target_list))

  # return automatic rescue outputs
  return(auto_rescue, rescue_candidates, rescue_targets)



#### MAIN ####

## Define main()
def main():
  
  ## Arguments and help
  parser = argparse.ArgumentParser(description = "Rescue artifacts discarded by \
  the SQANTI3 filter, i.e. find closest match for the artifacts in the reference \
  transcriptome and add them to the transcriptome.")
  
  ## Common arguments
  common = argparse.ArgumentParser(add_help = False)
  common.add_argument("sqanti_MLclassif", \
  help = "\t\tSQANTI ML output classification file.")
  common.add_argument("--isoforms", \
  help = "\t\tFASTA file output by SQANTI3 QC (*_corrected.fa).")
  common.add_argument("-g", "--refGTF", \
  help = "\t\tFull path to reference transcriptome GTF used when running SQANTI3 QC.")
  common.add_argument("-f", "--refGenome", \
  help = "\t\tFull path to reference genome FASTA used when running SQANTI3 QC.")
  common.add_argument("-e","--rescue_mono_exonic", \
  choices = ['all', 'fsm', 'none'], default = "all", \
  help='\t\tWhether or not to include mono-exonic artifacts in the rescue. Options include: none, fsm and all (default).')
  common.add_argument("-o","--output", \
  help = "\t\tPrefix for output files.", required = False)
  common.add_argument("-d","--dir", \
  help = "\t\tDirectory for output files. Default: Directory where the script was run.", \
  required = False)
  common.add_argument("--skip_report", action = "store_true", default = False, \
  help = '\t\tSkip creation of a report about the filtering')
  common.add_argument("-v", "--version", help="Display program version number.", \
  action='version', version='SQANTI3 '+str(__version__))

  subparsers = parser.add_subparsers(dest = 'subcommand')

  ## ML rescue arguments
  ml = subparsers.add_parser("ml", parents = [common], description = "Rescue for ML-filtered transcriptomes.")
  
  ml.add_argument("-r", "--randomforest", \
  help = "Full path to the randomforest.RData object obtained when running the SQANTI3 ML filter.")
  ml.add_argument("-j", "--threshold", type = float, default = 0.7, \
  help = "Default: 0.7. Machine learning probability threshold to filter elegible rescue targets (mapping hits).")
  
  ## Rules rescue arguments
  rules = subparsers.add_parser("rules", parents = [common], \
  description = "Rescue for rules-filtered transcriptomes.")
  
  # parse arguments
  args = parser.parse_args()
  
  
  ## Check that common arguments are valid
  args.sqanti_MLclassif = os.path.abspath(args.sqanti_MLclassif)
  if not os.path.isfile(args.sqanti_MLclassif):
    print("ERROR: {0} doesn't exist. Abort!".format(args.sqanti_MLclassif), file=sys.stderr)
    sys.exit(-1)
    
  if not os.path.isfile(args.isoforms):
    print("ERROR: {0} doesn't exist. Abort!".format(args.isoforms), file=sys.stderr)
    sys.exit(-1)
    
  if not os.path.isfile(args.refGTF):
    print("ERROR: {0} doesn't exist. Abort!".format(args.refGTF), file=sys.stderr)
    sys.exit(-1)
    
  if not os.path.isfile(args.refGenome):
    print("ERROR: {0} doesn't exist. Abort!".format(args.refGenome), file=sys.stderr)
    sys.exit(-1)
  
  ## Check that ML-specific args are valid
  if args.subcommand == "ML":
    if not os.path.isfile(args.randomforest):
    print("ERROR: {0} doesn't exist. Abort!".format(args.randomforest), file=sys.stderr)
    sys.exit(-1)
    
    if args.threshold < 0 or args.threshold > 1.:
      print("ERROR: --threshold must be between 0-1, value of {0} was supplied! Abort!".format(args.threshold), file=sys.stderr)
      sys.exit(-1)
  
  
  ## Run automatic rescue 
  # this part is run for both rules and ML and if all arg tests passed
  print("\nRunning automatic rescue...\n", file = sys.stdout)
  auto_result, candidates, targets = run_automatic_rescue(args)
  
  
  ## Convert reference transcriptome GTF to FASTA
  # make FASTA file name
  pre, ext = os.path.splitext(args.refGTF)
  refGTF_rename = pre + ".fasta"
  
  # build gffread command
  ref_cmd = "gffread -w {w} -g {g} {a}".format(w = refGTF_rename, g = args.refGenome, \
  a = args.refGTF)
  
  print(ref_cmd, "\n")
  
  # run gffread
  subprocess.call(ref_cmd, shell = True)
  
  


## Run main()
if __name__ == "__main__":
    main()
