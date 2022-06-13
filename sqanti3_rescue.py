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
gffread_path = distutils.spawn.find_executable('gffread')
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
  print("\nAutomatic rescue run via the following command:\n")
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
  auto_result, candidates, targets = run_automatic_rescue(args)
  

  #### PREPARATION OF FILES FOR MINIMAP2 ####
  
  print("\n-------------------------------------------------------\n")
  print("\n\tPREPARATION OF FILES FOR ARTIFACT MAPPING:\n")
  print("\n-------------------------------------------------------\n")

  ## Convert reference transcriptome GTF to FASTA
  
  print("\nCreating reference transcriptome FASTA from provided GTF (--refGTF)...\n")

  # make FASTA file name
  pre, ext = os.path.splitext(args.refGTF)
  refGTF_rename = pre + ".fasta"
  
  # build gffread command
  ref_cmd = "gffread -w {w} -g {g} {a}".format(w = refGTF_rename, g = args.refGenome, \
  a = args.refGTF)
  
  print("\n\tgffread command used:\n")
  print(ref_cmd, "\n")
  
  # run gffread
  subprocess.call(ref_cmd, shell = True)
  
  
  ## Filter reference transcriptome FASTA to only include target ref transcripts
  
  print("\nFiltering reference transcriptome FASTA to only rescue targets...\n")

  # make file names
  target_file =	args.dir + "/" + args.output + "_rescue_targets.tsv"
  ref_target_fasta = args.dir +	"/" + args.output + "_rescue_targets.ref.fasta"

  # make command
  fasta_cmd = "seqtk subseq {i}	{t} > {f}".format(i = refGTF_rename, \
  t = target_file, f = ref_target_fasta)

  print("\n\tseqtk command used:\n")
  print(fasta_cmd, "\n")

  # run	
  subprocess.call(fasta_cmd, shell = True)


  ## Filter SQ3	transcriptome FASTA to only include target LR transcripts
  
  print("\nFiltering supplied long read transcriptome FASTA (--isoforms) to only include rescue targets...\n")

  # make file names
  LR_target_fasta = args.dir +  "/" + args.output + "_rescue_targets.LR.fasta"

  # make command
  fasta_cmd = "seqtk subseq {i} {t} > {f}".format(i = args.isoforms, \
  t = target_file, f = LR_target_fasta)

  print("\n\tseqtk command used:\n")
  print(fasta_cmd, "\n")

  # run
  subprocess.call(fasta_cmd, shell = True)  


  ## join both FASTA files

  print("\nJoining reference and LR rescue target FASTA files...\n")

  target_fasta = args.dir +  "/" + args.output + "_rescue_targets.fasta"
  cat_cmd = "cat {r} {l} > {f}".format(r = ref_target_fasta, l = LR_target_fasta, \
  f = target_fasta)

  subprocess.call(cat_cmd, shell = True)

  print("\nRescue target FASTA was saved to ", target_fasta, "\n")
  print("\nCommand used:")
  print(cat_cmd, "\n")


  ## remove intermediate target FASTA files (LR and ref)
  
  print("\nRemoving intermediate target FASTA files...\n")

  rm_cmd = "rm {r} {l}".format(r = ref_target_fasta, l = LR_target_fasta)

  subprocess.call(rm_cmd, shell = True)


  ## Filter SQ3 FASTA to include rescue candidates

  print("\nCreating rescue candidate FASTA from supplied long read transcriptome fasta (--isoforms)...\n ")

  # make file names
  candidate_file = args.dir + "/" + args.output + "_rescue_candidates.tsv"
  candidate_fasta = args.dir + "/" + args.output + "_rescue_candidates.fasta"

  # make command
  fasta_cmd = "seqtk subseq {i} {t} > {f}".format(i = args.isoforms, \
  t = candidate_file, f = candidate_fasta)

  print(fasta_cmd, "\n")

  # run
  subprocess.call(fasta_cmd, shell = True)

  print("\nRescue candidate FASTA was saved to ", candidate_fasta, "\n")
  print("\n\tseqtk command used:\n")
  print(fasta_cmd, "\n")


  #### MAPPING ARTIFACTS (CANDIDATES) WITH MINIMAP2 ####

  print("\n-------------------------------------------------------\n")
  print("\n\tARTIFACT MAPPING (CANDIDATES VS TARGETS):\n")
  print("\n-------------------------------------------------------\n")

  ## Mapping

  print("\nMapping rescue candidates to rescue targets with minimap2...\n")

  # make file names
  sam_file = args.dir + "/" + args.output + "_mapped_rescue.sam"

  # make command
  minimap_cmd = "minimap2 --secondary=yes -ax map-hifi {t} {c} > {s}".format( \
  t = target_fasta, c = candidate_fasta, s = sam_file)

  # run
  subprocess.call(minimap_cmd, shell = True)

  print("\nMinimap2 results were saved to ", sam_file, "\n")
  print("\n\tminimap2 command used:\n")
  print(minimap_cmd, "\n")


  ## Filter mapping results (select SAM columns)

  print("\nBuilding candidate-target table of mapping hits...\n")

  # remove header from SAM
  sam_tmp_file = args.dir + "/" + args.output + "_mapped_rescue_noheader.sam"
  sam_cmd = "grep -v '@' {s} > {t}".format(s = sam_file, t = sam_tmp_file)
  subprocess.call(sam_cmd, shell = True)

  # get cols with candidate-target pairs + alignment type
  hits_file = args.dir + "/" + args.output + "_rescue_mapping_hits.tsv"
  hits_cmd = "cut -f1-3 {t} > {h}".format(t = sam_tmp_file, h = hits_file)
  subprocess.call(hits_cmd, shell = True)

  print("\nMapping hit table was saved to ", hits_file, "\n")

  # delete altered SAM file
  rm_cmd = "rm {t}".format(t = sam_tmp_file)
  subprocess.call(rm_cmd, shell = True)



## Run main()
if __name__ == "__main__":
    main()
