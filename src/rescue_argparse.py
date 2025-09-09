import argparse
from src.config import __version__, default_json

def rescue_argparse():

# Arguments and help
  ## Common arguments
  common = argparse.ArgumentParser(add_help = False)
  cr = common.add_argument_group("Required arguments")
  cr.add_argument("--filter_class", 
                  required=True,
                  help = "SQANTI filter (ML or rules) output classification file.")
  cr.add_argument("--refGTF",
                  required=True,
                  help = "Full path to reference transcriptome GTF used when running SQANTI3 QC.")
  cr.add_argument("--refFasta", 
                  required=True,
                  help = "Full path to reference genome FASTA used when running SQANTI3 QC.")

  # Specific input options
  ci = common.add_argument_group("Input options")
  ci.add_argument("--rescue_isoforms", 
                      help = "FASTA file output by SQANTI3 QC (*_corrected.fasta), i.e. the full long read transcriptome.")
  ci.add_argument("--rescue_gtf", 
                      help = "GTF file output by SQANTI3 filter (*.filtered.gtf).")
  ci.add_argument("-k", "--refClassif", 
                      help = "Full path to the classification file obtained when running SQANTI3 QC on the reference transcriptome.\
                        \nMandatory when running the rescue on full mode")
  ci.add_argument("--counts",
                      help = 'Isoforms abundancy values: "Isoform" \t "Count". Column names may differ')

  # Customization options
  cc = common.add_argument_group("Customization options")
  cc.add_argument("-e","--rescue_mono_exonic", 
                  choices = ['all', 'fsm', 'none'], 
                  default = "all", 
                  help='Whether or not to include mono-exonic artifacts in the rescue.\
                    \nDefault: %(default)s')
  cc.add_argument("--mode", 
                  choices = ["automatic", "full"],
                  default = "automatic", 
                  help = "If 'automatic' (default), only automatic rescue of FSM artifacts will be performed.\
                     \nIf 'full', rescue will include mapping of ISM, NNC and NIC artifacts to find potential replacement isoforms.")
  cc.add_argument("-q","--requant",
                  action="store_true",
                  help = "Run requantification of the rescued isoforms.")
  # Output options
  co = common.add_argument_group("Output options")
  co.add_argument("-o","--output", 
                      help = "Prefix for output files.", 
                      required = False)
  co.add_argument("-d","--dir", 
                      help = "Directory for output files. Default: Directory where the script was run.", 
                      required = False)
  # Performance options
  cp = common.add_argument_group("Extra options")
  cp.add_argument("-c", "--cpus",
                      type = int, default = 4, 
                      help = "Number of CPUs to use. Default: 4")
  cp.add_argument("-v", "--version", 
                      help="Display program version number.", 
                      action='version', 
                      version=f"SQANTI3 v{__version__}")
  cp.add_argument("-l","--log_level", default="INFO",choices=["ERROR","WARNING","INFO","DEBUG"],
                      help="Set the logging level %(default)s")


  ## Parser creation
  parser = argparse.ArgumentParser(description = "Rescue artifacts discarded by \
    the SQANTI3 filter, i.e. find closest match for the artifacts in the reference \
    transcriptome and add them to the transcriptome. \
    \nChoose between the filter applied: using rules or the Machine-Learning approach.")
  # Subparsers
  subparsers = parser.add_subparsers(dest = 'subcommand')
  ## Rules rescue arguments
  rules = subparsers.add_parser("rules", 
                                parents = [common],
                                description = "Rescue for rules-filtered transcriptomes.",
                                formatter_class = argparse.RawTextHelpFormatter)
  # rules options
  rf = rules.add_argument_group("Rules specific options")
  rf.add_argument("-j", "--json_filter",
                  default = default_json,
                  help = "Full path to the JSON file including the rules used when running the SQANTI3 rules filter. \
                    \nDefault: %(default)s")

  ## ML rescue arguments
  machine_learning = subparsers.add_parser("ml", 
                             parents = [common],
                             description = "Rescue for ML-filtered transcriptomes.",
                             formatter_class = argparse.RawTextHelpFormatter)
  # ML options
  ml = machine_learning.add_argument_group("Machine Learning specific options")
  ml.add_argument("-r", "--random_forest",
                  help = "Full path to the randomforest.RData object obtained when running the SQANTI3 ML filter.")
  ml.add_argument("-j", "--threshold", 
                  type = float, default = 0.7, 
                  help = "Machine learning probability threshold to filter elegible rescue targets (mapping hits). \
                    \nDefault: %(default)s")


# parse arguments
  return parser