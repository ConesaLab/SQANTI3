import argparse
from src.config import __version__, default_json

def rescue_argparse():

# Arguments and help
  ## Common arguments
  ## Parser creation
  parser = argparse.ArgumentParser(description = "Rescue artifacts discarded by \
    the SQANTI3 filter, i.e. find closest match for the artifacts in the reference \
    transcriptome and add them to the transcriptome. \
    \nChoose between the filter applied: using rules or the Machine-Learning approach.")

  cr = parser.add_argument_group("Required arguments")
  cr.add_argument("--filter_class",
                  required=True,
                  help = "SQANTI filter (ML or rules) output classification file.")
  cr.add_argument("-rg","--refGTF",
                  required=True,
                  help = "Full path to reference transcriptome GTF used when running SQANTI3 QC.")
  cr.add_argument("-rf","--refFasta",
                  required=True,
                  help = "Full path to reference genome FASTA used when running SQANTI3 QC.")
  
  # Specific input options
  ci = parser.add_argument_group("Input options")
  ci.add_argument("--corrected_isoforms_fasta", 
                      help = "FASTA file output by SQANTI3 QC (*_corrected.fasta), i.e. the full long read transcriptome.")
  ci.add_argument("--filtered_isoforms_gtf", 
                      help = "GTF file output by SQANTI3 filter (*.filtered.gtf).")
  ci.add_argument("-k", "--refClassif", 
                      help = "Full path to the classification file obtained when running SQANTI3 QC on the reference transcriptome.\
                        \nMandatory when running the rescue on full mode")
  ci.add_argument("--counts",
                      help = 'Isoforms abundancy values: "Isoform" \t "Count". Column names may differ')
  
  # Customization options
  cc = parser.add_argument_group("Customization options")
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
  cc.add_argument("--skip_requant",
                  action="store_true",
                  help = "Skip the requantification of the rescued isoforms, which redistributes counts \
                    \nfrom discarded artifacts to their replacement transcripts and requires --counts.")
  cc.add_argument("-s","--strategy", 
                  choices = ["rules", "ml"],
                  default = "rules", 
                  help = "Filter strategy used.\
                    \nDefault: %(default)s")
  
  # rules options
  rf = parser.add_argument_group("Rules specific options")
  rf.add_argument("-j", "--json_filter",
                  default = default_json,
                  help = "Full path to the JSON file including the rules used when running the SQANTI3 rules filter. \
                    \nDefault: %(default)s")
  rf.add_argument("--counts_design",
                  help = "TSV/CSV file with the experimental groups given to SQANTI3 QC (--counts_design). \
                    \nNeeded when the rules use prevalence_<group> columns, so that the evidence check \
                    \ncomputes them on the counts that rescued reference transcripts receive from their artifacts.")
  rf.add_argument("--min_expression", type = float, default = 0,
                  help = "Minimum count for a transcript to be expressed in a sample, as given to SQANTI3 QC (--min_expression): \
                    \n0 means any count above 0, any other value a count greater than or equal to it. \
                    \nUsed by the evidence check to compute prevalence and prevalence_<group> on the counts \
                    \nthat rescued reference transcripts receive from their artifacts. \
                    \nDefault: %(default)s")
  rf.add_argument("--skip_evidence_check", action="store_true",
                  help = "Do not apply the requisites on long-read counts (FL, prevalence, prevalence_<group>) \
                    \nto the rescued reference transcripts (behaviour of previous versions).")
  rf.add_argument("--map_automatic_fsm", action="store_true",
                  help = "Map the FSM artifacts of automatic rescue (needs --corrected_isoforms_fasta), so that \
                    \nthey are reassigned to a long-read isoform if their reference transcript fails the evidence check. \
                    \nBy default their counts go to the gene residual.")

  # ML options
  ml = parser.add_argument_group("Machine Learning specific options")
  ml.add_argument("-r", "--random_forest",
                  help = "Full path to the randomforest.RData object obtained when running the SQANTI3 ML filter.")
  ml.add_argument("-t", "--threshold", 
                  type = float, default = 0.7, 
                  help = "Machine learning probability threshold to filter elegible rescue targets (mapping hits). \
                    \nDefault: %(default)s")
  # Output options
  co = parser.add_argument_group("Output options")
  co.add_argument("-o","--output", 
                      default = "isoform",
                      help = "Prefix for output files.", 
                      required = False)
  co.add_argument("-d","--dir", 
                      default = "sqanti3_output",
                      help = "Directory for output files. Default: Directory where the script was run.", 
                      required = False)
  # Performance options
  cp = parser.add_argument_group("Extra options")
  cp.add_argument("-c", "--cpus",
                      type = int, default = 4, 
                      help = "Number of CPUs to use. Default: 4")
  cp.add_argument("-v", "--version", 
                      help="Display program version number.", 
                      action='version', 
                      version=f"SQANTI3 v{__version__}")
  cp.add_argument("-l","--log_level", default="INFO",choices=["ERROR","WARNING","INFO","DEBUG"],
                      help="Set the logging level %(default)s")

# parse arguments
  return parser