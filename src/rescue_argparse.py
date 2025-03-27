import argparse
from .config import __version__

def rescue_argparse():
  
# Arguments and help
    parser = argparse.ArgumentParser(description = "Rescue artifacts discarded by \
    the SQANTI3 filter, i.e. find closest match for the artifacts in the reference \
    transcriptome and add them to the transcriptome.")

    ## Common arguments
    common = argparse.ArgumentParser(add_help = False)
    common.add_argument("--filter_class", 
                        help = "SQANTI filter (ML or rules) output classification file.")
    common.add_argument("--rescue_isoforms", 
                        help = "FASTA file output by SQANTI3 QC (*_corrected.fasta), i.e. the full long read transcriptome.")
    # TODO: check why here is the filtered, but not above??
    common.add_argument("--rescue_gtf", 
                        help = "GTF file output by SQANTI3 filter (*.filtered.gtf).")
    common.add_argument("--refGTF", 
                        help = "Full path to reference transcriptome GTF used when running SQANTI3 QC.")
    common.add_argument("--refFasta", 
                        help = "Full path to reference genome FASTA used when running SQANTI3 QC.")
    common.add_argument("-k", "--refClassif", 
                        help = "Full path to the classification file obtained when running SQANTI3 QC on the reference transcriptome.")
    common.add_argument("-e","--rescue_mono_exonic", 
                        choices = ['all', 'fsm', 'none'], 
                        default = "all", 
                        help='Whether or not to include mono-exonic artifacts in the rescue. Options include: none, fsm and all (default).')
    common.add_argument("--mode", 
                        choices = ["automatic", "full"],
                        default = "automatic", 
                        help = "If 'automatic' (default), only automatic rescue of FSM artifacts will be performed. If 'full', rescue will include mapping of ISM, NNC and NIC artifacts to find potential replacement isoforms.")
    common.add_argument("-o","--output", 
                        help = "Prefix for output files.", 
                        required = False)
    common.add_argument("-d","--dir", 
                        help = "Directory for output files. Default: Directory where the script was run.", 
                        required = False)
    common.add_argument("-c", "--cpus",
                        type = int, default = 4, 
                        help = "Number of CPUs to use. Default: 4")
    common.add_argument("-v", "--version", 
                        help="Display program version number.", 
                        action='version', 
                        version=f"SQANTI3 f{__version__}")
    common.add_argument("-l","--log_level", default="INFO",choices=["ERROR","WARNING","INFO","DEBUG"],
                        help="Set the logging level %(default)s")
    # Subparsers
    subparsers = parser.add_subparsers(dest = 'subcommand')
    ## Rules rescue arguments
    rules = subparsers.add_parser("rules", parents = [common], \
    description = "Rescue for rules-filtered transcriptomes.")

    rules.add_argument("-j", "--json_filter", \
    help = "Full path to the JSON file including the rules used when running the SQANTI3 rules filter.")

    ## ML rescue arguments
    ml = subparsers.add_parser("ml", parents = [common], \
    description = "Rescue for ML-filtered transcriptomes.")

    ml.add_argument("-r", "--randomforest",
                    help = "Full path to the randomforest.RData object obtained when running the SQANTI3 ML filter.")
    ml.add_argument("-j", "--threshold", 
                    type = float, default = 0.7, 
                    help = "Default: 0.7. Machine learning probability threshold to filter elegible rescue targets (mapping hits).")


  # parse arguments
    return parser