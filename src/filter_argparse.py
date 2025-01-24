import argparse
from .config import __version__,__author__
from .argparse_utils import *
from .commands import utilitiesPath

def filter_argparse():
    default_json = utilitiesPath + "/filter/filter_default.json"

### Common arguments for both modes
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('sqanti_class', help='SQANTI3 QC classification file.')
    common.add_argument('--isoAnnotGFF3', help='isoAnnotLite GFF3 file to be filtered')
    common.add_argument('--isoforms', help='fasta/fastq isoform file to be filtered')
    common.add_argument('--gtf', help='GTF file to be filtered')
    common.add_argument('--sam', help='SAM alignment of the input fasta/fastq')
    common.add_argument('--faa', help='ORF prediction faa file to be filtered by SQANTI3')
    common.add_argument('-o','--output', help='Prefix for output files.', required=False)
    common.add_argument('-d','--dir', help='Directory for output files. Default: Directory where the script was run.', required=False)
    common.add_argument("-e","--filter_mono_exonic", action="store_true", default=False, help='When TRUE, all mono-exonic transcripts are automatically filtered (default: False)')
    common.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    common.add_argument("-c", "--cpus", type=int, default=4, help="Number of CPUs to use. Default: 4")
    common.add_argument("--skip_report", action="store_true", default=False, help='Skip creation of a report about the filtering')

# Parser creations
    parser = argparse.ArgumentParser(description="Filtering of Isoforms based on SQANTI3 attributes.\
\nChoose between a rules filter or a Machine-Learning based filter.") 
    subparsers = parser.add_subparsers(dest='subcommand')

### Rules filter arguments
    rules = subparsers.add_parser('rules', 
                                  parents=[common],
                                  description="Rules filter selected")
    rules.add_argument('-j', "--json_filter", default=default_json, help='JSON file where filtering rules are expressed. Rules must be set taking into account that attributes described in the filter will be present in those isoforms that should be kept. Default: utilities/filter/filter_default.json')

### ML filter arguments
    ml = subparsers.add_parser('ml', 
                               parents=[common],
                               description='ML filter selected')
    ml.add_argument('-t','--percent_training', type=float, default=0.8, \
    help="Proportion of the data that goes to training (parameter p of the function createDataPartition). \
    \nDefault: 0.8")
    ml.add_argument('-p', '--TP', \
    help="Path to file containing the list of the TP transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.")
    ml.add_argument('-n', '--TN', \
    help="Path to file containing the list of the TN transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.")
    ml.add_argument('-j', '--threshold', type=float, default=0.7, \
    help="Machine Learning probability threshold to classify transcripts as positive isoforms. Default: 0.7.")
    ml.add_argument('-f', '--force_fsm_in', default=False, \
    help="When TRUE, forces retaining FMS transcripts regardless of ML filter result (FSM are threfore automatically classified as isoforms). Default: FALSE.")
    ml.add_argument('--intermediate_files', default=False, \
    help="When TRUE, outputs ML filter intermediate files. Default: FALSE.")
    ml.add_argument('-r','--remove_columns', \
    help="Path to single-column file (no header) containing the names of the columns in SQ3's classification.txt file that are to be excluded during random forest training (optional).")
    ml.add_argument('-z', '--max_class_size', type=int , default=3000, \
    help="Maximum number of isoforms to include in True Positive and True Negative sets (default: 3000). TP and TN sets will be downsized to this value if they are larger.")
    ml.add_argument('-i',"--intrapriming", type=float, default=60, help='Adenine percentage at genomic 3\' end to flag an isoform as intra-priming (default: 60 )')

    # for subparser in [rules, ml]:
    #     subparser.add_argument('sqanti_class', help='SQANTI3 QC classification file.')

    return parser