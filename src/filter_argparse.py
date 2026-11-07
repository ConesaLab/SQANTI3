import argparse
from src.config import __version__,__author__, default_json

# Subcommands of sqanti3_filter.py; the first one is the default.
SUBCOMMANDS = ("auto", "rules", "ml")

def filter_argparse():

    ### Common arguments for both modes
    common = argparse.ArgumentParser(add_help=False)
    # Required arguments
    cr = common.add_argument_group("Required arguments")
    cr.add_argument('--sqanti_class', required=True, help='SQANTI3 QC classification file.')
    # Input options
    ci = common.add_argument_group("Input options")
    ci.add_argument('--isoAnnotGFF3', help='isoAnnotLite GFF3 file to be filtered')
    ci.add_argument('--filter_isoforms', help='fasta/fastq isoform file to be filtered')
    ci.add_argument('--filter_gtf', help='GTF file to be filtered')
    ci.add_argument('--filter_sam', help='SAM alignment of the input fasta/fastq')
    ci.add_argument('--filter_faa', help='ORF prediction faa file to be filtered by SQANTI3')
    # Output options
    co = common.add_argument_group("Output options")
    co.add_argument('-o','--output', help='Prefix for output files.')
    co.add_argument('-d','--dir', default='./sqanti3_results', help='Directory for output files. Default: %(default)s')
    co.add_argument("--skip_report", action="store_true", help='Skip creation of a report about the filtering')
    # Filtering options
    cf = common.add_argument_group("Filtering options")
    cf.add_argument("-e","--filter_mono_exonic", action="store_true", help='All mono-exonic transcripts are automatically filtered')
    # Performance options
    cp = common.add_argument_group("Extra options")
    cp.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    cp.add_argument("-c", "--cpus", type=int, default=4, help="Number of CPUs to use. Default: %(default)s")
    cp.add_argument("-l","--log_level", default="INFO",choices=["ERROR","WARNING","INFO","DEBUG"],
                        help="Set the logging level.\nDefault: %(default)s")

# Parser creations
    parser = argparse.ArgumentParser(description="Filtering of Isoforms based on SQANTI3 attributes.\
\nChoose between the automatic filter (default), a rules filter or a Machine-Learning based filter.") 
    subparsers = parser.add_subparsers(dest='subcommand')

### Automatic filter arguments (default subcommand)
    auto = subparsers.add_parser('auto',
                                 parents=[common],
                                 description="Automatic filter selected (default). Multi-exon NIC and NNC are filtered with graded pubLR"
                                     "\n(kept with pubLR_prob >= 0.7; between 0.5 and 0.7, only if they also pass the rules; artifacts below 0.5)"
                                     "\nor, when the junctions have short-read coverage, with the random forest of the ML filter cut where it keeps"
                                     "\n90% of the NIC (and NNC) that are recurrent in public data. Intra-primed NIC/NNC are artifacts."
                                     "\nThe other isoforms get the rules filter. Needs a QC run with --pantranscriptome;"
                                     "\nwithout it, the rules filter is applied to every isoform.",
                                 formatter_class=argparse.RawTextHelpFormatter)
    af = auto.add_argument_group("Automatic filter specific options")
    af.add_argument('-j', "--json_filter", default=default_json,
                    help="JSON file with the rules applied outside multi-exon NIC/NNC and to the NIC/NNC with pubLR_prob between 0.5 and 0.7."
                        "\nDefault: %(default)s")
    af.add_argument('-t','--percent_training', type=float, default=0.8,
                    help="Proportion of the data that goes to training the random forest.\nDefault: %(default)s")
    af.add_argument('-p', '--TP',
                    help="Path to file containing the list of the TP transcripts for the random forest, one ID by line, no header (optional).")
    af.add_argument('-n', '--TN',
                    help="Path to file containing the list of the TN transcripts for the random forest, one ID by line, no header (optional).")
    af.add_argument('-z', '--max_class_size', type=int , default=3000,
                    help="Maximum number of isoforms to include in True Positive and True Negative sets.\nDefault: %(default)s")
    af.add_argument('-i',"--intrapriming", type=float, default=60,
                    help='Adenine percentage at genomic 3\' end from which a NIC/NNC is intra-primed. Default: %(default)s')

### Rules filter arguments
    rules = subparsers.add_parser('rules', 
                                  parents=[common],
                                  description="Rules filter selected",
                                  formatter_class=argparse.RawTextHelpFormatter)
    # filtering options
    rf = rules.add_argument_group("Rules specific options")
    rf.add_argument('-j', "--json_filter", default=default_json, 
                    help="JSON file where filtering rules are expressed. Rules must be set taking into account that attributes described in the filter will be present in those isoforms that should be kept."
                        "\nDefault: %(default)s")

### ML filter arguments
    machine_learning = subparsers.add_parser('ml', 
                               parents=[common],
                               description='ML filter selected',
                               formatter_class=argparse.RawTextHelpFormatter)
    # ML options
    ml = machine_learning.add_argument_group("Machine Learning specific options")
    ml.add_argument('-t','--percent_training', type=float, default=0.8, \
                    help="Proportion of the data that goes to training (parameter p of the function createDataPartition). \
                    \nDefault: %(default)s")
    ml.add_argument('-p', '--TP', \
                    help="Path to file containing the list of the TP transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.")
    ml.add_argument('-n', '--TN', \
                    help="Path to file containing the list of the TN transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.")
    ml.add_argument('-j', '--threshold', type=float, default=0.7, \
                    help="Machine Learning probability threshold to classify transcripts as positive isoforms. \
        \nDefault: %(default)s")
    ml.add_argument('-f', '--force_fsm_in', action="store_true", \
                    help="Forces retaining FMS transcripts regardless of ML filter result (FSM are threfore automatically classified as isoforms).")
    ml.add_argument('--intermediate_files', action="store_true", \
                    help="Outputs ML filter intermediate files.")
    ml.add_argument('-r','--remove_columns', \
                    help="Path to single-column file (no header) containing the names of the columns in SQ3's classification.txt file that are to be excluded during random forest training (optional).")
    ml.add_argument('-z', '--max_class_size', type=int , default=3000, \
                    help="Maximum number of isoforms to include in True Positive and True Negative sets. TP and TN sets will be downsized to this value if they are larger.\nDefault: %(default)s")
    ml.add_argument('-i',"--intrapriming", type=float, default=60, help='Adenine percentage at genomic 3\' end to flag an isoform as intra-priming. Default: %(default)s')

    return parser


def default_subcommand(argv: list[str]) -> list[str]:
    """Put the default subcommand first when none is given.

    Args:
        argv: Command-line arguments, without the program name.

    Returns:
        argv, preceded by "auto" unless it starts with a subcommand or
        asks for the help or the version.
    """
    result = list(argv)
    explicit = {*SUBCOMMANDS, "-h", "--help", "-v", "--version"}
    if not argv or argv[0] not in explicit:
        result = ["auto", *argv]

    return result
