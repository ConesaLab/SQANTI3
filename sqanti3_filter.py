#!/usr/bin/env python3
__author__  = "francisco.pardo.palacios@gmail.com"

# New SQANTI3 filter. It will serve as a wrapper for "rules" filter and "Machine-Learning" filter.

# RULES FILTER --> Now it can work with a JSON filter
# By default, it will only keep Iso-Seq isoforms if:
# The isoform is FSM and does not have intrapriming.
# The isoform is ISM, NIC or NNC, does not have intrapriming nor RT-switching, and all junctions are either all canonical or short-read-supported
# The isoform is antisense, intergenic, genic, does not have intrapriming nor RT-switching, and all junctions are either all canonical or short-read-supported

# If the user wants to define new rules, it can define them in a JSON file following the same format used in the filter_default.json

# ML FILTER
# It will take as input the classification file obtained from SQANTI3 QC and apply a Random Forest algorithm to distinguish between "true" isoforms and artifacts.

# Regardless of the strategy chosen, sqanti_filter.py can return a filtered FASTA, filtered GTF and an updated classification file that can be used to compare the effects of filtering out
# bad quality transcripts.


import os

from src.filter_argparse import filter_argparse
from src.module_logging import filter_logger, update_logger
from src.config import __version__
from src.logging_config import art_logger,filter_art, get_logger_info
from src.filter_steps import filter_files, run_ML, run_rules
from src.argparse_utils import filter_args_validation
from src.write_parameters import write_filter_parameters

def main():
    global filter_logger
    # Create the log directory if it does not exist
    art_logger.info(filter_art())
    args = filter_argparse().parse_args()
    update_logger(filter_logger,args.dir,"filter",args.log_level)
    # Create the log directory if it does not exist
    os.makedirs(f"{args.dir}/logs", exist_ok=True)
    filter_args_validation(args)
    if args.output is None:
        args.output = args.sqanti_class[args.sqanti_class.rfind("/")+1:args.sqanti_class.rfind("_classification.txt")]
        filter_logger.warning(f"Output name not defined. All the outputs will have the prefix {args.output}")

    write_filter_parameters(args)

### Checking presence of files for ML. Check arguments --> If ok run ML
    filter_logger.info("Running SQANTI3 filtering...")
    if args.subcommand == 'ml':
        ids, inclusion_file = run_ML(args)
### Checking presence of files for Rules. Check arguments --> If ok run Rules
    if args.subcommand == 'rules':
        ids, inclusion_file = run_rules(args)
    
    filter_files(args.filter_isoforms, args.filter_gtf, args.filter_sam, args.filter_faa,
                 args.isoAnnotGFF3, args.dir,args.output, ids, inclusion_file)

if __name__ == "__main__":
    main()
