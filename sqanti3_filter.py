

#!/usr/bin/env python3
__author__  = "francisco.pardo.palacios@gmail.com"

"""
New SQANTI3 filter. It will serve as a wrapper for "rules" filter and "Machine-Learning" filter.

RULES FILTER --> Now it can work with a JSON filter
By default, it will only keep Iso-Seq isoforms if:
The isoform is FSM and does not have intrapriming.
The isoform is ISM, NIC or NNC, does not have intrapriming nor RT-switching, and all junctions are either all canonical or short-read-supported
The isoform is antisense, intergenic, genic, does not have intrapriming nor RT-switching, and all junctions are either all canonical or short-read-supported

If the user wants to define new rules, it can define them in a JSON file following the same format used in the filter_default.json

ML FILTER
It will take as input the classification file obtained from SQANTI3 QC and apply a Random Forest algorithm to distinguish between "true" isoforms and artifacts.

Regardless of the strategy chosen, sqanti_filter.py can return a filtered FASTA, filtered GTF and an updated classification file that can be used to compare the effects of filtering out
bad quality transcripts.

"""

import os, sys

from src.filter_argparse import filter_argparse
from src.module_logging import filter_logger
from src.config import __version__
from src.logging_config import art_logger,filter_art
from src.filter_steps import filter_files
from src.commands import utilitiesPath

def main():
    art_logger.info(filter_art())
    args = filter_argparse().parse_args()

    ### Checking presence of files. Common arguments
    args.sqanti_class = os.path.abspath(args.sqanti_class)
    if not os.path.isfile(args.sqanti_class):
        filter_logger.error(f"{args.sqanti_class} doesn't exist. Abort!")
        sys.exit(-1)

    if args.filter_isoforms is not None and not os.path.exists(args.filter_isoforms):
        filter_logger.error(f"{args.filter_isoforms} doesn't exist. Abort!")
        sys.exit(-1)

    if args.filter_gtf is not None and not os.path.exists(args.filter_gtf):
        filter_logger.error(f"{args.filter_gtf} doesn't exist. Abort!")
        sys.exit(-1)

    if args.filter_sam is not None and not os.path.exists(args.filter_sam):
        filter_logger.error(f"{args.filter_sam} doesn't exist. Abort!")
        sys.exit(-1)

    if args.filter_faa is not None and not os.path.exists(args.filter_faa):
        filter_logger.error(f"{args.filter_faa} doesn't exist. Abort!")
        sys.exit(-1)

    ### Define output dir and output name in case it was not defined
    if args.dir is None:
        args.dir = os.path.dirname(args.sqanti_class)
        filter_logger.warning(f"Output directory not defined. All the outputs will be stored at {args.dir} directory")
    else:
        if not os.path.exists(args.dir):
            os.makedirs(args.dir)

    if args.output is None:
        args.output = args.sqanti_class[args.sqanti_class.rfind("/")+1:args.sqanti_class.rfind("_classification.txt")]
        filter_logger.warning(f"Output name not defined. All the outputs will have the prefix {args.output}")

### Print out parameters so can be reproduced the same SQ run
    args.doc = os.path.join(os.path.abspath(args.dir), args.output+"_params.txt")
    filter_logger.info(f"Write arguments to {args.doc}...")
    with open(args.doc, 'w') as f:
      f.write("Version\t" + __version__ + "\n")
      f.write("Mode\t" + args.subcommand + "\n")
      f.write("ClassificationFile\t" + str(args.sqanti_class) + "\n")
      f.write("Isoforms\t" + (str(args.filter_isoforms) if args.filter_isoforms is not None else "NA")+ "\n")
      f.write("GTF\t" + (str(args.filter_gtf) if args.filter_gtf is not None else "NA") + "\n")
      f.write("SAM\t" + (str(args.filter_sam) if args.filter_sam is not None else "NA") + "\n")
      f.write("FAA\t" + (str(args.filter_faa) if args.filter_faa is not None else "NA") + "\n")
      f.write("isoAnnotGFF3\t" + (str(args.isoAnnotGFF3) if args.isoAnnotGFF3 is not None else "NA") + "\n")
      f.write("OutputPrefix\t" + str(args.output) + "\n")
      f.write("OutputDirectory\t" + os.path.abspath(args.dir) + "\n")
      f.write("FilterMonoexonic\t" + str(args.filter_mono_exonic) + "\n")
      f.write("SkipReport\t" + str(args.skip_report) + "\n")
      if args.subcommand == 'rules':
          f.write("JSON\t" + str(args.json_filter) + "\n")
      if args.subcommand == 'ml':
          f.write("PercentTraining\t" + str(args.percent_training) + "\n")
          f.write("TP\t" + (str(args.TP) if args.TP is not None else "NA") + "\n")
          f.write("TN\t" + (str(args.TN) if args.TN is not None else "NA") + "\n")
          f.write("Threshold\t" + str(args.threshold) + "\n")
          f.write("ForceFSM\t" + str(args.force_fsm_in) + "\n")
          f.write("KeepIntermediate\t" + str(args.intermediate_files) + "\n")
          f.write("ColumnsRemoved\t" + (str(args.remove_columns) if args.remove_columns is not None else "NA") + "\n")
          f.write("MaxClassSize\t" + str(args.max_class_size) + "\n")
          f.write("Intrapriming\t" + str(args.intrapriming) + "\n")

### Checking presence of files for ML. Check arguments --> If ok run ML

    filter_logger.info("Running SQANTI3 filtering...")

    if args.subcommand == 'ml':
        if args.TP is not None and not os.path.exists(args.TP):
            filter_logger.error(f"{args.TP} doesn't exist. Abort!")
            sys.exit(-1)
        if args.TN is not None and not os.path.exists(args.TN):
            filter_logger.error(f"{args.TN} doesn't exist. Abort!")
            sys.exit(-1)
        if args.remove_columns is not None and not os.path.exists(args.remove_columns):
            filter_logger.error(f"{args.remove_columns} doesn't exist. Abort!")
            sys.exit(-1)
        if args.percent_training < 0 or args.percent_training > 1.:
            filter_logger.error(f"--percent_training must be between 0-1, instead given {args.percent_training}! Abort!")
            sys.exit(-1)
        if args.intrapriming < 25 or args.intrapriming > 100:
            filter_logger.error(f"--intrapriming must be between 25-100, instead given {args.intrapriming}! Remember to use the percentage value. Abort!")
            sys.exit(-1)


        ids, inclusion_file = run_ML(args)
        filter_files(args, ids, inclusion_file)

### Checking presence of files for Rules. Check arguments --> If ok run Rules

    if args.subcommand == 'rules':
        if args.json_filter is not None and not os.path.exists(args.json_filter):
            filter_logger.error(f"{args.json_filter} doesn't exist. Abort!")
            sys.exit(-1)

        ids, inclusion_file = run_rules(args)
        filter_files(args, ids, inclusion_file)

if __name__ == "__main__":
    main()
