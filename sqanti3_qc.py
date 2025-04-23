#!/usr/bin/env python3
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila
# Modified by Liz (etseng@pacb.com) as SQANTI2/3 versions
# Modified by Fran (francisco.pardo.palacios@gmail.com) currently as SQANTI3 version (05/15/2020)
# Modified by Pablo (pabloatienzalo@gmail.com)

import os
import shutil
# Import SQANTI3 modules
from src.qc_argparse import qc_argparse
from src.qc_pipeline import run
from src.parallel import split_input_run, combine_split_runs, get_split_dir
from src.config import __version__
from src.argparse_utils import qc_args_validation
from src.logging_config import qc_art,art_logger
from src.module_logging import qc_logger
from src.write_parameters import write_qc_parameters
def main():

    art_logger.info(qc_art())
    args = qc_argparse().parse_args()
    # Check if the output directory exists, if not create it
    os.makedirs(f"{args.dir}/logs", exist_ok=True)
    args = qc_args_validation(args)
    if not os.path.exists(os.path.join(args.dir,'logs')):
        os.makedirs(os.path.join(args.dir,'logs'))
    write_qc_parameters(args)

    # Running functionality based on the chunks
    qc_logger.info(f"Initialising QC pipeline.")
    
    if args.chunks == 1:
        run(args)

    else:
        SPLIT_ROOT_DIR = get_split_dir(args.dir,args.output)
        split_dirs = split_input_run(args, SPLIT_ROOT_DIR)
        combine_split_runs(args,split_dirs)
        shutil.rmtree(SPLIT_ROOT_DIR)

    if args.isoAnnotLite:
        from src.helpers import get_corr_filenames, get_class_junc_filenames
        from src.qc_pipeline import run_isoAnnotLite 
        corrGTF, _, _, _ , _ = get_corr_filenames(args.dir,args.output)
        outputClassPath, outputJuncPath = get_class_junc_filenames(args.dir,args.output)
        run_isoAnnotLite(corrGTF, outputClassPath, outputJuncPath, args.output, args.gff3)

if __name__ == "__main__":
    main()