#!/usr/bin/env python3
__author__  = "francisco.pardo.palacios@gmail.com"
__version__ = '1.0'   # Python 3.7 syntax!

"""
SQANTI3 Comparison Tool. It will compare the results of different SQANTI3 runs.

It is a wrapper to run two modes of comparing results:
* Normal or sample or comparison: It will compare the two or more transcriptomes that were analyzed through SQ3.
This can be specially useful if there is of interest the differences between transcriptomes 
Isoform IDs can be different and are not related within classification files

* Filter comparison: This would be used to compare different ML results depending on the set of TP or FP defined

"""

import os, sys, argparse, subprocess
import distutils.spawn
from csv import DictReader, DictWriter
from Bio import SeqIO
from cupcake.io.BioReaders import GMAPSAMReader
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
RSCRIPTPATH = distutils.spawn.find_executable('Rscript')
RSCRIPT_SAMPLE_REPORT = 'SQANTI3_comparison_report.R'
RSCRIPT_SAMPLE_MARKDOWN = 'SQANTI3_comparison_report.Rmd'

if os.system(RSCRIPTPATH + " --version")!=0:
    print("Rscript executable not found! Abort!", file=sys.stderr)
    sys.exit(-1)

def main():
    parser = argparse.ArgumentParser(description="Comparison of transcriptomes based on SQANTI3 QC results.\
\nChoose between a 'samples' comparison or a comparative analysis about the effects of different filtering strategies.")
### Common arguments for both modes
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('sqanti_results', help='\t\tSQANTI3 output folder where classification and junction files to compare are stored.')
## It could be interesting to add options so the comparison can return let's say the "common" isoforms or the "specific ones
#    common.add_argument('--isoforms', help='\t\tfasta/fastq isoform file to be filtered')
#    common.add_argument('--gtf', help='\t\tGTF file to be filtered')
#    common.add_argument('--sam', help='\t\tSAM alignment of the input fasta/fastq')
#    common.add_argument('--faa', help='\t\tORF prediction faa file to be filtered by SQANTI3')
    common.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
    common.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
    subparsers = parser.add_subparsers(dest='subcommand')
### Rules filter arguments
    samples = subparsers.add_parser('samples', parents=[common], description="Samples comparison selected")
    samples.add_argument('-e', "--exclude_mono_exonic", action="store_true", default=False, help='\t\tMono-exonic transcripts will be excluded from any comparison.')
### ML filter arguments
    filters = subparsers.add_parser('filters', parents=[common],  description='Filters comparison selected')
    filters.add_argument('-t','--percent_training', type=float, default=0.8, \
    help="Proportion of the data that goes to training (parameter p of the function createDataPartition). \
    \nDefault: 0.8")

    args = parser.parse_args()

### Checking presence of files. Common arguments
    args.sqanti_results = os.path.abspath(args.sqanti_results)
    if not os.path.exists(args.sqanti_results):
        print("ERROR: {0} doesn't exist. Abort!".format(args.sqanti_results), file=sys.stderr)
        sys.exit(-1)

#    if args.isoforms is not None and not os.path.exists(args.isoforms):
#        print("ERROR: {0} doesn't exist. Abort!".format(args.isoform), file=sys.stderr)
#        sys.exit(-1)

#    if args.gtf is not None and not os.path.exists(args.gtf):
#        print("ERROR: {0} doesn't exist. Abort!".format(args.gtf_file), file=sys.stderr)
#        sys.exit(-1)

#    if args.sam is not None and not os.path.exists(args.sam):
#        print("ERROR: {0} doesn't exist. Abort!".format(args.sam), file=sys.stderr)
#        sys.exit(-1)

#    if args.faa is not None and not os.path.exists(args.faa):
#        print("ERROR: {0} doesn't exist. Abort!".format(args.faa), file=sys.stderr)
#        sys.exit(-1)
#    if args.intrapriming < 0.25 or args.intrapriming > 1.:
#        print("ERROR: --intrapriming must be between 0.25-1, instead given {0}! Abort!".format(args.intrapriming), file=sys.stderr)
### Define output dir and output name in case it was not defined
    if args.dir is None:
        args.dir=os.path.dirname(args.sqanti_class)
        print("Output directory not defined. All the outputs will be stored at {0} directory".format(args.dir), file=sys.stderr)
    else:
        if not os.path.exists(args.dir):
            os.makedirs(args.dir)
    if args.output is None:
        args.output=os.path.basename(args.sqanti_results)
        print("Output name not defined. All the outputs will have the prefix {0}".format(args.output), file=sys.stderr)

### Check arguments samples mode --> If ok run comparison scripts
    if args.subcommand == 'samples':
#        if args.TP is not None and not os.path.exists(args.TP):
#            print("ERROR: {0} doesn't exist. Abort!".format(args.TP), file=sys.stderr)
#            sys.exit(-1)
        print("\nRunning SQANTI3 comparison with Samples mode...", file=sys.stdout)
        run_samples_comparison(args)
#       filter_files(args, ids)


### Check arguments --> If ok run filters comparison

    if args.subcommand == 'filters':
#        if args.runAlength < 4 or args.runAlength > 20:
#            print("ERROR: --runAlength must be between 4-20, instead given {0}! Abort!".format(args.runAlength), file=sys.stderr)
#            sys.exit(-1)

        run_filters_comparison(args)
#        filter_files(args, ids)

if __name__ == "__main__":
    main()



