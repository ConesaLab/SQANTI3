#!/usr/bin/env python3
__author__  = "francisco.pardo.palacios@gmail.com"
__version__ = '5.0'   # Python 3.7 syntax!

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

import os, sys, argparse, subprocess
import distutils.spawn
from csv import DictReader, DictWriter
from Bio import SeqIO
from cupcake.io.BioReaders import GMAPSAMReader
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
RSCRIPTPATH = distutils.spawn.find_executable('Rscript')
RSCRIPT_REPORT = 'report_filter/SQANTI3_filter_report.R'
RSCRIPT_ML = 'filter/SQANTI3_MLfilter.R'
RSCRIPT_RULES = 'filter/SQANTI3_rules_filter.R'
default_json = utilitiesPath + "/filter/filter_default.json"

if os.system(RSCRIPTPATH + " --version")!=0:
    print("Rscript executable not found! Abort!", file=sys.stderr)
    sys.exit(-1)

def filter_files(args, ids_to_keep, inclusion_f):
    prefix = args.dir + "/" + args.output
    # filter FASTA/FASTQ file
    if args.isoforms is not None:
        fafq_type = 'fasta'
        with open(args.isoforms) as h:
            if h.readline().startswith('@'): fafq_type = 'fastq'
        fout=open(prefix + '.filtered.' + fafq_type, 'w')
        for r in SeqIO.parse(open(args.isoforms), fafq_type):
            if r.id in ids_to_keep:
                SeqIO.write(r, fout, fafq_type)
        fout.close()
        print("Output written to: {0}".format(fout.name), file=sys.stdout)

    # filter GTF
    if args.gtf is not None:       
        outputGTF = prefix + '.filtered.gtf'
        with open(outputGTF, 'w') as f:
            for r in collapseGFFReader(args.gtf):
                if r.seqid in ids_to_keep:
                    write_collapseGFF_format(f, r)
            print("Output written to: {0}".format(f.name), file=sys.stdout)

    # filter SAM
    if args.sam is not None:
        outputSam = prefix + '.filtered.sam'
        with open(outputSam, 'w') as f:
            reader = GMAPSAMReader(args.sam, True)
            f.write(reader.header)
            for r in reader:
                if r.qID in ids_to_keep:
                    f.write(r.record_line + '\n')
            print("Output written to: {0}".format(f.name), file=sys.stdout)

    # filter FAA 
    if args.faa is not None:
        outputFAA = prefix + '.filtered.faa'
        with open(outputFAA, 'w') as f:
            for r in SeqIO.parse(open(args.faa), 'fasta'):
                if r.id in ids_to_keep:
                    f.write(">{0}\n{1}\n".format(r.description, r.seq))
        print("Output written to: {0}".format(f.name), file=sys.stdout)

    # filter isoAnnot GFF3
    if args.isoAnnotGFF3 is not None:
        outputGFF3 = prefix + '.filtered.gff3'
        awk_cmd = """awk 'FNR==NR {{ a[$1]; next }} ($1 in a)' {l} {g} > {o}""".format(l=inclusion_f, g=args.isoAnnotGFF3, o=outputGFF3)
        subprocess.call(awk_cmd, shell=True)
        print("Output written to: {0}".format(outputGFF3), file=sys.stdout)

def run_ML(args):
    cmd = RSCRIPTPATH + " {u}/{s} -c {c} -o {o} -d {d} -t {t} -j {j} -i {i} -f {f} \
    -e {e} -m {m} -z {z}".format(u=utilitiesPath, s=RSCRIPT_ML, c=args.sqanti_class, \
    o=args.output, d=args.dir, t=args.percent_training, j=args.threshold, i=args.intrapriming ,\
    f=args.force_fsm_in, e=args.filter_mono_exonic, m=args.intermediate_files, r=args.remove_columns, z=args.max_class_size)
    
    report_cmd=RSCRIPTPATH + " {u}/{s} -d {d} -o {o} -u {u} -f ml ".format(u=utilitiesPath, s=RSCRIPT_REPORT, \
    o=args.output, d=args.dir)

    if args.TP is not None:
        if not os.path.isfile(args.TP):
            print("ERROR: {0} doesn't exist. Abort!".format(args.TP), file=sys.stderr)
            sys.exit(-1)
        else:
            cmd = cmd + " -p {0}".format(args.TP)
    if args.TN is not None:
        if not os.path.isfile(args.TN):
            print("ERROR: {0} doesn't exist. Abort!".format(args.TN), file=sys.stderr)
            sys.exit(-1)
        else:
            cmd = cmd + " -n {0}".format(args.TN)
    if args.remove_columns is not None:
        if not os.path.isfile(args.remove_columns):
            print("ERROR: {0} doesn't exist. Abort!".format(args.remove_columns), file=sys.stderr)
            sys.exit(-1)
        else:
            cmd = cmd + " -r {0}".format(args.remove_columns)
    print(cmd)
    subprocess.call(cmd, shell=True)
    if not args.skip_report:
      subprocess.call(report_cmd, shell=True)
    # After running ML R code, an inclusion list will be generated. Those IDs must be passed to the filter files function
    inclusion_list = args.dir + "/" + args.output + "_inclusion-list.txt"
    seqs_to_keep = set(line.strip() for line in open(inclusion_list))
    return(seqs_to_keep, inclusion_list)

def run_rules(args):
    cmd = RSCRIPTPATH + " {u}/{s} -c {c} -o {o} -d {d} -j {j}".format(u=utilitiesPath, \
     s=RSCRIPT_RULES, c=args.sqanti_class, o=args.output, d=args.dir, j=args.json_filter)

    report_cmd=RSCRIPTPATH + " {u}/{s} -d {d} -o {o} -u {u} -f rules ".format(u=utilitiesPath, s=RSCRIPT_REPORT, \
    o=args.output, d=args.dir)
    
    if args.json_filter is not None:
        if not os.path.isfile(args.json_filter):
            print("ERROR: {0} doesn't exist. Abort!".format(args.json_filter), file=sys.stderr)
            sys.exit(-1)
    print(cmd)
    subprocess.call(cmd, shell=True)
    if not args.skip_report:
      subprocess.call(report_cmd, shell=True)
    # After running Rules Filter code, an inclusion list will be generated. Those IDs must be passed to the filter files function
    inclusion_list = args.dir + "/" + args.output + "_inclusion-list.txt"
    seqs_to_keep = set(line.strip() for line in open(inclusion_list))
    return(seqs_to_keep, inclusion_list)


def main():
    parser = argparse.ArgumentParser(description="Filtering of Isoforms based on SQANTI3 attributes.\
\nChoose between a rules filter or a Machine-Learning based filter.")
### Common arguments for both modes
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('sqanti_class', help='\t\tSQANTI3 QC classification file.')
    common.add_argument('--isoAnnotGFF3', help='\t\tisoAnnotLite GFF3 file to be filtered')
    common.add_argument('--isoforms', help='\t\tfasta/fastq isoform file to be filtered')
    common.add_argument('--gtf', help='\t\tGTF file to be filtered')
    common.add_argument('--sam', help='\t\tSAM alignment of the input fasta/fastq')
    common.add_argument('--faa', help='\t\tORF prediction faa file to be filtered by SQANTI3')
    common.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
    common.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
    common.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    common.add_argument("--skip_report", action="store_true", default=False, help='\t\tSkip creation of a report about the filtering')
    subparsers = parser.add_subparsers(dest='subcommand')

### Rules filter arguments
    rules = subparsers.add_parser('rules', parents=[common], description="Rules filter selected")
    rules.add_argument('-j', "--json_filter", default=default_json, help='\t\tJSON file where filtering rules are expressed. Rules must be set taking into account that attributes described in the filter will be present in those isoforms that should be kept. Default: utilities/filter/filter_default.json')

### ML filter arguments
    ml = subparsers.add_parser('ML', parents=[common],  description='ML filter selected')
    ml.add_argument('-t','--percent_training', type=float, default=0.8, \
    help="Proportion of the data that goes to training (parameter p of the function createDataPartition). \
    \nDefault: 0.8")
    ml.add_argument('-p', '--TP', \
    help="Path to file containing the list of the TP transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.")
    ml.add_argument('-n', '--TN', \
    help="Path to file containing the list of the TN transcripts, one ID by line, no header (optional). If not supplied, it will be generated from input data.")
    ml.add_argument('-j', '--threshold', type=float, default=0.7, \
    help="Machine Learning probability threshold to classify transcripts as positive isoforms.")
    ml.add_argument('-f', '--force_fsm_in', default=False, \
    help="If activated, forces retaining only multi-exon transcripts, all mono-exon isoforms will be automatically removed.")
    ml.add_argument('--intermediate_files', default=False, \
    help="If activated, outputs ML filter intermediate files.")
    ml.add_argument('-r','--remove_columns', \
    help="Path to single-column file (no header) containing the names of the columns in SQ3's classification.txt file that are to be excluded during random forest training (optional).")
    ml.add_argument('-z', '--max_class_size', type=int , default=3000, \
    help="Maximum number of isoforms to include in True Positive and True Negative sets. TP and TN sets will be downsized to this value if they are larger.")
    ml.add_argument('-i',"--intrapriming", type=float, default=60, help='\t\tAdenine percentage at genomic 3\' end to flag an isoform as intra-priming (default: 60 )')
    ml.add_argument("-e","--filter_mono_exonic", action="store_true", default=False, help='\t\tFilter out all mono-exonic transcripts (default: OFF)')

    args = parser.parse_args()

### Checking presence of files. Common arguments
    args.sqanti_class = os.path.abspath(args.sqanti_class)
    if not os.path.isfile(args.sqanti_class):
        print("ERROR: {0} doesn't exist. Abort!".format(args.sqanti_class), file=sys.stderr)
        sys.exit(-1)

    if args.isoforms is not None and not os.path.exists(args.isoforms):
        print("ERROR: {0} doesn't exist. Abort!".format(args.isoform), file=sys.stderr)
        sys.exit(-1)

    if args.gtf is not None and not os.path.exists(args.gtf):
        print("ERROR: {0} doesn't exist. Abort!".format(args.gtf_file), file=sys.stderr)
        sys.exit(-1)

    if args.sam is not None and not os.path.exists(args.sam):
        print("ERROR: {0} doesn't exist. Abort!".format(args.sam), file=sys.stderr)
        sys.exit(-1)

    if args.faa is not None and not os.path.exists(args.faa):
        print("ERROR: {0} doesn't exist. Abort!".format(args.faa), file=sys.stderr)
        sys.exit(-1)
### Define output dir and output name in case it was not defined
    if args.dir is None:
        args.dir=os.path.dirname(args.sqanti_class)
        print("Output directory not defined. All the outputs will be stored at {0} directory".format(args.dir), file=sys.stderr)
    else:
        if not os.path.exists(args.dir):
            os.makedirs(args.dir)
    if args.output is None:
        args.output=args.sqanti_class[args.sqanti_class.rfind("/")+1:args.sqanti_class.rfind("_classification.txt")]
        print("Output name not defined. All the outputs will have the prefix {0}".format(args.output), file=sys.stderr)

### Define TRUE or FALSE for boolean arguments

### Checking presence of files for ML. Check arguments --> If ok run ML
    
    print("\nRunning SQANTI3 filtering...\n", file=sys.stdout)
    
    if args.subcommand == 'ML':
        if args.TP is not None and not os.path.exists(args.TP):
            print("ERROR: {0} doesn't exist. Abort!".format(args.TP), file=sys.stderr)
            sys.exit(-1)
        if args.TN is not None and not os.path.exists(args.TN):
            print("ERROR: {0} doesn't exist. Abort!".format(args.TN), file=sys.stderr)
            sys.exit(-1)
        if args.remove_columns is not None and not os.path.exists(args.remove_columns):
            print("ERROR: {0} doesn't exist. Abort!".format(args.remove_columns), file=sys.stderr)
            sys.exit(-1)
        if args.percent_training < 0 or args.percent_training > 1.:
            print("ERROR: --percent_training must be between 0-1, instead given {0}! Abort!".format(args.intrapriming), file=sys.stderr)
            sys.exit(-1)
        if args.intrapriming < 25 or args.intrapriming > 100:
            print("ERROR: --intrapriming must be between 25-100, instead given {0}! Remember to use the percentage value. Abort!".format(args.intrapriming), file=sys.stderr)


        ids, inclusion_file = run_ML(args)
        filter_files(args, ids, inclusion_file)

### Checking presence of files for Rules. Check arguments --> If ok run Rules

    if args.subcommand == 'rules':
        if args.json_filter is not None and not os.path.exists(args.json_filter):
            print("ERROR: {0} doesn't exist. Abort!".format(args.json_filter), file=sys.stderr)
            sys.exit(-1)

        ids, inclusion_file = run_rules(args)
        filter_files(args, ids, inclusion_file)

if __name__ == "__main__":
    main()
