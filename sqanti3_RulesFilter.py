#!/usr/bin/env python3
__author__  = "etseng@pacb.com"
__version__ = '3.0'   # Python 3.7 syntax!

"""
Lightweight filtering of SQANTI by using .classification.txt output

Only keep Iso-Seq isoforms if:
The isoform is FSM, ISM, or NIC and (does not have intrapriming or has polyA_motif)
The isoform is NNC, does not have intrapriming/or polyA motif, not RT-switching, and all junctions are either all canonical or short-read-supported
The isoform is antisense, intergenic, genic, does not have intrapriming/or polyA motif, not RT-switching, and all junctions are either all canonical or short-read-supported
"""

import os, sys, argparse, subprocess
import distutils.spawn
from csv import DictReader, DictWriter
from Bio import SeqIO
from cupcake.io.BioReaders import GMAPSAMReader
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

utilitiesPath =  os.path.dirname(os.path.realpath(__file__))+"/utilities/"
RSCRIPTPATH = distutils.spawn.find_executable('Rscript')
RSCRIPT_REPORT = 'SQANTI3_report.R'

if os.system(RSCRIPTPATH + " --version")!=0:
    print("Rscript executable not found! Abort!", file=sys.stderr)
    sys.exit(-1)


CATEGORY_DICT = {'full-splice_match': 'FSM',
                 'incomplete-splice_match': 'ISM',
                 'novel_in_catalog': 'NIC',
                 'novel_not_in_catalog': 'NNC',
                 'antisense': 'AS',
                 'intergenic': 'intergenic',
                 'genic_intron': 'intron',
                 'genic': 'genic',
                 'fusion': 'fusion'}

def sqanti_filter_lite(args):

    fafq_type = 'fasta'
    with open(args.isoforms) as h:
        if h.readline().startswith('@'): fafq_type = 'fastq'

    prefix = args.sqanti_class[:args.sqanti_class.rfind('.')]

    fcsv = open(prefix + '.filtered_lite_reasons.txt', 'w')
    fcsv.write("# classification: {0}\n".format(args.sqanti_class))
    fcsv.write("# isoform: {0}\n".format(args.isoforms))
    fcsv.write("# intrapriming cutoff: {0}\n".format(args.intrapriming))
    fcsv.write("# min_cov cutoff: {0}\n".format(args.min_cov))
    fcsv.write("filtered_isoform,reason\n")

    fout = open(prefix + '.filtered_lite.' + fafq_type, 'w')

    seqids_to_keep = set()
    total_count = 0
    for r in DictReader(open(args.sqanti_class), delimiter='\t'):
        total_count += 1
        filter_flag, filter_msg = False, ""
        percA = float(r['perc_A_downstream_TTS']) / 100
        assert 0 <= percA <= 1
        runA = 0
        while runA < len(r['seq_A_downstream_TTS']):
            if r['seq_A_downstream_TTS'][runA]!='A':
                break
            runA += 1
        min_cov = float(r['min_cov']) if r['min_cov']!='NA' else None
        num_exon = int(r['exons'])
        is_RTS = r['RTS_stage'] == 'TRUE'
        is_canonical = r['all_canonical']=='canonical'
        is_monoexonic = (num_exon == 1)

        cat = CATEGORY_DICT[r['structural_category']]

        potential_intrapriming = (percA >= args.intrapriming or runA >= args.runAlength) and \
                                 r['polyA_motif'] == 'NA' and \
                                 (r['diff_to_gene_TSS'] == 'NA' or abs(
                                     int(r['diff_to_gene_TTS'])) > args.max_dist_to_known_end)


        if cat in ['FSM']:
            if potential_intrapriming:
                filter_flag, filter_msg = True, "IntraPriming"
            elif args.filter_mono_exonic and is_monoexonic:
                filter_flag, filter_msg = True, "Mono-Exonic"
        else:
            if potential_intrapriming:
                filter_flag, filter_msg = True, "IntraPriming"
            elif args.filter_mono_exonic and is_monoexonic:
                filter_flag, filter_msg = True, "Mono-Exonic"
            elif is_RTS:
                filter_flag, filter_msg = True, "RTSwitching"
            elif (not is_canonical) and (min_cov is None or (min_cov is not None and min_cov < args.min_cov)):
                filter_flag, filter_msg = True, "LowCoverage/Non-Canonical"

        if not filter_flag:
            seqids_to_keep.add(r['isoform'])
        else:
            fcsv.write("{0},{1}\n".format(r['isoform'], filter_msg))

    print("{0} isoforms read from {1}. {2} to be kept.".format(total_count, args.sqanti_class, len(seqids_to_keep)), file=sys.stdout)

    if not args.skipFaFq:
        for r in SeqIO.parse(open(args.isoforms), fafq_type):
            if r.id in seqids_to_keep:
                SeqIO.write(r, fout, fafq_type)
        fout.close()
        print("Output written to: {0}".format(fout.name), file=sys.stdout)

    # write out a new .classification.txt, .junctions.txt
    outputClassPath = prefix + '.filtered_lite_classification.txt'
    with open(outputClassPath, 'w') as f:
        reader = DictReader(open(args.sqanti_class), delimiter='\t')
        writer = DictWriter(f, reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for r in reader:
            if r['isoform'] in seqids_to_keep:
                writer.writerow(r)
        print("Output written to: {0}".format(f.name), file=sys.stdout)

    if not args.skipJunction:
        outputJuncPath = prefix + '.filtered_lite_junctions.txt'
        with open(outputJuncPath, 'w') as f:
            reader = DictReader(open(args.sqanti_class.replace('_classification', '_junctions')), delimiter='\t')
            writer = DictWriter(f, reader.fieldnames, delimiter='\t')
            writer.writeheader()
            for r in reader:
                if r['isoform'] in seqids_to_keep:
                    writer.writerow(r)
            print("Output written to: {0}".format(f.name), file=sys.stdout)

    if not args.skipGTF:
        outputGTF = prefix + '.filtered_lite.gtf'
        with open(outputGTF, 'w') as f:
            for r in collapseGFFReader(args.gtf_file):
                if r.seqid in seqids_to_keep:
                    write_collapseGFF_format(f, r)
            print("Output written to: {0}".format(f.name), file=sys.stdout)

    if args.sam is not None:
        outputSam = prefix + '.filtered_lite.sam'
        with open(outputSam, 'w') as f:
            reader = GMAPSAMReader(args.sam, True)
            f.write(reader.header)
            for r in reader:
                if r.qID in seqids_to_keep:
                    f.write(r.record_line + '\n')
            print("Output written to: {0}".format(f.name), file=sys.stdout)

    if args.faa is not None:
        outputFAA = prefix + '.filtered_lite.faa'
        with open(outputFAA, 'w') as f:
            for r in SeqIO.parse(open(args.faa), 'fasta'):
                if r.id in seqids_to_keep:
                    f.write(">{0}\n{1}\n".format(r.description, r.seq))
        print("Output written to: {0}".format(f.name), file=sys.stdout)



    if args.report != 'skip':
        print("**** Generating SQANTI3 report....", file=sys.stderr)
        cmd = RSCRIPTPATH + " {d}/{f} {c} {j} {p} {d} {a} {b}".format(d=utilitiesPath, f=RSCRIPT_REPORT, c=outputClassPath, j=outputJuncPath, p="mock", a=args.saturation, b=args.report)
        if subprocess.check_call(cmd, shell=True)!=0:
            print("ERROR running command: {0}".format(cmd), file=sys.stderr)
            sys.exit(-1)


def main():
    parser = argparse.ArgumentParser(description="Filtering of Isoforms based on SQANTI3 attributes")
    parser.add_argument('sqanti_class', help='\t\tSQANTI classification output file.')
    parser.add_argument('isoforms', help='\t\tfasta/fastq isoform file to be filtered by SQANTI3')
    parser.add_argument('gtf_file', help='\t\tGTF of the input fasta/fastq')
    parser.add_argument('--sam', help='\t\t(Optional) SAM alignment of the input fasta/fastq')
    parser.add_argument('--faa', help="\t\t(Optional) ORF prediction faa file to be filtered by SQANTI3")
    parser.add_argument('-a',"--intrapriming", type=float, default=0.6, help='\t\tAdenine percentage at genomic 3\' end to flag an isoform as intra-priming (default: 0.6)')
    parser.add_argument('-r', "--runAlength", type=int, default=6, help='\t\tContinuous run-A length at genomic 3\' end to flag an isoform as intra-priming (default: 6)')
    parser.add_argument('-m',"--max_dist_to_known_end", type=int, default=50, help="\t\tMaximum distance to an annotated 3' end to preserve as a valid 3' end and not filter out (default: 50bp)")
    parser.add_argument("-c", "--min_cov", type=int, default=3, help="\t\tMinimum junction coverage for each isoform (only used if min_cov field is not 'NA'), default: 3")
    parser.add_argument("--filter_mono_exonic", action="store_true", default=False, help='\t\tFilter out all mono-exonic transcripts (default: OFF)')
    parser.add_argument("--skipGTF", action="store_true", default=False, help='\t\tSkip output of GTF')
    parser.add_argument("--skipFaFq", action="store_true", default=False, help='\t\tSkip output of isoform fasta/fastq')
    parser.add_argument("--skipJunction", action="store_true", default=False, help='\t\tSkip output of junctions file')
    parser.add_argument("--saturation", action="store_true", default=False, help='\t\tInclude saturation curves into report')
    parser.add_argument("--report", choices=['html', 'pdf', 'both', 'skip'], default='html', help='\t\tselect report format\t\t--html\t\t--pdf\t\t--both\t\t--skip')
    #parser.add_argument("--always_keep_canonical", default=False, action="store_true", help="Always keep isoforms with all canonical junctions, regardless of other criteria. (default: False)")
    parser.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))

    args = parser.parse_args()

    if args.intrapriming < 0.25 or args.intrapriming > 1.:
        print("ERROR: --intrapriming must be between 0.25-1, instead given {0}! Abort!".format(args.intrapriming), file=sys.stderr)
        sys.exit(-1)
    if args.runAlength < 4 or args.runAlength > 20:
        print("ERROR: --runAlength must be between 4-20, instead given {0}! Abort!".format(args.runAlength), file=sys.stderr)
        sys.exit(-1)

    args.sqanti_class = os.path.abspath(args.sqanti_class)
    if not os.path.isfile(args.sqanti_class):
        print("ERROR: {0} doesn't exist. Abort!".format(args.sqanti_class), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.isoforms):
        print("ERROR: {0} doesn't exist. Abort!".format(args.isoform), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.gtf_file):
        print("ERROR: {0} doesn't exist. Abort!".format(args.gtf_file), file=sys.stderr)
        sys.exit(-1)

    if args.sam is not None and not os.path.exists(args.sam):
        print("ERROR: {0} doesn't exist. Abort!".format(args.sam), file=sys.stderr)
        sys.exit(-1)

    if args.faa is not None and not os.path.exists(args.faa):
        print("ERROR: {0} doesn't exist. Abort!".format(args.faa), file=sys.stderr)
        sys.exit(-1)

    print("\nRunning SQANTI3 filtering...\n", file=sys.stdout)

    sqanti_filter_lite(args)



if __name__ == "__main__":
    main()
