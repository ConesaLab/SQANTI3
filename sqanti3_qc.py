#!/usr/bin/env python3
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila
# Modified by Liz (etseng@pacb.com) as SQANTI2/3 versions
# Modified by Fran (francisco.pardo.palacios@gmail.com) currently as SQANTI3 version (05/15/2020)
# Modified by Pablo (pabloatienzalo@gmail.com)

import os, sys
import shutil

# Import SQANTI3 modules
from src.qc_argparse import qc_argparse, args_validation
from src.qc_pipeline import run
from src.parallel import split_input_run, combine_split_runs, get_split_dir
from src.config import __version__
# TODO: Change import errors to propper logging
# try:
#     from Bio import SeqIO
# except ImportError:
#     print(f"ImportError: {e}", file=sys.stderr)
#     print("Unable to import Biopython! Please make sure Biopython is installed.", file=sys.stderr)
#     sys.exit(1)

# try:
#     from bx.intervals import Interval, IntervalTree
# except ImportError:
#     print(f"ImportError: {e}", file=sys.stderr)
#     print("Unable to import bx-python! Please make sure bx-python is installed.", file=sys.stderr)
#     sys.exit(1)

# try:
#     from BCBio import GFF as BCBio_GFF
# except ImportError:
#     print(f"ImportError: {e}", file=sys.stderr)
#     print("Unable to import BCBio! Please make sure bcbiogff is installed.", file=sys.stderr)
#     sys.exit(1)


def main():

    args = qc_argparse()
    args_validation(args)


    # path and prefix for output files
    if args.output is None:
        args.output = os.path.splitext(os.path.basename(args.isoforms))[0]

    if args.dir is None:
        args.dir = os.getcwd()
    else:
        args.dir = os.path.abspath(args.dir)
        if os.path.isdir(args.dir):
            print("WARNING: output directory {0} already exists. Overwriting!".format(args.dir), file=sys.stderr)
        else:
            os.makedirs(args.dir)

    args.genome = os.path.abspath(args.genome)
    if not os.path.isfile(args.genome):
        print("ERROR: genome fasta {0} doesn't exist. Abort!".format(args.genome), file=sys.stderr)
        sys.exit()

    args.isoforms = os.path.abspath(args.isoforms)
    if not os.path.isfile(args.isoforms):
        print("ERROR: Input isoforms {0} doesn't exist. Abort!".format(args.isoforms), file=sys.stderr)
        sys.exit()

    if args.fasta:
        if args.aligner_choice == 'gmap':
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                print("GMAP index {0} doesn't exist! Abort.".format(args.gmap_index), file=sys.stderr)
                sys.exit()
        elif args.aligner_choice == 'deSALT':
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                print("deSALT index {0} doesn't exist! Abort.".format(args.gmap_index), file=sys.stderr)
                sys.exit()

        print("Cleaning up isoform IDs...", file=sys.stderr)
        from src.helpers import rename_isoform_seqids
        args.isoforms = rename_isoform_seqids(args.isoforms, args.force_id_ignore)
        print("Cleaned up isoform fasta file written to: {0}".format(args.isoforms), file=sys.stderr)


    args.annotation = os.path.abspath(args.annotation)
    if not os.path.isfile(args.annotation):
        print("ERROR: Annotation doesn't exist. Abort!".format(args.annotation), file=sys.stderr)
        sys.exit()

    #if args.aligner_choice == "gmap":
    #    args.sense = "sense_force" if args.sense else "auto"
    #elif args.aligner_choice == "minimap2":
    #    args.sense = "f" if args.sense else "b"
    ## (Liz) turned off option for --sense, always TRUE
    if args.aligner_choice == "gmap":
        args.sense = "sense_force"
    elif args.aligner_choice == "minimap2":
        args.sense = "f"
    #elif args.aligner_choice == "deSALT":  #deSALT does not support this yet
    #    args.sense = "--trans-strand"


    args.novel_gene_prefix = None
    # Print out parameters so can be put into report PDF later
    args.doc = os.path.join(os.path.abspath(args.dir), args.output+".params.txt")
    print("Write arguments to {0}...".format(args.doc, file=sys.stdout))
    with open(args.doc, 'w') as f:
        f.write("Version\t" + __version__ + "\n")
        f.write("Input\t" + os.path.abspath(args.isoforms) + "\n")
        f.write("Annotation\t" + os.path.abspath(args.annotation) + "\n")
        f.write("Genome\t" + os.path.abspath(args.genome) + "\n")
        f.write("MinRefLength\t"+ str(args.min_ref_len) + "\n")
        f.write("ForceIdIgnore\t"+str(args.force_id_ignore) + "\n")
        f.write("Aligner\t" + str(args.aligner_choice) + "\n")
        f.write("FLCount\t" + (os.path.abspath(args.fl_count) if args.fl_count is not None else "NA") + "\n")
        f.write("Expression\t" + (os.path.abspath(args.expression) if args.expression is not None else "NA") + "\n")
        f.write("Junction\t" + (os.path.abspath(args.coverage) if args.coverage is not None else "NA") + "\n")
        f.write("CAGEPeak\t" + (os.path.abspath(args.CAGE_peak)  if args.CAGE_peak is not None else "NA") + "\n")
        f.write("PolyAMotif\t" + (os.path.abspath(args.polyA_motif_list) if args.polyA_motif_list is not None else "NA") + "\n")
        f.write("PolyAPeak\t" + (os.path.abspath(args.polyA_peak)  if args.polyA_peak is not None else "NA") + "\n")
        f.write("IsFusion\t" + str(args.is_fusion) + "\n")
        f.write("PhyloP\t" + (os.path.abspath(args.phyloP_bed)  if args.phyloP_bed is not None else "NA") + "\n")
        f.write("SkipORF\t" + str(args.skipORF) + "\n")
        f.write("ORFInput\t" + (os.path.abspath(args.orf_input) if args.orf_input is not None else "NA" ) + "\n" )
        f.write("FASTAused\t" + str(args.fasta) +"\n")
        f.write("Expression\t" + (os.path.abspath(args.expression) if args.expression is not None else "NA" ) + "\n")
        f.write("GMAPindex\t" + (os.path.abspath(args.gmap_index) if args.gmap_index is not None else "NA" ) + "\n")
        f.write("OutputPrefix\t" + str(args.output) + "\n")
        f.write("OutputDirectory\t" + os.path.abspath(args.dir) + "\n")
        f.write("Coverage\t" + (os.path.abspath(args.coverage) if args.coverage is not None else "NA") + "\n" )
        f.write("CanonicalSites\t" + str(args.sites) + "\n")
        f.write("PostTTSWindow\t" + str(args.window) + "\n")
        f.write("GeneName\t" + str(args.genename) + "\n")
        f.write("ReportType\t" + str(args.report) + "\n")
        f.write("RunIsoAnnotLite\t" + str(args.isoAnnotLite) + "\n")
        f.write("isoAnnotGFF3\t" + (os.path.abspath(args.gff3) if args.gff3 is not None else "NA") + "\n")
        f.write("ShortReads\t" + (os.path.abspath(args.short_reads) if args.short_reads is not None else "NA") + "\n")
        f.write("ShortReadsBAMs\t" + (os.path.abspath(args.SR_bam) if args.SR_bam is not None else "NA") + "\n")
        f.write("ratioTSSmetric\t" + str(args.ratio_TSS_metric) + "\n")

    # Running functionality based on the chunks
    print("**** Running SQANTI3...", file=sys.stdout)
    
    if args.chunks == 1:
        run(args)

    else:
        split_dirs = split_input_run(args)
        combine_split_runs(args, split_dirs)
        SPLIT_ROOT_DIR = get_split_dir(args)
        shutil.rmtree(SPLIT_ROOT_DIR)

    if args.isoAnnotLite:
        from src.helpers import get_corr_filenames, get_class_junc_filenames
        from src.qc_pipeline import run_isoAnnotLite 
        corrGTF, corrSAM, corrFASTA, corrORF , corrCDS_GTF_GFF = get_corr_filenames(args)
        outputClassPath, outputJuncPath = get_class_junc_filenames(args)
        run_isoAnnotLite(corrGTF, outputClassPath, outputJuncPath, args.output, args.gff3)

if __name__ == "__main__":
    main()
