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
from src.argparse_utils import args_validation

def main():

    args = qc_argparse().parse_args()
    args = args_validation(args)


    # path and prefix for output files


    # Print out parameters so can be put into report PDF later
    args.doc = os.path.join(os.path.abspath(args.dir), args.output+".params.txt")
    print("Write arguments to {0}...".format(args.doc, file=sys.stdout))
    with open(args.doc, 'w') as f:
        f.write("Version\t" + __version__ + "\n")
        f.write("Input\t" + os.path.abspath(args.isoforms) + "\n")
        f.write("Annotation\t" + os.path.abspath(args.refGTF) + "\n")
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
