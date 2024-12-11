#!/usr/bin/env python3
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila
# Modified by Liz (etseng@pacb.com) as SQANTI2/3 versions
# Modified by Fran (francisco.pardo.palacios@gmail.com) currently as SQANTI3 version (05/15/2020)
# Modified by Pablo (pabloatienzalo@gmail.com)
__author__  = "etseng@pacb.com"
__version__ = '5.3.0'  # Python 3.7

import os, sys
import shutil
import argparse

# Import SQANTI3 modules
from src.qc_pipeline import run
from src.parallel import split_input_run, combine_split_runs, get_split_dir
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
    global utilitiesPath

    # TODO: Take the parser out of here
    # The arguments are divided into categories, based on their functionality to SQANTI3
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    parser.add_argument('isoforms', help='\tIsoforms (FASTA/FASTQ) or GTF format. It is recommended to provide them in GTF format, but if it is needed to map the sequences to the genome use a FASTA/FASTQ file with the --fasta option.')
    parser.add_argument('annotation', help='\t\tReference annotation file (GTF format)')
    parser.add_argument('genome', help='\t\tReference genome (Fasta format)')
    parser.add_argument("--min_ref_len", type=int, default=0, help="\t\tMinimum reference transcript length (default: 0 bp)")
    parser.add_argument("--force_id_ignore", action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    parser.add_argument("--aligner_choice", choices=['minimap2', 'deSALT', 'gmap', "uLTRA"], default='minimap2')
    parser.add_argument('--CAGE_peak', help='\t\tFANTOM5 Cage Peak (BED format, optional)')
    parser.add_argument("--polyA_motif_list", help="\t\tRanked list of polyA motifs (text, optional)")
    parser.add_argument("--polyA_peak", help='\t\tPolyA Peak (BED format, optional)')
    parser.add_argument("--phyloP_bed", help="\t\tPhyloP BED for conservation score (BED, optional)")
    parser.add_argument("--skipORF", default=False, action="store_true", help="\t\tSkip ORF prediction (to save time)")
    parser.add_argument("--is_fusion", default=False, action="store_true", help="\t\tInput are fusion isoforms, must supply GTF as input")
    parser.add_argument("--orf_input", help="\t\tInput fasta to run ORF on. By default, ORF is run on genome-corrected fasta - this overrides it. If input is fusion (--is_fusion), this must be provided for ORF prediction.")
    parser.add_argument('--fasta', help='\t\tUse when running SQANTI by using as input a FASTA/FASTQ with the sequences of isoforms', action='store_true')
    parser.add_argument('-e','--expression', help='\t\tExpression matrix (supported: Kallisto tsv)', required=False)
    parser.add_argument('-x','--gmap_index', help='\t\tPath and prefix of the reference index created by gmap_build. Mandatory if using GMAP unless -g option is specified.')
    parser.add_argument('-t', '--cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. (default: 10)')
    parser.add_argument('-n', '--chunks', default=1, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up (default: 1).')
    #parser.add_argument('-z', '--sense', help='\t\tOption that helps aligners know that the exons in you cDNA sequences are in the correct sense. Applicable just when you have a high quality set of cDNA sequences', required=False, action='store_true')
    parser.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
    parser.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
    parser.add_argument('-c','--coverage', help='\t\tJunction coverage files (provide a single file, comma-delmited filenames, or a file pattern, ex: "mydir/*.junctions").', required=False)
    parser.add_argument('-s','--sites', default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    parser.add_argument('-w','--window', default="20", help='\t\tSize of the window in the genomic DNA screened for Adenine content downstream of TTS', required=False, type=int)
    parser.add_argument('--genename', help='\t\tUse gene_name tag from GTF to define genes. Default: gene_id used to define genes', default=False, action='store_true')
    parser.add_argument('-fl', '--fl_count', help='\t\tFull-length PacBio abundance file', required=False)
    parser.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    parser.add_argument("--saturation", action="store_true", default=False, help='\t\tInclude saturation curves into report')
    parser.add_argument("--report", choices=['html', 'pdf', 'both', 'skip'], default='html', help='\t\tselect report format\t\t--html\t\t--pdf\t\t--both\t\t--skip')
    parser.add_argument('--isoAnnotLite' , help='\t\tRun isoAnnot Lite to output a tappAS-compatible gff3 file',required=False, action='store_true' , default=False)
    parser.add_argument('--gff3' , help='\t\tPrecomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes',required=False)
    parser.add_argument('--short_reads', help='\t\tFile Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.', required=False)
    parser.add_argument('--SR_bam' , help='\t\t Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome', required=False)
    parser.add_argument('--isoform_hits' , help='\t\t Report all FSM/ISM isoform hits in a separate file', required=False, default = False, action='store_true')
    parser.add_argument('--ratio_TSS_metric' , help='\t\t Define which statistic metric should be reported in the ratio_TSS column', choices=['max', 'mean', 'median', '3quartile'], default='max')

    args = parser.parse_args()
    # Arguments checks
    if args.is_fusion:
        if args.orf_input is None:
            print("WARNING: Currently if --is_fusion is used, no ORFs will be predicted. Supply --orf_input if you want ORF to run!", file=sys.stderr)
            args.skipORF = True
        if args.fasta:
            print("ERROR: if --is_fusion is on, must supply GTF as input", file=sys.stderr)
            sys.exit(1)

    if args.gff3 is not None:
        args.gff3 = os.path.abspath(args.gff3)
        if not os.path.isfile(args.gff3):
            print("ERROR: Precomputed tappAS GFF3 annoation file {0} doesn't exist. Abort!".format(args.genome), file=sys.stderr)
            sys.exit(1)

    if args.expression is not None:
        if os.path.isdir(args.expression)==True:
            print("Expression files located in {0} folder".format(args.expression), file=sys.stderr)
        else:
            for f in args.expression.split(','):
                if not os.path.exists(f):
                        print("Expression file {0} not found. Abort!".format(f), file=sys.stderr)
                        sys.exit(1)


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
