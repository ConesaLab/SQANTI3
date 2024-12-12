import argparse, sys, os
from .config import __version__,__author__


def qc_argparse():
    ap = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    # Required arguments
    apr = ap.add_argument_group("Required arguments")
    apr.add_argument('isoforms', help='\tIsoforms (FASTA/FASTQ) or GTF format. It is recommended to provide them in GTF format, but if it is needed to map the sequences to the genome use a FASTA/FASTQ file with the --fasta option.')
    apr.add_argument('annotation', help='\t\tReference annotation file (GTF format)')
    apr.add_argument('genome', help='\t\tReference genome (Fasta format)')

    # Customization and filtering args
    apc = ap.add_argument_group("Customization and filtering")
    apc.add_argument("--min_ref_len", type=int, default=0, help="\t\tMinimum reference transcript length (default: 0 bp)")
    apc.add_argument("--force_id_ignore", action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    apc.add_argument('--fasta', help='\t\tUse when running SQANTI by using as input a FASTA/FASTQ with the sequences of isoforms', action='store_true')
    apc.add_argument('--genename', help='\t\tUse gene_name tag from GTF to define genes. Default: gene_id used to define genes', default=False, action='store_true')
    apc.add_argument('--short_reads', help='\t\tFile Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.', required=False)
    apc.add_argument('--SR_bam' , help='\t\t Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome', required=False)

    # Aligner and mapping options
    apa = ap.add_argument_group("Aligner and mapping options")
    #TODO: set a default aligner
    apa.add_argument("--aligner_choice", choices=['minimap2', 'deSALT', 'gmap', "uLTRA"], default='minimap2')
    apa.add_argument('-x','--gmap_index', help='\t\tPath and prefix of the reference index created by gmap_build. Mandatory if using GMAP unless -g option is specified.')
    apa.add_argument('-s','--sites', default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)

    # ORF prediction
    apo = ap.add_argument_group("ORF prediction")
    apo.add_argument("--skipORF", default=False, action="store_true", help="\t\tSkip ORF prediction (to save time)")
    apo.add_argument("--orf_input", help="\t\tInput fasta to run ORF on. By default, ORF is run on genome-corrected fasta - this overrides it. If input is fusion (--is_fusion), this must be provided for ORF prediction.")

    # Functional annotation
    apf = ap.add_argument_group("Functional annotation")
    apf.add_argument('--CAGE_peak', help='\t\tFANTOM5 Cage Peak (BED format, optional)')
    apf.add_argument("--polyA_motif_list", help="\t\tRanked list of polyA motifs (text, optional)")
    apf.add_argument("--polyA_peak", help='\t\tPolyA Peak (BED format, optional)')
    apf.add_argument("--phyloP_bed", help="\t\tPhyloP BED for conservation score (BED, optional)")

    # Output options
    apout = ap.add_argument_group("Output options")
    apout.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
    apout.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
    apout.add_argument("--saturation", action="store_true", default=False, help='\t\tInclude saturation curves into report')
    apout.add_argument("--report", choices=['html', 'pdf', 'both', 'skip'], default='html', help='\t\tselect report format\t\t--html\t\t--pdf\t\t--both\t\t--skip')
    apout.add_argument('--isoform_hits' , help='\t\t Report all FSM/ISM isoform hits in a separate file', required=False, default = False, action='store_true')
    apout.add_argument('--ratio_TSS_metric' , help='\t\t Define which statistic metric should be reported in the ratio_TSS column', choices=['max', 'mean', 'median', '3quartile'], default='max')

    # Performance options
    app = ap.add_argument_group("Performance options")
    app.add_argument('-t', '--cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. (default: 10)')
    app.add_argument('-n', '--chunks', default=1, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up (default: 1).')

    # Optional arguments
    apm = ap.add_argument_group("Optional arguments")
    apm.add_argument("--is_fusion", default=False, action="store_true", help="\t\tInput are fusion isoforms, must supply GTF as input")
    apm.add_argument('-e','--expression', help='\t\tExpression matrix (supported: Kallisto tsv)', required=False)
    apm.add_argument('-c','--coverage', help='\t\tJunction coverage files (provide a single file, comma-delmited filenames, or a file pattern, ex: "mydir/*.junctions").', required=False)
    apm.add_argument('-w','--window', default="20", help='\t\tSize of the window in the genomic DNA screened for Adenine content downstream of TTS', required=False, type=int)
    apm.add_argument('-fl', '--fl_count', help='\t\tFull-length PacBio abundance file', required=False)
    apm.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    apm.add_argument('--isoAnnotLite' , help='\t\tRun isoAnnot Lite to output a tappAS-compatible gff3 file',required=False, action='store_true' , default=False)
    apm.add_argument('--gff3' , help='\t\tPrecomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes',required=False)

    return ap.parse_args()

def args_validation(args):
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
