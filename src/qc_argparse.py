import argparse
from src.config import __version__,__author__

def qc_argparse():
    ap = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms",
                                 exit_on_error=True)
    # Required arguments
    apr = ap.add_argument_group("Required arguments")
    apr.add_argument('--isoforms',  required= True,help='Isoforms (FASTA/FASTQ) or GTF format. It is recommended to provide them in GTF format, but if it is needed to map the sequences to the genome use a FASTA/FASTQ file with the --fasta option.')
    apr.add_argument('--refGTF',  required= True, help='Reference annotation file (GTF format)')
    apr.add_argument('--refFasta',  required=True, help='Reference genome (Fasta format)')

    # Customization and filtering args
    apc = ap.add_argument_group("Customization and filtering")
    apc.add_argument("--min_ref_len", type=int, default=0, help="Minimum reference transcript length (default: 0 bp)")
    apc.add_argument("--force_id_ignore", action="store_true", help=" Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    apc.add_argument('--fasta', action='store_true', help='Use when running SQANTI by using as input a FASTA/FASTQ with the sequences of isoforms')
    apc.add_argument('--genename', action='store_true' ,help='Use gene_name tag from GTF to define genes. Default: gene_id used to define genes',)
    apc.add_argument('--novel_gene_prefix', default=None, help='Prefix for novel isoforms (default: None)')
    apc.add_argument('-s','--sites', default="ATAC,GCAG,GTAG", help='Set of splice sites to be considered as canonical, in a comma separated list. (default: %(default)s)')
    apc.add_argument('-w','--window', default=20, type=int, help='Size of the window in the genomic DNA screened for Adenine content downstream of TTS (default: %(default)s)')

    # Aligner and mapping options
    apa = ap.add_argument_group("Aligner and mapping options")
    #TODO: set a default aligner
    apa.add_argument("--aligner_choice", choices=['minimap2', 'deSALT', 'gmap', "uLTRA"], default='minimap2', help="Select your aligner of choice: minimap2, deSALT, gmap, uLTRA (default: %(default)s)")
    apa.add_argument('-x','--gmap_index', help='Path and prefix of the reference index created by gmap_build. Mandatory if using GMAP .')

    # ORF prediction
    apo = ap.add_argument_group("ORF prediction")
    apo.add_argument("--skipORF", action="store_true", help="Skip ORF prediction (to save time)")
    apo.add_argument("--orf_input",  help="Input fasta to run ORF on. By default, ORF is run on genome-corrected fasta - this overrides it. If input is fusion (--is_fusion), this must be provided for ORF prediction.")

    # Orthogonal data inputs
    apod = ap.add_argument_group("Orthogonal data inputs")
    apod.add_argument('--short_reads',  help='File Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.')
    apod.add_argument('--SR_bam',  help=' Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome')
    apod.add_argument('--CAGE_peak',  help='FANTOM5 Cage Peak (BED format, optional)')
    apod.add_argument("--polyA_motif_list",  help="Ranked list of polyA motifs (text, optional)")
    apod.add_argument("--polyA_peak",  help='PolyA Peak (BED format, optional)')
    apod.add_argument("--phyloP_bed",  help="PhyloP BED for conservation score (BED, optional)")
    apod.add_argument('-e','--expression',  help='Expression matrix (supported: Kallisto tsv)')
    apod.add_argument('-c','--coverage', help='Junction coverage files (provide a single file, comma-delmited filenames, or a file pattern, ex: "mydir/*.junctions").')
    apod.add_argument('-fl', '--fl_count', help='Full-length PacBio abundance file')

    # Functional annotation
    apf = ap.add_argument_group("Functional annotation")
    apf.add_argument('--isoAnnotLite', action='store_true', help='Run isoAnnot Lite to output a tappAS-compatible gff3 file')
    apf.add_argument('--gff3' , help='Precomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes')
    
    # Output options
    apout = ap.add_argument_group("Output options")
    apout.add_argument('-o','--output',default= "isoforms", help='Prefix for output files')
    apout.add_argument('-d','--dir',  default='./sqanti3_results', help='Directory for output files. (Default: Directory where the script was run.)')
    apout.add_argument("--saturation", action="store_true", default=False, help='Include saturation curves into report')
    apout.add_argument("--report", choices=['html', 'pdf', 'both', 'skip'], default='html', help=f"Select report format: {', '.join(['html', 'pdf', 'both', 'skip'])} (default: %(default)s)")
    apout.add_argument('--isoform_hits' , action='store_true', help=' Report all FSM/ISM isoform hits in a separate file')
    apout.add_argument('--ratio_TSS_metric' , choices=['max', 'mean', 'median', '3quartile'], default='max', help=' Define which statistic metric should be reported in the ratio_TSS column (default: %(default)s)')

    # Performance options
    app = ap.add_argument_group("Performance options")
    app.add_argument('-t', '--cpus', default=8, type=int, help='Number of threads used during alignment by aligners. (default: 10)')
    app.add_argument('-n', '--chunks', default=1, type=int, help='Number of chunks to split SQANTI3 analysis in for speed up (default: 1).')
    app.add_argument("-l","--log_level", default="INFO",choices=["ERROR","WARNING","INFO","DEBUG"],
                        help="Set the logging level %(default)s")

    # Optional arguments
    apoa = ap.add_argument_group("Optional arguments")
    apoa.add_argument("--is_fusion", action="store_true", help="Input are fusion isoforms, must supply GTF as input")
    apoa.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    apoa.add_argument('--tusco', choices=['human','mouse'], help='Generate a TUSCO benchmarking report for the given species (human or mouse)', required=False)

    return ap
