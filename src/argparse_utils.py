import os, sys
import subprocess

from src.commands import GFFREAD_PROG
from src.module_logging import qc_logger, filter_logger, rescue_logger

def valid_file(filename,logger):
    if not os.path.isfile(filename):
        logger.error(f"File {filename} not found. Abort!")
        sys.exit(1)
    return filename

def valid_fasta(filename,logger):
    valid_file(filename,logger)
    extension = filename.split('.')[-1]
    valid_extensions = ['fasta', 'fa', 'fastq', 'fq', 'faa', 'fna']
    if extension not in valid_extensions:
        logger.error(f"File {filename} is not a FASTA file. Abort!")
        sys.exit(1)
    return filename

def valid_gtf(filename,logger):
    valid_file(filename,logger)
    if not filename.endswith('.gtf'):
        if not filename.endswith('.gff') or not filename.endswith('.gff3'):
            qc_logger.error(f"File {filename} is not a GTF file. Abort!")
            sys.exit(1)
        else:
            logger.warning("GTF file is in GFF3 format. Converting to GTF format.")
            # GFF to GTF (in case the user provides gff instead of gtf)
            gtf_name = filename.replace('.gff', '.gtf')
            # Use the run command?
            try:
                subprocess.call([GFFREAD_PROG, filename , '-T', '-o', gtf_name])
            except (RuntimeError, TypeError, NameError):
                logger.error(f'File {filename} without GTF/GFF format.')
                raise SystemExit(1)
            logger.info(f"GFF file converted to GTF format. New file: {gtf_name}")
            filename = gtf_name

    # Check if the GTF file is in the correct format
    ind = 0
    with open(filename) as isoforms_gtf:
        for line in isoforms_gtf:
            if line[0] != "#" and len(line.split("\t"))!=9:
                logger.error("Input isoforms file not in expected GTF format.")
                sys.exit()
            elif len(line.split("\t"))==9:
                ind += 1
        if ind == 0:
            logger.warning(f"GTF has {filename} no annotation lines.")

    return filename

def valid_bed(filename,logger):
    valid_file(filename,logger)
    if not filename.endswith('.bed'):
        logger.error(f"File {filename} is not a BED file. Abort!")
        sys.exit(1)
    return filename

def valid_dir(dirname,logger):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    else:
        contents = os.listdir(dirname)
        non_log_items = [item for item in contents if item != "logs" or not os.path.isdir(os.path.join(dirname, item))]
        if non_log_items:
            logger.warning(f"Output directory {dirname} already exists. Overwriting!")
    return dirname

# TODO: Add specific validation for Kallisto output
def valid_matrix(filename,logger):
    valid_file(filename,logger)
    if not filename.endswith('.tsv'):
        logger.error(f"File {filename} is not a TSV file. Abort!")
        sys.exit(1)

#TODO: Get a condition to see if it is a pacbio file
def valid_PacBio_abund(filename,logger):
    valid_file(filename,logger)
    if not filename.endswith('abundance.tsv'):
        logger.error(f"File {filename} is not a PacBio abundance file. Abort!")
        sys.exit(1)
    return filename

def valid_gff3(filename,logger):
    valid_file(filename,logger)
    if not filename.endswith('.gff3'):
        logger.error(f"File {filename} is not a GFF3 file. Abort!")
        sys.exit(1)
    return filename

def valid_sr(filename,logger):
    if not os.path.isdir(filename):
        if not filename.endswith('.fofn') and not filename.endswith('.bam'):
            logger.error(f"File {filename} is not a BAM file, a directory, or a FOFN file. Abort!")
            sys.exit(1)
        else:
            if not os.path.isfile(filename):
                logger.error(f"File {filename} not found. Abort!")
                sys.exit(1)
    return filename

def contains_gene_name(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            if 'gene_name' in line:
                return True
    return False

### Validation for the arguments

def qc_args_validation(args):
    # Required arguments
    valid_file(args.isoforms,qc_logger)
    if args.isoforms.endswith('.gtf') or args.isoforms.endswith('.gff'):
        args.fasta = False
    elif valid_fasta(args.isoforms,qc_logger):
        args.fasta = True
    else:
        qc_logger.error("Input isoforms must be in GTF, FASTA, or FASTQ format. Abort!")
        sys.exit(1)
    
    valid_gtf(args.refGTF,qc_logger)
    valid_fasta(args.refFasta,qc_logger)

    # Customization validation
    if args.short_reads is not None:
        valid_file(args.short_reads,qc_logger)
    if args.SR_bam is not None:
        valid_sr(args.SR_bam,qc_logger)

    # Mappers checks
    if args.fasta:
        if args.aligner_choice == 'gmap':
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                qc_logger.error(f"GMAP index {args.gmap_index} doesn't exist! Abort.")
                sys.exit()
        elif args.aligner_choice == 'deSALT':
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                qc_logger.error(f"deSALT index {args.gmap_index} doesn't exist! Abort.")
                sys.exit()

        qc_logger.info("Cleaning up isoform IDs...")
        from src.helpers import rename_isoform_seqids
        args.isoforms = rename_isoform_seqids(args.isoforms, args.force_id_ignore)
        qc_logger.info(f"Cleaned up isoform fasta file written to: {args.isoforms}")
    
    # ORF prediction checks
    if args.orf_input is not None:
        valid_fasta(args.orf_input,qc_logger)

    # Functional annotation checks
    if args.CAGE_peak is not None:
        valid_bed(args.CAGE_peak,qc_logger)
    if args.polyA_motif_list is not None:
        valid_file(args.polyA_motif_list,qc_logger)
    if args.polyA_peak is not None:
        valid_bed(args.polyA_peak,qc_logger)
    if args.phyloP_bed is not None:
        valid_bed(args.phyloP_bed,qc_logger)

    # Output options checks
    valid_dir(args.dir,qc_logger)


    if args.gff3 is not None:
        valid_gff3(args.gff3,qc_logger)

    # Check for gene_name attribute in GTF
    if args.isoAnnotLite or args.genename:
        if not contains_gene_name(args.isoforms):
            option = "--isoAnnotLite" if args.isoAnnotLite else "--genename"
            qc_logger.error("The 'gene_name' tag was not found in the input GTF file. SQANTI3 requires this tag to function properly.")
            qc_logger.error(f"Please include the 'gene_name' tag in the GTF, or omit the {option} option.")
            sys.exit(1)


    # Fusion isoforms checks
    if args.is_fusion:
        if args.orf_input is None:
            qc_logger.warning("Currently if --is_fusion is used, no ORFs will be predicted. Supply --orf_input if you want ORF to run!")
            args.skipORF = True
        if args.fasta:
            qc_logger.error("If --is_fusion is on, must supply GTF as input")
            sys.exit(1)
    # Expression checks TODO: Improve this checks
    if args.expression is not None:
        if os.path.isdir(args.expression):
            qc_logger.info(f"Expression files located in {args.expression} folder")
        else:
            for f in args.expression.split(','):
                valid_matrix(f,qc_logger)
    # Output prefix checks
    if args.output is None:
        args.output = os.path.splitext(os.path.basename(args.isoforms))[0]
    
    return args
## This sense argument is no longer present
   #if args.aligner_choice == "gmap":
    #    args.sense = "sense_force" if args.sense else "auto"
    #elif args.aligner_choice == "minimap2":
    #    args.sense = "f" if args.sense else "b"
    ## (Liz) turned off option for --sense, always TRUE
    # if args.aligner_choice == "gmap":
    #     args.sense = "sense_force"
    # elif args.aligner_choice == "minimap2":
    #     args.sense = "f"
    #elif args.aligner_choice == "deSALT":  #deSALT does not support this yet
    #    args.sense = "--trans-strand"

def filter_args_validation(args):
    # Mandatory + possible inputs
    valid_file(args.sqanti_class, filter_logger)
    if args.filter_isoforms is not None:
        valid_fasta(args.filter_isoforms, filter_logger)
    if args.filter_gtf is not None:
        valid_gtf(args.filter_gtf, filter_logger)
    if args.filter_sam is not None:
        valid_file(args.filter_sam, filter_logger)
    if args.filter_faa is not None:
        valid_fasta(args.filter_faa, filter_logger)

    valid_dir(args.dir,filter_logger)
    if args.subcommand == 'rules':
        valid_file(args.json_filter, filter_logger)
    if args.subcommand == 'ml':
        if args.TP is not None:
            valid_file(args.TP, filter_logger)
        if args.TN is not None:
            valid_file(args.TN, filter_logger)
        if args.remove_columns is not None:
            valid_file(args.remove_columns, filter_logger)
        if args.percent_training < 0 or args.percent_training > 1.:
            filter_logger.error(f"--percent_training must be between 0-1, instead given {args.percent_training}! Abort!")
            sys.exit(-1)
        if args.threshold < 0 or args.threshold > 1.:
            filter_logger.error(f"--threshold must be between 0-1, instead given {args.threshold}! Abort!")
            sys.exit(-1)
        if args.intrapriming < 25 or args.intrapriming > 100:
            filter_logger.error(f"--intrapriming must be between 25-100, instead given {args.intrapriming}! Remember to use the percentage value. Abort!")
            sys.exit(-1)
        if args.remove_columns is not None:
            valid_file(args.remove_columns, filter_logger)
    

def rescue_args_validation(args):
    valid_file(args.filter_class, rescue_logger)
    valid_dir(args.dir,rescue_logger)
    valid_gtf(args.refGTF,rescue_logger)
    valid_fasta(args.refFasta,rescue_logger)
    if args.rescue_gtf is not None:
        valid_gtf(args.rescue_gtf,rescue_logger)
    if args.mode == "full":
        try:
            valid_file(args.refClassif,rescue_logger)
        except:
            rescue_logger.error("When running the rescue in full mode, the reference classification file is mandatory.")
            sys.exit(1)
        try:
            valid_fasta(args.rescue_isoforms,rescue_logger)
        except:
            rescue_logger.error("When running the rescue in full mode, the corrected isoforms FASTA file is mandatory.")
            sys.exit(1)
        if args.subcommand == 'rules':
            valid_file(args.json_filter, rescue_logger)
        if args.subcommand == "ml":
            if args.threshold < 0 or args.threshold > 1.:
                rescue_logger.error(f"--threshold must be between 0-1, instead given {args.threshold}! Abort!")
                sys.exit(-1)
            valid_file(args.random_forest,rescue_logger)