import os, sys
import subprocess

from .commands import GFFREAD_PROG
from .logging_config import logger

def valid_file(filename):
    if not os.path.isfile(filename):
        logging.error(f"File {filename} not found. Abort!")
        sys.exit(1)
    return filename

def valid_fasta(filename):
    valid_file(filename)
    if not filename.endswith('.fasta') and not filename.endswith('.fa'):
        logging.error(f"File {filename} is not a FASTA file. Abort!")
        sys.exit(1)
    return filename

def valid_gtf(filename):
    valid_file(filename)
    if not filename.endswith('.gtf'):
        if not filename.endswith('.gff'):
            logging.error(f"File {filename} is not a GTF file. Abort!")
            sys.exit(1)
        else:
            logging.warning("GTF file is in GFF3 format. Converting to GTF format.")
            # GFF to GTF (in case the user provides gff instead of gtf)
            gtf_name = filename.replace('.gff', '.gtf')
            # Use the run command?
            try:
                subprocess.call([GFFREAD_PROG, filename , '-T', '-o', gtf_name])
            except (RuntimeError, TypeError, NameError):
                logging.error(f'File {filename} without GTF/GFF format.')
                raise SystemExit(1)
            logging.info(f"GFF file converted to GTF format. New file: {gtf_name}")
            filename = gtf_name

    # Check if the GTF file is in the correct format
    ind = 0
    with open(filename) as isoforms_gtf:
        for line in isoforms_gtf:
            if line[0] != "#" and len(line.split("\t"))!=9:
                logging.error("Input isoforms file with not GTF format.")
                sys.exit()
            elif len(line.split("\t"))==9:
                ind += 1
        if ind == 0:
            logging.warning(f"GTF has {filename} no annotation lines.")

    return filename

def valid_bed(filename):
    valid_file(filename)
    if not filename.endswith('.bed'):
        logging.error(f"File {filename} is not a BED file. Abort!")
        sys.exit(1)
    return filename

def valid_dir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    else:
        logging.warning(f"output directory {dirname} already exists. Overwriting!")
    return dirname

# TODO: Add specific validation for Kallisto output
def valid_matrix(filename):
    valid_file(filename)
    if not filename.endswith('.tsv'):
        logging.error(f"File {filename} is not a TSV file. Abort!")
        sys.exit(1)

#TODO: Get a condition to see if it is a pacbio file
def valid_PacBio_abund(filename):
    valid_file(filename)
    if not filename.endswith('abundance.tsv'):
        logging.error(f"File {filename} is not a PacBio abundance file. Abort!")
        sys.exit(1)
    return filename

def valid_gff3(filename):
    valid_file(filename)
    if not filename.endswith('.gff3'):
        logging.error(f"File {filename} is not a GFF3 file. Abort!")
        sys.exit(1)
    return filename

def valid_sr(filename):
    if not os.path.isdir(filename):
        if not filename.endswith('.fofn') and not filename.endswith('.bam'):
            logging.error(f"File {filename} is not a BAM file, a directory, or a FOFN file. Abort!")
            sys.exit(1)
        else:
            if not os.path.isfile(filename):
                logging.error(f"File {filename} not found. Abort!")
                sys.exit(1)
    return filename

### Validation for the arguments

def args_validation(args):
    # Required arguments
    if args.isoforms.endswith('.gtf'):
        args.fasta = False
    elif args.isoforms.endswith('.fasta') or args.isoforms.endswith('.fastq'):
        args.fasta = True
    else:
        logging.error("Input isoforms must be in GTF, FASTA, or FASTQ format. Abort!")
        sys.exit(1)

    # Mappers checks
    if args.fasta:
        if args.aligner_choice == 'gmap':
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                logging.error(f"GMAP index {args.gmap_index} doesn't exist! Abort.")
                sys.exit()
        elif args.aligner_choice == 'deSALT':
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                logging.error(f"deSALT index {args.gmap_index} doesn't exist! Abort.")
                sys.exit()

        logging.info("Cleaning up isoform IDs...")
        from src.helpers import rename_isoform_seqids
        args.isoforms = rename_isoform_seqids(args.isoforms, args.force_id_ignore)
        logging.info(f"Cleaned up isoform fasta file written to: {args.isoforms}")
    
    # Fusion isoforms checks
    if args.is_fusion:
        if args.orf_input is None:
            logging.warning("Currently if --is_fusion is used, no ORFs will be predicted. Supply --orf_input if you want ORF to run!")
            args.skipORF = True
        if args.fasta:
            logging.error("If --is_fusion is on, must supply GTF as input")
            sys.exit(1)
    # Expression checks
    if args.expression is not None:
        if os.path.isdir(args.expression)==True:
            logging.info(f"Expression files located in {args.expression} folder")
        else:
            for f in args.expression.split(','):
                if not os.path.exists(f):
                        logging.error(f"Expression file {f} not found. Abort!")
                        sys.exit(1)
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

