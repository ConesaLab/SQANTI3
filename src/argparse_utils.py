import os, sys
import subprocess
from src.commands import GFFREAD_PROG

def valid_file(filename):
    if not os.path.isfile(filename):
        print(f"ERROR: File {filename} not found. Abort!", file=sys.stderr)
        sys.exit(1)
    return filename

def valid_fasta(filename):
    valid_file(filename)
    if not filename.endswith('.fasta') and not filename.endswith('.fa') and not filename.endswith('.fastq') and not filename.endswith('.fq'):
        print(f"ERROR: File {filename} is not a FASTA file. Abort!", file=sys.stderr)
        sys.exit(1)
    return filename

def valid_gtf(filename):
    valid_file(filename)
    if not filename.endswith('.gtf'):
        if not filename.endswith('.gff'):
            print(f"ERROR: File {filename} is not a GTF file. Abort!", file=sys.stderr)
            sys.exit(1)
        else:
            print("WARNING: GTF file is in GFF3 format. Converting to GTF format...", file=sys.stderr)
            # GFF to GTF (in case the user provides gff instead of gtf)
            gtf_name = filename.replace('.gff', '.gtf')
            # Use the run command?
            try:
                subprocess.call([GFFREAD_PROG, filename , '-T', '-o', gtf_name])
            except (RuntimeError, TypeError, NameError):
                sys.stderr.write('ERROR: File %s without GTF/GFF format.\n' % filename)
                raise SystemExit(1)
            print("GFF file converted to GTF format. New file: {0}".format(gtf_name), file=sys.stderr)
            filename = gtf_name

    # Check if the GTF file is in the correct format
    ind = 0
    with open(filename) as isoforms_gtf:
        for line in isoforms_gtf:
            if line[0] != "#" and len(line.split("\t"))!=9:
                sys.stderr.write("\nERROR: input isoforms file with not GTF format.\n")
                sys.exit()
            elif len(line.split("\t"))==9:
                ind += 1
        if ind == 0:
            print("WARNING: GTF has {0} no annotation lines.".format(filename), file=sys.stderr)

    return filename

def valid_bed(filename):
    valid_file(filename)
    if not filename.endswith('.bed'):
        print(f"ERROR: File {filename} is not a BED file. Abort!", file=sys.stderr)
        sys.exit(1)
    return filename

def valid_dir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    else:
        print(f"WARNING: output directory {dirname} already exists. Overwriting!", file=sys.stderr)
    return dirname

# TODO: Add specific validation for Kallisto output
def valid_matrix(filename):
    valid_file(filename)
    if not filename.endswith('.tsv'):
        print(f"ERROR: File {filename} is not a TSV file. Abort!", file=sys.stderr)
        sys.exit(1)

#TODO: Get a condition to see if it is a pacbio file
def valid_PacBio_abund(filename):
    valid_file(filename)
    if not filename.endswith('abundance.tsv'):
        print(f"ERROR: File {filename} is not a PacBio abundance file. Abort!", file=sys.stderr)
        sys.exit(1)
    return filename

def valid_gff3(filename):
    valid_file(filename)
    if not filename.endswith('.gff3'):
        print(f"ERROR: File {filename} is not a GFF3 file. Abort!", file=sys.stderr)
        sys.exit(1)
    return filename

def valid_sr(filename):
    if not os.path.isdir(filename):
        if not filename.endswith('.fofn') and not filename.endswith('.bam'):
            print(f"ERROR: File {filename} is not a BAM file, a directory, or a FOFN file. Abort!", file=sys.stderr)
            sys.exit(1)
        else:
            if not os.path.isfile(filename):
                print(f"ERROR: File {filename} not found. Abort!", file=sys.stderr)
                sys.exit(1)
    return filename

### Validation for the arguments

def args_validation(args):
    # Required arguments
    valid_file(args.isoforms)
    if args.isoforms.endswith('.gtf') or args.isoforms.endswith('.gff'):
        args.fasta = False
    elif valid_fasta(args.isoforms):
        args.fasta = True
    else:
        print("ERROR: Input isoforms must be in GTF, FASTA, or FASTQ format. Abort!", file=sys.stderr)
        sys.exit(1)

    # Mappers checks
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
    
    # Fusion isoforms checks
    if args.is_fusion:
        if args.orf_input is None:
            print("WARNING: Currently if --is_fusion is used, no ORFs will be predicted. Supply --orf_input if you want ORF to run!", file=sys.stderr)
            args.skipORF = True
        if args.fasta:
            print("ERROR: if --is_fusion is on, must supply GTF as input", file=sys.stderr)
            sys.exit(1)
    # Expression checks
    if args.expression is not None:
        if os.path.isdir(args.expression)==True:
            print("Expression files located in {0} folder".format(args.expression), file=sys.stderr)
        else:
            for f in args.expression.split(','):
                if not os.path.exists(f):
                        print("Expression file {0} not found. Abort!".format(f), file=sys.stderr)
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

