import os, sys, json
import gzip

from src.commands import GFFREAD_PROG, run_command
from src.helpers import find_isoform_id_collisions
from src.module_logging import qc_logger, filter_logger, rescue_logger
from src.parsers import read_counts_design, check_counts_design, count_file_samples
from src.utilities.rescue.evidence_check import target_count_rules, check_design_for_rules

def valid_file(filename,logger):
    if not os.path.isfile(filename):
        logger.error(f"File {filename} not found. Abort!")
        sys.exit(1)
    return filename

def valid_fasta(filename,logger,allow_gz=False):
    """
    :param allow_gz: accept gzipped files (only for inputs read with the gzip-aware helpers)
    """
    valid_file(filename,logger)
    extension = filename.removesuffix('.gz').split('.')[-1]
    valid_extensions = ['fasta', 'fa', 'fastq', 'fq', 'faa', 'fna']
    if extension not in valid_extensions:
        logger.error(f"File {filename} is not a FASTA file. Abort!")
        sys.exit(1)
    if filename.endswith('.gz') and not allow_gz:
        logger.error(f"File {filename} is gzipped, but SQANTI3 needs this file uncompressed. Please decompress it. Abort!")
        sys.exit(1)
    return filename

def valid_unique_ids(filename,logger):
    """
    SQANTI3 shortens sequence IDs to the part before '|' or the first space. Abort if two sequences
    end up with the same ID, since they would become duplicated isoforms.
    """
    collisions = find_isoform_id_collisions(filename)
    if collisions:
        examples = "\n".join(f"  {new_id} <- {', '.join(ids)}" for new_id, ids in list(collisions.items())[:10])
        logger.error(f"{len(collisions)} sequence IDs in {filename} are not unique once shortened to the part "
                     f"before '|' or the first space (SQANTI3 needs unique IDs). "
                     f"{'First 10 cases' if len(collisions) > 10 else 'Cases'}:\n{examples}\n"
                     f"Please make the IDs unique before running SQANTI3. Abort!")
        sys.exit(1)
    return filename

def valid_gtf(filename,logger,out_dir):
    """
    Accepts a GTF, or a GFF/GFF3 that is converted to GTF with gffread.
    The converted file is written to out_dir, because the input folder may be read-only (issue #611).

    :return: path of the GTF to use downstream (the converted one for GFF input)
    """
    valid_file(filename,logger)
    converted = False
    if not filename.endswith('.gtf'):
        if not filename.endswith('.gff') and not filename.endswith('.gff3'):
            logger.error(f"File {filename} is not a GTF file. Abort!")
            sys.exit(1)
        else:
            logger.warning("GTF file is in GFF3 format. Converting to GTF format.")
            gtf_name = os.path.join(out_dir, os.path.splitext(os.path.basename(filename))[0] + ".converted.gtf")
            run_command(f"{GFFREAD_PROG} {filename} -T -o {gtf_name}", logger, f"{out_dir}/logs/gffread.log", "GFF to GTF conversion")
            logger.info(f"GFF file converted to GTF format. New file: {gtf_name}")
            filename = gtf_name
            converted = True

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
            # gffread exits without error on a file that is not GFF, and just writes an empty GTF
            if converted:
                logger.error(f"No annotation lines could be read from {filename}. Please check that the input is a valid GFF file. Abort!")
                sys.exit(1)
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

def valid_abundance_file(filename,logger):
    valid_file(filename,logger)
    if not filename.endswith('.tsv') and not filename.endswith('.csv'):
        logger.error(f"File {filename} is not a TSV/CSV file. Abort!")
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

def index_sr_bams(SR_bam, index_dir, cpus, logger):
    """
    The TSS ratio counts reads through the BAM index. The --SR_bam files are checked and indexed here,
    before the run is split, so that the chunks of a -n > 1 run share one index.
    """
    from src.utils import get_files_from_dir
    from src.utilities.short_reads import get_bam_index
    bams = [SR_bam] if SR_bam.endswith('.bam') else get_files_from_dir(SR_bam, ".bam")
    for bam in bams:
        valid_file(bam,logger)
        get_bam_index(bam, index_dir, cpus)

def contains_gene_name(filepath):
    open_function = gzip.open if filepath.endswith('.gz') else open
    with open_function(filepath, 'rt') as f:
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
    elif valid_fasta(args.isoforms,qc_logger,allow_gz=True):
        args.fasta = True
    else:
        qc_logger.error("Input isoforms must be in GTF, FASTA, or FASTQ format. Abort!")
        sys.exit(1)
    if args.fasta:
        valid_unique_ids(args.isoforms, qc_logger)

    valid_fasta(args.refFasta,qc_logger)

    # Customization validation
    if args.short_reads is not None:
        valid_file(args.short_reads,qc_logger)
    # BAM indexes made by SQANTI3 for the TSS ratio go here
    args.bam_index_dir = os.path.join(os.path.abspath(args.dir), "bam_index")
    if args.SR_bam is not None:
        valid_sr(args.SR_bam,qc_logger)
        index_sr_bams(args.SR_bam, args.bam_index_dir, args.cpus, qc_logger)

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
    # After valid_dir: a GFF reference is converted into the output directory
    args.refGTF_input = args.refGTF # the parameters file records what the user gave, not the converted GTF
    args.refGTF = valid_gtf(args.refGTF,qc_logger,args.dir)
    
    # FL counts given
    if args.fl_count is not None:
        valid_file(args.fl_count,qc_logger)
    if args.counts_design is not None:
        if args.fl_count is None:
            qc_logger.error("--counts_design requires a multi-sample --fl_count file.")
            sys.exit(1)
        valid_file(args.counts_design,qc_logger)
        # Checked against the samples of the count file here, before the long steps run,
        # and not again when the counts are read
        check_counts_design(read_counts_design(args.counts_design), count_file_samples(args.fl_count))
    if args.min_expression < 0:
        qc_logger.error(f"--min_expression must be 0 or greater, instead given {args.min_expression}.")
        sys.exit(1)

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
    # FL count file check
    if args.fl_count is not None:
        valid_abundance_file(args.fl_count,qc_logger)
    # Output prefix checks
    if args.output is None:
        args.output = os.path.splitext(os.path.basename(args.isoforms))[0]
    
    return args

def filter_args_validation(args):
    # Mandatory + possible inputs
    valid_file(args.sqanti_class, filter_logger)
    if args.filter_isoforms is not None:
        valid_fasta(args.filter_isoforms, filter_logger)
    if args.filter_sam is not None:
        valid_file(args.filter_sam, filter_logger)
    if args.filter_faa is not None:
        valid_fasta(args.filter_faa, filter_logger)

    valid_dir(args.dir,filter_logger)
    # After valid_dir: a GFF is converted into the output directory
    if args.filter_gtf is not None:
        args.filter_gtf_input = args.filter_gtf # the parameters file records what the user gave, not the converted GTF
        args.filter_gtf = valid_gtf(args.filter_gtf, filter_logger, args.dir)
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
    """
    Validate the rescue inputs before any step runs.

    :return: whether the evidence check of the rescued reference transcripts has to be run.
    """
    valid_file(args.filter_class, rescue_logger)
    valid_dir(args.dir,rescue_logger)
    args.refGTF_input = args.refGTF # the parameters file records what the user gave, not the converted GTF
    args.refGTF = valid_gtf(args.refGTF,rescue_logger,args.dir)
    valid_fasta(args.refFasta,rescue_logger)
    if args.filtered_isoforms_gtf is not None:
        args.filtered_isoforms_gtf_input = args.filtered_isoforms_gtf
        args.filtered_isoforms_gtf = valid_gtf(args.filtered_isoforms_gtf,rescue_logger,args.dir)
    if args.mode == "full":
        try:
            valid_file(args.refClassif,rescue_logger)
        except:
            rescue_logger.error("When running the rescue in full mode, the reference classification file is mandatory.")
            sys.exit(1)
        try:
            valid_fasta(args.corrected_isoforms_fasta,rescue_logger)
        except:
            rescue_logger.error("When running the rescue in full mode, the corrected isoforms FASTA file is mandatory.")
            sys.exit(1)
        if args.strategy == 'rules':
            valid_file(args.json_filter, rescue_logger)
        if args.strategy == "ml":
            if args.threshold < 0 or args.threshold > 1.:
                rescue_logger.error(f"--threshold must be between 0-1, instead given {args.threshold}! Abort!")
                sys.exit(-1)
            valid_file(args.random_forest,rescue_logger)
    if args.counts_design is not None:
        valid_file(args.counts_design, rescue_logger)
    if args.map_automatic_fsm and args.corrected_isoforms_fasta is None:
        rescue_logger.error("--map_automatic_fsm needs the corrected isoforms FASTA file (--corrected_isoforms_fasta).")
        sys.exit(1)
    if args.min_expression < 0:
        rescue_logger.error(f"--min_expression must be 0 or greater, instead given {args.min_expression}.")
        sys.exit(1)
    if args.counts is not None:
        valid_file(args.counts, rescue_logger)
    elif not args.skip_requant:
        rescue_logger.error("Requantification is run by default and needs the --counts file. "
                            "Provide --counts, or --skip_requant to rescue without requantification.")
        sys.exit(1)
    # The design is checked against the samples of --counts here, before mapping, and not again later
    design = None
    if args.counts_design is not None and args.counts is not None:
        design = check_counts_design(read_counts_design(args.counts_design, rescue_logger),
                                     count_file_samples(args.counts), rescue_logger)
    # Evidence check of the rescued reference transcripts: only for the rules strategy, when the
    # rules of reference targets have requisites on long-read counts (FL, prevalence, prevalence_<group>).
    # It also reads --counts, so the file is needed even with --skip_requant.
    evidence_check = False
    if args.skip_evidence_check:
        rescue_logger.info("Evidence check of reference targets skipped (--skip_evidence_check).")
    elif args.strategy == "rules":
        valid_file(args.json_filter, rescue_logger)
        with open(args.json_filter) as f:
            rule_sets = target_count_rules(json.load(f))
        evidence_check = bool(rule_sets)
        if evidence_check and args.counts is None:
            rescue_logger.error("The rules have requisites on long-read counts: the evidence check of the "
                                "rescued reference transcripts needs the --counts file used for requantification. "
                                "Provide --counts, or --skip_evidence_check to reintroduce them without the check.")
            sys.exit(1)
        if evidence_check:
            check_design_for_rules(rule_sets, design)
    return evidence_check
