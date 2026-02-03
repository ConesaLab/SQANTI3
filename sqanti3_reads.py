#!/usr/bin/env python3
import subprocess, os, sys, glob
import pandas as pd
import hashlib
import shutil

# SQANTI_Reads: Structural and Quality Annotation of Novel Transcripts in reads
# Author: Carolina Monzo

from src.reads_argparse import reads_argparser
from src.module_logging import reads_logger, update_logger
import logging

__author__  = "carolina.monzo@csic.es"
__version__ = '1.0'  # Python 3.7


utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "src/utilities")
sys.path.insert(0, utilitiesPath)
sqantiqcPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))

def fill_design_table(args):
    df = pd.read_csv(args.inDESIGN, sep = ",")
    # If number of columns is less than 2, probably wrongly formatted
    if df.shape[1] < 2:
        reads_logger.error(f"ERROR: {args.inDESIGN} is incorrectly formatted, is it not separated by commas?")
        sys.exit(-1)

    # Create the new columns
    # We always overwrite these columns to ensure they match the current run arguments
    # regardless of what might be in the input CSV (which could have stale paths).
    
    # Classification file is the OUTPUT file we will create
    df['classification_file'] = args.OUTPUT + '/' + df['file_acc'] + '/' + df['sampleID'] + '_reads_classification.txt'
    
    # Junction file is the INPUT file we need to read
    # If using sqanti_dirs (fast mode), it's there. If running from scratch, it will be in OUTPUT.
    # We default to pointing to sqanti_dirs here; if running from scratch, get_method_runSQANTI3 logic 
    # implies we'll find it eventually or create it.
    if args.sqanti_dirs:
         df['junction_file'] = args.sqanti_dirs + '/' + df['file_acc'] + '/' + df['sampleID'] + '_junctions.txt'
    else:
         # If no sqanti_dirs provided, we assume we are generating it in OUTPUT
         df['junction_file'] = args.OUTPUT + '/' + df['file_acc'] + '/' + df['sampleID'] + '_junctions.txt'
        
    return(df)

def get_method_runSQANTI3(args, df):

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        
        # Only check for existing SQANTI3 output if sqanti_dirs is provided
        if args.sqanti_dirs:
            classification_file = os.path.join(args.sqanti_dirs, file_acc, f"{sampleID}_classification.txt")
            junction_file = row["junction_file"]

            # Check for directory containing classification and junction file
            directory_path = os.path.join(args.sqanti_dirs, file_acc)

            if args.verbose:
                reads_logger.debug(f"Checking for directory: {directory_path}")
                reads_logger.debug(f"Checking for classification file: {classification_file}")
                reads_logger.debug(f"Checking for junction file: {junction_file}")

            if os.path.isdir(directory_path):
                classification_file_path = os.path.join(classification_file)
                junction_file_path = os.path.join(junction_file)
                
                if os.path.isfile(classification_file_path) and os.path.isfile(junction_file_path):
                    if args.verbose:
                        reads_logger.debug(f"[INFO] You inputted SQANTI3 directories, we will run sqanti_reads in fast mode for sample {directory_path}")
                    continue
                else:
                     if args.verbose:
                         reads_logger.debug(f"Directory found but files missing:\n  Class: {os.path.isfile(classification_file_path)}\n  Junc: {os.path.isfile(junction_file_path)}")
        
        # Check for .gtf or .gff file
        gtf_pattern = os.path.join(args.input_dir, f"{file_acc}*.g*f")
        reads_logger.debug(gtf_pattern)
        try:
            gtf_files = glob.glob(gtf_pattern)[0]
        except IndexError:
            pass
        else:
            if os.path.isfile(gtf_files):
                if os.path.isfile(args.refFasta) is False:
                    reads_logger.error(f'[ERROR] You inputted gtf files to run SQANTI3 but no reference genome FASTA')
                    sys.exit(-1)
                if os.path.isfile(args.refGTF) is False:
                    reads_logger.error(f'[ERROR] You inputted gtf files to run SQANTI3 but no reference annotation GTF')
                    sys.exit(-1)
                if args.verbose:
                    reads_logger.debug(f'[INFO] You inputted gtfs, we will run sqanti_reads in simple mode for sample {gtf_files}')
                
                # Update output directory to be args.OUTPUT instead of args.dir
                cmd_sqanti = (
                    f"python {sqantiqcPath}/sqanti3_qc.py "
                    f"--isoforms {gtf_files} "
                    f"--refGTF {args.refGTF} "
                    f"--refFasta {args.refFasta} "
                    f"--min_ref_len {args.min_ref_len} "
                    f"--aligner_choice {args.aligner_choice} "
                    f"-t {args.cpus} "
                    f"-d {args.OUTPUT}/{file_acc} "
                    f"-o {sampleID} "
                    f"-s {args.sites}"
                )

                subprocess.call(cmd_sqanti, shell = True)
                continue

        # Check for .fastq files
        fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fastq")
        try:
            fastq_files = glob.glob(fastq_pattern)[0]
        except IndexError:
            pass
        else:
            if os.path.isfile(fastq_files):
                if os.path.isfile(args.refFasta) is False:
                    reads_logger.error(f'[ERROR] You inputted fastq files to map but no reference genome FASTA')
                    sys.exit(-1)
                if os.path.isfile(args.refGTF) is False:
                    reads_logger.error(f'[ERROR] You inputted fastq files to map but no reference annotation GTF')
                    sys.exit(-1)
                if args.verbose:
                    reads_logger.debug(f'[INFO] You inputted reads, we will run sqanti_reads in simple mode for sample {fastq_files}')

                # Update output directory to be args.OUTPUT instead of args.dir
                cmd_sqanti = (
                    f"python {sqantiqcPath}/sqanti3_qc.py "
                    f"--isoforms {fastq_files} "
                    f"--refGTF {args.refGTF} "
                    f"--refFasta {args.refFasta} "
                    f"--min_ref_len {args.min_ref_len} "
                    f"--aligner_choice {args.aligner_choice} "
                    f"-t {args.cpus} "
                    f"-d {args.OUTPUT}/{file_acc} "
                    f"-o {sampleID} "
                    f"-s {args.sites} "
                    f"-n {args.chunks} "
                    f"--fasta"
                )

                reads_logger.debug(cmd_sqanti)
                subprocess.call(cmd_sqanti, shell = True)
                continue

        # If none of the conditions are met, raise an error
        reads_logger.error(f"ERROR: The file_acc you included in your design file ({file_acc}) does not correspond to .fastq, .gtf or directories with junctions and classification files in the {args.sqanti_dirs} or {args.input_dir} directory")
        sys.exit(-1)

def make_UJC_hash(args, df):

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        
        # Check if we processed this sample just now (in get_method_runSQANTI3)
        # If so, the output is already in args.OUTPUT.
        # If not (fast mode), the output is in args.dir.
        possible_input_path_new = os.path.join(args.OUTPUT, file_acc, sampleID)
        
        # We check for the classification file to confirm presence
        if os.path.exists(f"{possible_input_path_new}_classification.txt"):
             inputPathPrefix = possible_input_path_new
        elif args.sqanti_dirs and os.path.exists(os.path.join(args.sqanti_dirs, file_acc, f"{sampleID}_classification.txt")):
             inputPathPrefix = os.path.join(args.sqanti_dirs, file_acc, sampleID)
        else:
             reads_logger.error(f"Could not find SQANTI3 output for sample {sampleID} in {args.OUTPUT} or {args.sqanti_dirs}")
             sys.exit(-1)
        
        # Output: Always args.OUTPUT (where we write modified files)
        outputPathPrefix = os.path.join(args.OUTPUT, file_acc, sampleID)
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(outputPathPrefix), exist_ok=True)

<<<<<<< HEAD
        reads_logger.info("**** Calculating UJCs...")
                
=======
        print("**** Calculating UJCs...", file = sys.stdout)

>>>>>>> master
        # Ensure the corrected GTF contains gene_id attributes on every exon/CDS line so that
        # downstream `gtftools` does not fail with `IndexError: list index out of range`.
        # `gffread` will rewrite the file adding the missing attributes. We write to a
        # temporary file and then move it back so the final filename remains unchanged.
        
        # We need to work on a copy of the GTF to avoid modifying input
        input_gtf = f"{inputPathPrefix}_corrected.gtf"
        output_gtf = f"{outputPathPrefix}_corrected.gtf"
        
        if not os.path.exists(output_gtf):
             if os.path.exists(input_gtf):
                 shutil.copy(input_gtf, output_gtf)
             else:
                 reads_logger.error(f"Input GTF not found: {input_gtf}")
                 sys.exit(-1)

        gffread_tmp = f"{outputPathPrefix}_corrected.gtf.gffread_tmp"
        gffread_cmd = f"gffread {output_gtf} -T -o {gffread_tmp} && mv {gffread_tmp} {output_gtf}"

        try:
            subprocess.check_call(gffread_cmd, shell=True)
        except subprocess.CalledProcessError:
            reads_logger.error(f"ERROR running command: {gffread_cmd}\n Missing or failed gffread")
            sys.exit(-1)

        ## Take the corrected GTF
        # TODO: Change the file naming to use the standards of SQANTI3 via the helper function
        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {output_gtf} | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {output_gtf}"""
        ujc_cmd = f"""awk -F'\t' -v OFS="\t" '{{print $5,"chr"$1,$4,$2+1"_"$3}}' {outputPathPrefix}tmp_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{{print $1,$2"_"$3"_"$4}}' > {outputPathPrefix}tmp_UJC.txt"""

        if subprocess.check_call(introns_cmd, shell=True)!=0:
            reads_logger.error(f"ERROR running command: {introns_cmd}\n Missing GTFTOOLS")
            sys.exit(-1)

        if os.path.exists(f"{outputPathPrefix}_corrected.gtf.ensembl"):
            os.remove(f"{outputPathPrefix}_corrected.gtf.ensembl")

        if subprocess.check_call(ujc_cmd, shell=True)!=0:
            reads_logger.error(f"ERROR running command: {introns_cmd}\n Missing BEDTOOLS")
            sys.exit(-1)
        os.remove(f"{outputPathPrefix}tmp_introns.bed")

        ## Pandas merge to the left
        input_classfile = f"{inputPathPrefix}_classification.txt"
        if not os.path.exists(input_classfile):
             reads_logger.error(f"Input classification file not found: {input_classfile}")
             sys.exit(-1)
             
        clas_df = pd.read_csv(input_classfile, sep = "\t", usecols = [0, 1, 2, 7], dtype = "str")
        clas_df.columns = ["isoform", "chr", "strand", "associated_transcript"]
        ujc_df = pd.read_csv(f"{outputPathPrefix}tmp_UJC.txt", sep = "\t", names = ["isoform", "jxn_string"], dtype = "str")

        merged_df = pd.merge(clas_df, ujc_df, on = "isoform", how = "left")
        # Fill missing values in UJC column using the transcript ID
        merged_df["jxn_string"] = merged_df.apply(lambda row: row["chr"] + "_" + row["strand"] + "_" + "monoexon" + "_" + row["associated_transcript"] if pd.isna(row["jxn_string"]) else row["jxn_string"], axis=1)

        merged_df['jxnHash'] = merged_df['jxn_string'].apply(
                    lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest())

        merged_df.to_csv(f"{outputPathPrefix}_temp.txt", index = False, sep = "\t")
<<<<<<< HEAD
        
        cmd_paste = f"""bash -c 'paste <(cat {input_classfile} | tr -d '\r') <(cut -f 5,6 {outputPathPrefix}_temp.txt | tr -d '\r') > {outputPathPrefix}_reads_classification.txt'"""
=======

        cmd_paste = f"""bash -c 'paste <(cat {classfile} | tr -d '\r') <(cut -f 5,6 {outputPathPrefix}_temp.txt | tr -d '\r') > {outputPathPrefix}_reads_classification.txt'"""
>>>>>>> master
        subprocess.call(cmd_paste, shell = True)

        os.remove(f"{outputPathPrefix}tmp_UJC.txt")
        os.remove(f"{outputPathPrefix}_temp.txt")

def main():
    global utilitiesPath
    global sqantiqcPath

<<<<<<< HEAD
    args = reads_argparser().parse_args()
=======
    #arguments
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    parser.add_argument('--genome', type=str, help='\t\tReference genome (Fasta format).', default = False, required = False)
    parser.add_argument('--annotation', type=str, help='\t\tReference annotation file (GTF format).', default = False, required = True)
    parser.add_argument('-de', '--design', type=str, dest="inDESIGN" ,required=True, help='Path to design file, must have sampleID and file_acc column.')
    parser.add_argument('-i', '--input_dir', type=str, default = './', help = '\t\tPath to directory where fastq/GTF files are stored. Or path to parent directory with children directories of SQANTI3 runs. Default: Directory where the script was run.')
    parser.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='This is the column name that plots are to be faceted by. Default: None')
    parser.add_argument('-p','--prefix', type=str, dest="PREFIX", required=False, help='SQANTI-reads output filename prefix. Default: sqantiReads')
    parser.add_argument('-d','--dir', type=str, help='\t\tDirectory for output sqanti_reads files. Default: Directory where the script was run.', default = "./", required=False)
    parser.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length. Default: 0 bp")
    parser.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    parser.add_argument('--aligner_choice', type=str, choices=['minimap2', "uLTRA"], default='minimap2', help="\t\tDefault: minimap2")
    parser.add_argument('-t', '--cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. Default: 10')
    parser.add_argument('-n', '--chunks', default=1, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up. Default: 1')
    parser.add_argument('-s','--sites', type=str, default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    parser.add_argument('-ge','--gene_expression', type=int, dest="ANNOTEXP", required=False, help='Expression cut off level for determining underannotated genes. Default = 100', default = 100)
    parser.add_argument('-je','--jxn_expression', type=int, dest="JXNEXP", required=False, help='Coverage threshold for detected reference donors and acceptor. Default = 10', default = 10)
    parser.add_argument('-pc','--perc_coverage', type=int, dest="PERCCOV", required=False, help='Percent gene coverage of UJC for determining well-covered unannotated transcripts. Default = 20', default = 20)
    parser.add_argument('-pj','--perc_junctions', type=int, dest="PERCMAXJXN", required=False, help='Percent of the max junctions in gene for determining near full-length putative novel transcripts. Default = 80', default = 80)
    parser.add_argument('-fl','--factor_level', type=str, dest="FACTORLVL", required=False, help='Factor level to evaluate for underannotation', default = None)
    parser.add_argument('--all_tables', dest="ALLTABLES", action='store_true', help='Export all output tables. Default tables are gene counts, ujc counts, length_summary, cv and and underannotated gene tables')
    parser.add_argument('--pca_tables', dest="PCATABLES", action='store_true', help='Export table for making PCA plots')
    parser.add_argument('--skip_hash', dest="SKIPHASH", action='store_true', help='Skip the hashing step')
    parser.add_argument('--skip_plots', dest="SKIPPLOTS", action='store_true', help='Skip the plotting step')
    parser.add_argument('--report', type=str, choices = ["pdf", "html", "both"], default = 'pdf', help = "\t\tDefault: pdf")
    parser.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true")
    parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-reads '+str(__version__))
>>>>>>> master

    # Expand user paths (handle ~/) and absolute paths
    if args.sqanti_dirs:
        args.sqanti_dirs = os.path.abspath(os.path.expanduser(args.sqanti_dirs))
    if args.input_dir:
        args.input_dir = os.path.abspath(os.path.expanduser(args.input_dir))
    if args.OUTPUT:
        args.OUTPUT = os.path.abspath(os.path.expanduser(args.OUTPUT))
    if args.inDESIGN:
        args.inDESIGN = os.path.abspath(os.path.expanduser(args.inDESIGN))
    if args.refGTF:
        args.refGTF = os.path.abspath(os.path.expanduser(args.refGTF))
    if args.refFasta:
        args.refFasta = os.path.abspath(os.path.expanduser(args.refFasta))

    # Ensure output directory exists
    if not os.path.exists(args.OUTPUT):
        os.makedirs(args.OUTPUT)

    # Set up logger
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    
    update_logger(reads_logger, args.OUTPUT, "reads", log_level)

    # Check and read design file
    df = fill_design_table(args)
    
    # Save the processed design file to OUTPUT directory
    design_basename = os.path.basename(args.inDESIGN)
    new_design_path = os.path.join(args.OUTPUT, f"processed_{design_basename}")
    df.to_csv(new_design_path, index=False)
    
    # Update args to point to the new design file for downstream tools
    args.inDESIGN = new_design_path

    # Check method and run SQANTI3
    get_method_runSQANTI3(args, df)

    # Make UJC and hash
    if not args.SKIPHASH:
        make_UJC_hash(args, df)

<<<<<<< HEAD
    # Run plotting script directly as a function call
    reads_logger.info("Running SQANTI-reads tables and plots generation...")
    
    prefix = args.PREFIX if args.PREFIX else "sqantiReads"
    
    from src.utilities.sqanti_reads_tables_and_plots_02ndk import run_reads_plots
    run_reads_plots(
        ref_gtf=args.refGTF,
        design_file=args.inDESIGN,
        out_dir=args.OUTPUT,
        prefix=prefix,
        factor=args.inFACTOR,
        gene_expression=args.ANNOTEXP,
        jxn_expression=args.JXNEXP,
        perc_coverage=args.PERCCOV,
        perc_junctions=args.PERCMAXJXN,
        factor_level=args.FACTORLVL,
        all_tables=args.ALLTABLES,
        pca_tables=args.PCATABLES,
        report=args.report
    )
=======
    # Run plotting script
    if not args.SKIPPLOTS:
        plotting_script_path = os.path.join(os.path.dirname(__file__), 'src/utilities', 'sqanti_reads_tables_and_plots_02ndk.py')
        print(__file__)
        cmd_plotting = f"python {plotting_script_path} --ref {args.annotation} --design {args.inDESIGN} -o {args.dir} --gene-expression {args.ANNOTEXP} --jxn-expression {args.JXNEXP} --perc-coverage {args.PERCCOV} --perc-junctions {args.PERCMAXJXN} --report {args.report}"
        if args.inFACTOR:
            cmd_plotting = cmd_plotting + f" --factor {args.inFACTOR}"
        if args.FACTORLVL != None:
            cmd_plotting = cmd_plotting + f" --factor-level {args.FACTORLVL}"
        if args.PREFIX:
            cmd_plotting = cmd_plotting + f" --prefix {args.PREFIX}"
        else:
            cmd_plotting = cmd_plotting + " --prefix sqantiReads"
        if args.ALLTABLES:
            cmd_plotting = cmd_plotting + " --all-tables"
        if args.PCATABLES:
            cmd_plotting = cmd_plotting + " --pca-tables"
        print(cmd_plotting)

        subprocess.call(cmd_plotting, shell = True)
>>>>>>> master


if __name__ == "__main__":
    main()
