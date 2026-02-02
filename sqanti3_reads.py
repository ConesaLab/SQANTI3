#!/usr/bin/env python3
import subprocess, os, sys, glob
import pandas as pd
import hashlib

# SQANTI_Reads: Structural and Quality Annotation of Novel Transcripts in reads
# Author: Carolina Monzo

from src.reads_argparse import reads_argparser
from src.module_logging import reads_logger, update_logger
from src.utilities.sqanti_reads_tables_and_plots_02ndk import run_reads_plots
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
    df['classification_file'] = args.dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_reads_classification.txt'
    df['junction_file'] = args.dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_junctions.txt'
    df.to_csv(args.inDESIGN, sep = ',', index = False)
    return(df)

def get_method_runSQANTI3(args, df):

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        classification_file = os.path.join(args.dir, file_acc, f"{sampleID}_classification.txt")
        junction_file = row["junction_file"]

        # Check for directory containing classification and junction file
        directory_path = os.path.join(args.dir, file_acc)

        if os.path.isdir(directory_path):
            classification_file_path = os.path.join(classification_file)
            junction_file_path = os.path.join(junction_file)
            if os.path.isfile(classification_file_path) and os.path.isfile(junction_file_path):
                if args.verbose:
                    reads_logger.debug(f"[INFO] You inputted SQANTI3 directories, we will run sqanti_reads in fast mode for sample {directory_path}")
                continue

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
                cmd_sqanti = (
                    f"python {sqantiqcPath}/sqanti3_qc.py "
                    f"--isoforms {gtf_files} "
                    f"--refGTF {args.refGTF} "
                    f"--refFasta {args.refFasta} "
                    f"--min_ref_len {args.min_ref_len} "
                    f"--aligner_choice {args.aligner_choice} "
                    f"-t {args.cpus} "
                    f"-d {args.dir}/{file_acc} "
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

                cmd_sqanti = (
                    f"python {sqantiqcPath}/sqanti3_qc.py "
                    f"--isoforms {fastq_files} "
                    f"--refGTF {args.refGTF} "
                    f"--refFasta {args.refFasta} "
                    f"--min_ref_len {args.min_ref_len} "
                    f"--aligner_choice {args.aligner_choice} "
                    f"-t {args.cpus} "
                    f"-d {args.dir}/{file_acc} "
                    f"-o {sampleID} "
                    f"-s {args.sites} "
                    f"-n {args.chunks} "
                    f"--fasta"
                )

                reads_logger.debug(cmd_sqanti)
                subprocess.call(cmd_sqanti, shell = True)
                continue

        # If none of the conditions are met, raise an error
        reads_logger.error(f"ERROR: The file_acc you included in your design file does not correspond to .fastq, .gtf or directories with junctions and classification files in the {args.input_dir} directory")
        sys.exit(-1)

def make_UJC_hash(args, df):

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        # Input dir, sqanti3 dir, samplename
        outputPathPrefix = os.path.join(args.dir, file_acc, sampleID)

        reads_logger.info("**** Calculating UJCs...")
                
        # Ensure the corrected GTF contains gene_id attributes on every exon/CDS line so that
        # downstream `gtftools` does not fail with `IndexError: list index out of range`.
        # `gffread` will rewrite the file adding the missing attributes. We write to a
        # temporary file and then move it back so the final filename remains unchanged.
        gffread_tmp = f"{outputPathPrefix}_corrected.gtf.gffread_tmp"
        gffread_cmd = f"gffread {outputPathPrefix}_corrected.gtf -T -o {gffread_tmp} && mv {gffread_tmp} {outputPathPrefix}_corrected.gtf"

        try:
            subprocess.check_call(gffread_cmd, shell=True)
        except subprocess.CalledProcessError:
            reads_logger.error(f"ERROR running command: {gffread_cmd}\n Missing or failed gffread")
            sys.exit(-1)

        ## Take the corrected GTF
        # TODO: Change the file naming to use the standards of SQANTI3 via the helper function
        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {outputPathPrefix}_corrected.gtf | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {outputPathPrefix}_corrected.gtf"""
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
        classfile = f"{outputPathPrefix}_classification.txt"
        clas_df = pd.read_csv(classfile, sep = "\t", usecols = [0, 1, 2, 7], dtype = "str")
        clas_df.columns = ["isoform", "chr", "strand", "associated_transcript"]
        ujc_df = pd.read_csv(f"{outputPathPrefix}tmp_UJC.txt", sep = "\t", names = ["isoform", "jxn_string"], dtype = "str")

        merged_df = pd.merge(clas_df, ujc_df, on = "isoform", how = "left")
        # Fill missing values in UJC column using the transcript ID
        merged_df["jxn_string"] = merged_df.apply(lambda row: row["chr"] + "_" + row["strand"] + "_" + "monoexon" + "_" + row["associated_transcript"] if pd.isna(row["jxn_string"]) else row["jxn_string"], axis=1)

        merged_df['jxnHash'] = merged_df['jxn_string'].apply(
                    lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest())

        merged_df.to_csv(f"{outputPathPrefix}_temp.txt", index = False, sep = "\t")

        cmd_paste = f"""bash -c 'paste <(cat {classfile} | tr -d '\r') <(cut -f 5,6 {outputPathPrefix}_temp.txt | tr -d '\r') > {outputPathPrefix}_reads_classification.txt'"""
        subprocess.call(cmd_paste, shell = True)

        os.remove(f"{outputPathPrefix}tmp_UJC.txt")
        os.remove(f"{outputPathPrefix}_temp.txt")

def main():
    global utilitiesPath
    global sqantiqcPath

    args = reads_argparser().parse_args()

    # Set up logger
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    
    update_logger(reads_logger, args.dir, "reads", log_level)

    # Check and read design file
    df = fill_design_table(args)

    # Check method and run SQANTI3
    get_method_runSQANTI3(args, df)

    # Make UJC and hash
    if not args.SKIPHASH:
        make_UJC_hash(args, df)

    # Run plotting script directly as a function call
    if not args.SKIPPLOTS:
        reads_logger.info("Running SQANTI-reads tables and plots generation...")
        
        prefix = args.PREFIX if args.PREFIX else "sqantiReads"
        
        run_reads_plots(
            ref_gtf=args.refGTF,
            design_file=args.inDESIGN,
            out_dir=args.dir,
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


if __name__ == "__main__":
    main()
