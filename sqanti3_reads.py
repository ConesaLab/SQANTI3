#!/usr/bin/env python3
import subprocess, os, re, sys, glob
import argparse
import pandas as pd
import shutil
import hashlib
#!/usr/bin/env python3
# SQANTI_Reads: Structural and Quality Annotation of Novel Transcripts in reads
# Author: Carolina Monzo

__author__  = "carolina.monzo@csic.es"
__version__ = '1.0'  # Python 3.7


utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "src/utilities")
sys.path.insert(0, utilitiesPath)
sqantiqcPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))


FIELDS_JUNC = ['isoform', 'chrom', 'strand', 'junction_number', 'genomic_start_coord',
                   'genomic_end_coord', 'junction_category',
                   'diff_to_Ref_start_site', 'diff_to_Ref_end_site', 'canonical']

FIELDS_CLASS = ['isoform', 'chrom', 'strand', 'length',  'exons',  'structural_category',
                'associated_gene', 'associated_transcript',  'ref_length', 'ref_exons',
                'subcategory', 'RTS_stage', 'all_canonical',
                'predicted_NMD', 'perc_A_downstream_TTS', "jxn_string"]

RSCRIPTPATH = shutil.which('Rscript')

def fill_design_table(args):
    df = pd.read_csv(args.inDESIGN, sep = ",")
    # If number of columns is less than 2, probably wrongly formatted
    if df.shape[1] < 2:
        print("ERROR: is incorrectly formatted, is it not separated by commas?".format(args.inDESIGN), file=sys.stderr)
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
                    print(f"[INFO] You inputted SQANTI3 directories, we will run sqanti_reads in fast mode for sample {directory_path}", file=sys.stdout)
                continue
        
        # Check for .gtf or .gff file
        gtf_pattern = os.path.join(args.input_dir, f"{file_acc}*.g*f")
        print(gtf_pattern)
        try:
            gtf_files = glob.glob(gtf_pattern)[0]
        except IndexError:
            pass
        else:
            if os.path.isfile(gtf_files):
                if os.path.isfile(args.genome) is False:
                    print(f'[ERROR] You inputted gtf files to run SQANTI3 but no reference genome FASTA', file=sys.stdout)
                    sys.exit(-1)
                if os.path.isfile(args.annotation) is False:
                    print(f'[ERROR] You inputted gtf files to run SQANTI3 but no reference annotation GTF', file=sys.stdout)
                    sys.exit(-1)
                if args.verbose:
                    print(f'[INFO] You inputted gtfs, we will run sqanti_reads in simple mode for sample {gtf_files}', file=sys.stdout)
                cmd_sqanti = f"python {sqantiqcPath}/sqanti3_qc.py --isoforms {gtf_files} --refGTF {args.annotation} --refFasta {args.genome} --skipORF --min_ref_len {args.min_ref_len} --aligner_choice {args.aligner_choice} -t {args.cpus} -d {args.dir}/{file_acc} -o {sampleID} -s {args.sites}"

                if args.force_id_ignore:
                    cmd_sqanti = cmd_sqanti + " --force_id_ignore"
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
                if os.path.isfile(args.genome) is False:
                    print(f'[ERROR] You inputted fastq files to map but no reference genome FASTA', file=sys.stdout)
                    sys.exit(-1)
                if os.path.isfile(args.annotation) is False:
                    print(f'[ERROR] You inputted fastq files to map but no reference annotation GTF', file=sys.stdout)
                    sys.exit(-1)
                if args.verbose:
                    print(f'[INFO] You inputted reads, we will run sqanti_reads in simple mode for sample {fastq_files}', file=sys.stdout)

                cmd_sqanti = f"python {sqantiqcPath}/sqanti3_qc.py \
                                {fastq_files} {args.annotation} {args.genome} \
                                --skipORF --min_ref_len {args.min_ref_len} \
                                --aligner_choice {args.aligner_choice} \
                                -t {args.cpus} -d {args.dir}/{file_acc} \
                                -o {sampleID} -s {args.sites} -n {args.chunks} \
                                --fasta"
                if args.force_id_ignore:
                    cmd_sqanti = cmd_sqanti + " --force_id_ignore"

                print(cmd_sqanti, file=sys.stdout)
                subprocess.call(cmd_sqanti, shell = True)
                continue
        
        # If none of the conditions are met, raise an error
        print(f"ERROR: The file_acc you included in your design file does not correspond to .fastq, .gtf or directories with junctions and classification files in the {args.input_dir} directory", file=sys.stdout)
        sys.exit(-1)

def make_UJC_hash(args, df):

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        # Input dir, sqanti3 dir, samplename
        outputPathPrefix = os.path.join(args.dir, file_acc, sampleID)

        print("**** Calculating UJCs...", file = sys.stdout)
                
        ## Take the corrected GTF
        # TODO: Change the file naming to use the standards of SQANTI3 via the helper function
        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {outputPathPrefix}_corrected.cds.gff3 | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {outputPathPrefix}_corrected.cds.gff3"""
        ujc_cmd = f"""awk -F'\t' -v OFS="\t" '{{print $5,"chr"$1,$4,$2+1"_"$3}}' {outputPathPrefix}tmp_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{{print $1,$2"_"$3"_"$4}}' > {outputPathPrefix}tmp_UJC.txt"""
            
        if subprocess.check_call(introns_cmd, shell=True)!=0:
            print("ERROR running command: {0}\n Missing GTFTOOLS".format(introns_cmd), file=sys.stderr)
            sys.exit(-1)
            
        if os.path.exists(f"{outputPathPrefix}_corrected.gtf.ensembl"):
            os.remove(f"{outputPathPrefix}_corrected.gtf.ensembl")
            
        if subprocess.check_call(ujc_cmd, shell=True)!=0:
            print("ERROR running command: {0}\n Missing BEDTOOLS".format(introns_cmd), file=sys.stderr)
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
    parser.add_argument('--report', type=str, choices = ["pdf", "html", "both"], default = 'pdf', help = "\t\tDefault: pdf")
    parser.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true")
    parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-reads '+str(__version__))

    args = parser.parse_args()

    # Check and read design file
    df = fill_design_table(args)

    # Check method and run SQANTI3
    get_method_runSQANTI3(args, df)

    # Make UJC and hash
    if not args.SKIPHASH:
        make_UJC_hash(args, df)

    # Run plotting script
    plotting_script_path = os.path.join(os.path.dirname(__file__), 'src/utilities', 'sqanti_reads_tables_and_plots_02ndk.py')

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
        cmd_plotting = cmd_plotting + "--pca-tables"
    print(cmd_plotting)

    subprocess.call(cmd_plotting, shell = True)


if __name__ == "__main__":
    main()

