#!usr/bin/python
import subprocess, os, re, sys
import argparse
import pandas as pd
import distutils.spawn
#!/usr/bin/env python3
# SQANTI_Reads: Structural and Quality Annotation of Novel Transcripts in reads
# Author: Carolina Monzo

__author__  = "carolina.monzo@csic.es"
__version__ = '5.2'  # Python 3.7


utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sqantiqcPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))


FIELDS_JUNC = ['read', 'chrom', 'strand', 'junction_number', 'genomic_start_coord',
                   'genomic_end_coord', 'junction_category',
                   'diff_to_Ref_start_site', 'diff_to_Ref_end_site', 'canonical']

FIELDS_CLASS = ['read', 'chrom', 'strand', 'length',  'exons',  'structural_category',
                'associated_gene', 'associated_transcript',  'ref_length', 'ref_exons',
                'subcategory', 'RTS_stage', 'all_canonical',
                'predicted_NMD', 'perc_A_downstream_TTS', "Jxn_string"]

RSCRIPTPATH = distutils.spawn.find_executable('Rscript')

def run_sqantiReads(args):
    outputPathPrefix = os.path.join(args.dir, args.output)
    
    if args.classification is not None:
        # Reduce the junctions file
        cmd_reduceJunctions = f"cat {outputPathPrefix}_junctions.txt | cut -f 1,2,3,4,5,6,8,11,12,15 > {outputPathPrefix}_r_junctions.txt"

        if subprocess.check_call(cmd_reduceJunctions, shell = True) !=0:
            print("Problem reducing the junction file")
            sys.exit(-1)
        os.remove(f"{outputPathPrefix}_junctions.txt")
        cmd_sed = f"sed '1s/isoform/read/' {outputPathPrefix}_r_junctions.txt > {outputPathPrefix}_reads_junctions.txt"
        subprocess.call(cmd_sed, shell = True)
        os.remove(f"{outputPathPrefix}_r_junctions.txt")


    print("**** Calculating UJCs...", file = sys.stdout)
            
    ## Take the corrected GTF
    introns_cmd = f"gtftools -i {outputPathPrefix}tmp_introns.bed {outputPathPrefix}_corrected.gtf"
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

    # Reduce the classification file
    if args.classification is not None:
        classfile = f"{outputPathPrefix}_r_classification.txt"
        cmd_reduceClassification = f"cat {outputPathPrefix}_classification.txt | cut -f 1,2,3,4,5,6,7,8,9,10,15,16,17,37,38 | sed '1s/isoform/read/' > {outputPathPrefix}_r_classification.txt"
        if subprocess.check_call(cmd_reduceClassification, shell = True) !=0:
            print("Problem reducing the classification file")
            sys.exit(-1)
    else:
        classfile = f"{outputPathPrefix}_classification.txt"
        

    ## Pandas merge to the left
    clas_df = pd.read_csv(classfile, sep = "\t", usecols = [0, 7])
    clas_df.columns = ["read", "associated_transcript"]
    ujc_df = pd.read_csv(f"{outputPathPrefix}tmp_UJC.txt", sep = "\t", names = ["read", "jxn_string"])
    
    merged_df = pd.merge(clas_df, ujc_df, on = "read", how = "left")
    
    # Fill missing values in UJC column using the transcript ID
    merged_df["jxn_string"] = merged_df.apply(lambda row: 'Mono_' + row["associated_transcript"] if pd.isna(row["jxn_string"]) else row["jxn_string"], axis=1)
    
    merged_df.to_csv(f"{outputPathPrefix}_temp.txt", index = False, sep = "\t")
    
    cmd_paste = ['bash -c "paste <(cat ', classfile, ') <(cut -f 3 ', outputPathPrefix, '_temp.txt) > ', outputPathPrefix, '_reads_classification.txt"']
    subprocess.call("".join(cmd_paste), shell = True)
    
    if args.classification is not None:
        os.remove(f"{outputPathPrefix}_r_classification.txt")
    
    os.remove(f"{outputPathPrefix}tmp_UJC.txt")
    os.remove(f"{outputPathPrefix}_temp.txt")
    os.remove(f"{outputPathPrefix}_classification.txt")





def main():
    global utilitiesPath
    global sqantiqcPath

    #arguments
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    parser.add_argument('--reads', default = None, help='\Reads (FASTA/FASTQ) or GTF format. It is recommended to provide them in GTF format, but if it is needed to map the sequences to the genome use a FASTA/FASTQ file with the --fasta option.')
    parser.add_argument('--annotation', help='\t\tReference annotation file (GTF format)', default = None)
    parser.add_argument('--genome', help='\t\tReference genome (Fasta format)', default = None)
    parser.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length (default: 0 bp)")
    parser.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    parser.add_argument('--aligner_choice', choices=['minimap2', "uLTRA"], default='minimap2')
    parser.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
    parser.add_argument('-t', '--cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. (default: 10)')
    parser.add_argument('-n', '--chunks', default=1, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up (default: 1).')
    parser.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
    parser.add_argument('-s','--sites', default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='SQANTI3 '+str(__version__))
    parser.add_argument('--saturation', action="store_true", default=False, help='\t\tInclude saturation curves into report')
    parser.add_argument('--fasta', help='\t\tUse when running SQANTI by using as input a FASTA/FASTQ with the sequences of isoforms', action='store_true')
    args = parser.parse_args()

    ## If they give us a classification file, don't run sqanti, run only sqanti_reads
    if args.reads is not None and args.annotation is not None and args.genome is not None:
        cmd_sqanti = f"python3 {sqantiqcPath}sqanti_qc.py {args.reads} {args.annotation} {args.genome} --min_ref_len {args.min_ref_len} --force_id_ignore {args.force_id_ignore} --aligner_choice {args.aligner_choice} -t {args.cpus} -d {args.dir} -o {args.output} -s {args.sites} -v {args.version} --saturation {args.saturation} --fasta {args.fasta}"
            
        subprocess.call(cmd_sqanti, shell = True)
        args.classification = "caca"
    else:
        args.classification = None
    
    ## Check that there is a classification and a corrected.gtf file in the output directory
    if os.path.exists(os.path.join(args.dir, f"{args.output}_classification.txt")) and os.path.exists(os.path.join(args.dir, f"{args.output}_corrected.gtf")):
        ### Run SQANTI READS ###
        run_sqantiReads(args)
        print("SQANTI-reads run complete")
    else:
        print(f"ERROR: There are no classification.txt and corrected.gtf files in your output directory: {args.dir}")
    
        

if __name__ == "__main__":
    main()

