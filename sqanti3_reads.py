#!/usr/bin/env python
import subprocess, os, sys, glob
import pandas as pd
import hashlib

# SQANTI_Reads: Structural and Quality Annotation of Novel Transcripts in reads
# Author: Carolina Monzo

from src.reads_argparse import reads_argparser
from src.module_logging import reads_logger, update_logger
import logging

__author__  = "carolina.monzo@csic.es"
__version__ = '1.0'  # Python 3.7

# Streaming awk that emits each multi-exon transcript's unique junction chain (UJC)
# straight from the corrected GTF, replacing the old gffread -> gtftools -> bedtools
# pipeline (gffread scaled ~quadratically, ~11 h on a 6 GB GTF; gtftools then
# reparsed the whole file in pure Python, ~75 min/sample). A transcript's introns
# are just the gaps between its exons sorted by start, so the junction string
#   chr<chrom>_<strand>_<i1start>_<i1end>_<i2start>_<i2end>...
# (each intron = prev_exon_end+1 .. next_exon_start-1) is byte-identical to the old
# output. Exons are insertion-sorted per transcript (no dependency on gffread's
# normalisation); monoexons produce no line (jxn_string stays NaN, keyed separately).
_UJC_AWK = r'''BEGIN{ FS="\t" }
function emit(  i,j,ks,ke,chain){
  if(n<2){ n=0; return }
  for(i=2;i<=n;i++){ ks=es[i]; ke=ee[i]; j=i-1; while(j>=1 && es[j]>ks){ es[j+1]=es[j]; ee[j+1]=ee[j]; j-- } es[j+1]=ks; ee[j+1]=ke }
  chain=""; for(i=1;i<n;i++) chain=chain"_"(ee[i]+1)"_"(es[i+1]-1); print cur"\t""chr"c"_"st chain; n=0
}
$3=="exon"{
  p=index($0,"transcript_id \""); if(p==0) next; rest=substr($0,p+15); q=index(rest,"\""); tid=substr(rest,1,q-1);
  if(tid!=cur){ emit(); cur=tid; c=$1; sub(/^chr/,"",c); st=$7 }
  n++; es[n]=$4; ee[n]=$5
}
END{ emit() }'''


utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "src/utilities")
sys.path.insert(0, utilitiesPath)
sqantiqcPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))


def _run_parallel(func, items, jobs):
    """Run func(item) over items, in parallel with a thread pool when jobs>1.

    A thread pool is used (not processes) because the per-sample work is
    dominated by external subprocesses / file I/O, which release the GIL — so
    this gives real parallelism without pickling. Exceptions propagate.
    """
    jobs = max(1, int(jobs or 1))
    items = list(items)
    if jobs <= 1 or len(items) <= 1:
        for it in items:
            func(it)
        return
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=min(jobs, len(items))) as ex:
        # list() forces evaluation so any worker exception is re-raised here.
        list(ex.map(func, items))

def fill_design_table(args):
    df = pd.read_csv(args.inDESIGN, sep = ",")
    # If number of columns is less than 2, probably wrongly formatted
    if df.shape[1] < 2:
        reads_logger.error(f"ERROR: {args.inDESIGN} is incorrectly formatted, is it not separated by commas?")
        sys.exit(-1)

    # Create the new columns
    # We always overwrite these columns to ensure they match the current run arguments
    # regardless of what might be in the input CSV (which could have stale paths).
    
    # Classification file: normally the (jxnHash-augmented) file that the UJC
    # hashing step writes into OUTPUT. With --skip_hash that step is not run, so
    # the file is never created in OUTPUT; read the already-hashed classification
    # straight from the existing SQANTI3 dirs instead (otherwise the run fails
    # looking for a file that was never created).
    if getattr(args, "SKIPHASH", False) and args.sqanti_dirs:
        df['classification_file'] = args.sqanti_dirs + '/' + df['file_acc'] + '/' + df['sampleID'] + '_reads_classification.txt'
    else:
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

def _run_qc_cmd(cmd):
    subprocess.call(cmd, shell=True)


def get_method_runSQANTI3(args, df):

    qc_cmds = []
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
                if not args.refFasta or os.path.isfile(args.refFasta) is False:
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

                qc_cmds.append(cmd_sqanti)
                continue

        # Check for .fastq files
        fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fastq")
        try:
            fastq_files = glob.glob(fastq_pattern)[0]
        except IndexError:
            pass
        else:
            if os.path.isfile(fastq_files):
                if not args.refFasta or os.path.isfile(args.refFasta) is False:
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
                qc_cmds.append(cmd_sqanti)
                continue

        # If none of the conditions are met, raise an error
        reads_logger.error(f"ERROR: The file_acc you included in your design file ({file_acc}) does not correspond to .fastq, .gtf or directories with junctions and classification files in the {args.sqanti_dirs} or {args.input_dir} directory")
        sys.exit(-1)

    # Run the per-sample SQANTI3-QC commands (simple mode), in parallel if --jobs>1.
    if qc_cmds:
        reads_logger.info(f"Running SQANTI3-QC for {len(qc_cmds)} sample(s) with jobs={getattr(args, 'jobs', 1)}")
        _run_parallel(_run_qc_cmd, qc_cmds, getattr(args, 'jobs', 1))

def make_UJC_hash(args, df):

    # Per-sample UJC hashing (single awk pass over the corrected GTF + pandas). Each sample
    # writes only into its own output dir, so samples are independent and run in
    # parallel when --jobs>1 (thread pool: the work is subprocess/I/O dominated).
    def _one(row):
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

        reads_logger.info("**** Calculating UJCs...")
                
        # Build every multi-exon transcript's junction chain directly from the
        # corrected GTF in a single streaming awk pass (see _UJC_AWK above). This
        # replaces the gffread -> gtftools -> bedtools pipeline whose gffread step
        # alone took ~11 h on a 6 GB GTF and whose gtftools step took ~75 min/sample.
        # Nothing downstream reads the (previously gffread-rewritten) corrected GTF —
        # only the classification and junction files are loaded — so it is no longer
        # copied or rewritten.
        input_gtf = f"{inputPathPrefix}_corrected.gtf"
        if not os.path.exists(input_gtf):
            reads_logger.error(f"Input GTF not found: {input_gtf}")
            sys.exit(-1)

        # LC_ALL=C keeps gawk's byte ops (index/substr) out of slow multibyte-locale
        # paths; the awk is otherwise I/O-bound on the GTF read.
        ujc_cmd = "LC_ALL=C awk '" + _UJC_AWK + f"' {input_gtf} > {outputPathPrefix}tmp_UJC.txt"
        try:
            subprocess.check_call(ujc_cmd, shell=True)
        except subprocess.CalledProcessError:
            reads_logger.error(f"ERROR building UJC chains (awk) from {input_gtf}")
            sys.exit(-1)

        ## Pandas merge to the left
        input_classfile = f"{inputPathPrefix}_classification.txt"
        if not os.path.exists(input_classfile):
             reads_logger.error(f"Input classification file not found: {input_classfile}")
             sys.exit(-1)
             
        # Classification columns: 0=isoform, 1=chrom, 2=start, 3=end, 4=strand,
        # 7=structural_category, 11=associated_transcript. (Earlier code mislabeled
        # col 2=start as "strand" and col 7=structural_category as "associated_transcript",
        # which corrupted the monoexon UJC — every monoexon got a unique, unshareable hash.)
        clas_df = pd.read_csv(input_classfile, sep = "\t", usecols = [0, 1, 2, 3, 4, 7, 11], dtype = "str")
        clas_df.columns = ["isoform", "chr", "start", "end", "strand", "structural_category", "associated_transcript"]
        ujc_df = pd.read_csv(f"{outputPathPrefix}tmp_UJC.txt", sep = "\t", names = ["isoform", "jxn_string"], dtype = "str")

        merged_df = pd.merge(clas_df, ujc_df, on = "isoform", how = "left")

        # Fill in the UJC key for monoexons (multi-exon isoforms already have a
        # junction-chain jxn_string; monoexons have none). Vectorised over just the
        # monoexon subset instead of a per-row Python apply across every read — on
        # multi-million-read samples that apply dominated the "Calculating UJCs" step.
        # Output is identical to the old row-wise keying:
        #  - novel monoexon (no real associated_transcript) -> keyed by genomic span,
        #    so distinct loci stay distinct instead of collapsing to one novel UJC;
        #  - annotated monoexon -> keyed by reference transcript, so the same
        #    monoexonic transcript is one shared UJC across samples.
        mono = merged_df["jxn_string"].isna()
        if mono.any():
            sub = merged_df.loc[mono]
            at = sub["associated_transcript"]
            is_novel = at.isna() | at.isin(("novel", ""))
            novel_key = (sub["chr"] + "_" + sub["strand"] + "_" + sub["start"] + "_"
                         + sub["end"] + "_monoexon_" + sub["structural_category"])
            annot_key = (sub["chr"] + "_" + sub["strand"] + "_monoexon_"
                         + at.fillna("") + "_" + sub["structural_category"])
            merged_df.loc[mono, "jxn_string"] = novel_key.where(is_novel, annot_key)

        # Hash only the DISTINCT junction strings and broadcast back to every read.
        # Reads that share a UJC (the common case) collapse to a single sha256 call
        # instead of one per read — the dominant cost on large samples. Identical output.
        _hmap = {s: hashlib.sha256(s.encode("utf-8")).hexdigest()
                 for s in merged_df["jxn_string"].unique()}
        merged_df["jxnHash"] = merged_df["jxn_string"].map(_hmap)

        # Write only the two columns the paste actually appends (jxn_string, jxnHash);
        # the merge preserves the classification's row order, so a positional paste is
        # correct. Dropping the other four columns (incl. the long per-read isoform IDs)
        # avoids writing millions of discarded values to the temp file.
        merged_df[["jxn_string", "jxnHash"]].to_csv(f"{outputPathPrefix}_temp.txt", index = False, sep = "\t")

        cmd_paste = f"""bash -c 'paste <(cat {input_classfile} | tr -d '\r') <(cut -f 1,2 {outputPathPrefix}_temp.txt | tr -d '\r') > {outputPathPrefix}_reads_classification.txt'"""
        subprocess.call(cmd_paste, shell = True)

        os.remove(f"{outputPathPrefix}tmp_UJC.txt")
        os.remove(f"{outputPathPrefix}_temp.txt")

    _run_parallel(_one, [row for _, row in df.iterrows()], getattr(args, 'jobs', 1))

def main():
    global utilitiesPath
    global sqantiqcPath

    args = reads_argparser().parse_args()

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

    # Run plotting script directly as a function call
    if args.SKIPPLOTS:
        reads_logger.info("--skip_plots set: skipping SQANTI-reads tables and plots generation.")
        return

    reads_logger.info("Running SQANTI-reads tables and plots generation...")

    prefix = args.PREFIX if args.PREFIX else "sqantiReads"

    from src.utilities.sqanti_reads_plots import run_reads_plots
    from src.utilities.sqanti_reads_config import load_config
    config = load_config(args.CONFIG)
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
        report=args.report,
        config=config,
        jobs=args.jobs,
    )


if __name__ == "__main__":
    main()
