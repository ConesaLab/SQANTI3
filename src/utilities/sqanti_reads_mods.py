import sys
import pysam
import os
import glob
import pandas as pd
import numpy as np
from src.commands import MINIMAP2_MODS_CMD, run_command
from src.helpers import get_corr_filenames

# For translation from one mod name to another
ACCEPTED = ["m5C","m6A","inosine","pseU","2OmeC","2OmeA","2OmeG","2OmeU"]
mod_lookup = {
    # m5C
    "5-Methylcytosine": "m5C","C+m": "m5C","27551": "m5C",
    # m6A
    "N(6)-Methyladenosine": "m6A","A+a": "m6A","21891": "m6A",
    # inosine
    "Inosine": "inosine","A+17596": "inosine","17596": "inosine",
    # pseU
    "Pseudouridine": "pseU","T+17802": "pseU","17802": "pseU",
    # 2OmeC
    "2'-O-methylcytidine": "2OmeC","C+19228": "2OmeC","19228": "2OmeC",
    # 2OmeA
    "2'-O-methyladenosine": "2OmeA","A+69426": "2OmeA","69426": "2OmeA",
    # 2OmeG
    "2'-O-methylguanosine": "2OmeG","G+19229": "2OmeG","19229": "2OmeG",
    # 2OmeU
    "2'-O-methyluridine": "2OmeU","T+19227": "2OmeU","19227": "2OmeU"
    }

CHEBI_to_acc = {"m":"m5C","a":"m6A", "19228":"2OmeC", 
                 "17802":"pseU","19227":"2OmeU", "17596":"inosine" , 
                 "69426":"2OmeA", "19229":"2OmeG"}
acc_to_CHEBI = {v: k for k, v in CHEBI_to_acc.items()}

def lookup_mods(args):
    # Function for mod_types parsing
    types = args.MODTYPES
    all_correct = True
    for i in range(len(types)):
        if types[i] not in ACCEPTED:
            # We look if it's in our dictionary
            if types[i] in mod_lookup.keys():
                types[i] = mod_lookup[types[i]]
            else:
                all_correct = False
    return types, all_correct
                
def validate_MODS(args,reads_logger):
    len_probs = len(args.MODPROB)
    len_mods = len(args.MODTYPES)
    if len_probs != 1 and len_probs != len_mods:
        reads_logger.error(f"[ERROR] The length of the probabilities vector must be equal to the number of modification types or 1 (global probability).")
        sys.exit(-1)
    for i in args.MODPROB:
        if i < 0.5 or i > 1:
            reads_logger.error("[ERROR] The probabilities have to be a number between 0.5 and 1")
            sys.exit(-1)
    if args.MODTYPES == []:
        reads_logger.warning("[WARNING] No modification types given. Every modifications in the uBAMs will be included")
    if args.sqanti_dirs:
        reads_logger.error("[ERROR] Modifications mode cannot be run in Fast mode.")
        sys.exit(-1)
    types, all_correct = lookup_mods(args)
    if all_correct:
        args.MODTYPES = types
    else:
        reads_logger.error("[ERROR] At least one of the given modifications is not accepted (accepted: m5C,m6A,inosine, pseU,2OmeC,2OmeA,2OmeG,2OmeU or their CHEBI code)")
        sys.exit(-1)
    if len(set(args.MODTYPES)) != len(args.MODTYPES):
        reads_logger.error("[ERROR] The modification types given are not unique")
        sys.exit(-1)

def check_uBAM(path, reads_logger,n=1000):
    '''
    Function to check if first n reads of a ubam have MM/ML tags.
    It also checks if the bam is aligned and returns the modifications found in the first n reads
    '''
    with pysam.AlignmentFile(path,"rb", check_sq=False) as bam:
        if bam.nreferences > 0:
            reads_logger.error("[ERROR] BAM seems to have been aligned. uBAM is required" )
            sys.exit(-1)
        types_found = set()
        types_new = set()
        for i, read in enumerate(bam):
            if i >= n:
                break
            try:
                mods = read.modified_bases
            except (ValueError, RuntimeError):
                continue
            if not mods:
                continue
            # for every modification found in the read, we check if it's in our dictionary and add it to the set of found modifications
            for base, strand, code in mods:
                key =f"{base}+{code}"
                name = mod_lookup.get(key, None)
                if name:
                    types_found.add(name)
                else:
                    types_new.add(key)
        if types_new:
            reads_logger.warning(f"[WARNING] The following modifications were found in the uBAM but are not supported, they will be ignored: {', '.join(sorted(types_new))}") 
    return types_found

def mod_files_check(args, file_acc, reads_logger):
    """Checks that all the files needed for mods mode are in order"""
    # Check for uBAM
    uBAM_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
    try:
        uBAM_files = glob.glob(uBAM_pattern)[0]
    except IndexError:
        reads_logger.error(f"[ERROR] {uBAM_pattern} does not exist")
        sys.exit(-1)
    else:
        if os.path.isfile(args.refFasta) is False:
            reads_logger.error(f'[ERROR] You inputted uBAM files to map SQANTI3 but no reference genome FASTA')
            sys.exit(-1)
        if os.path.isfile(args.refGTF) is False:
            reads_logger.error(f'[ERROR] You inputted uBAM files to map but no reference annotation GTF')
            sys.exit(-1)
        if args.verbose:
            reads_logger.debug(f'[INFO] You inputted bam files, we will run sqanti_reads in modifications mode for sample {uBAM_files}')
        mods = check_uBAM(uBAM_files,reads_logger)
        if not mods:
            reads_logger.error(f"[ERROR] While checking first 1000 reads, no MM/ML tags were found")
            sys.exit(-1)
        if not set(args.MODTYPES).issubset(mods):
            reads_logger.error(f"[ERROR] The following modifications were requested but not found in the uBAM: {', '.join(sorted(set(args.MODTYPES) - mods))}")
            sys.exit(-1)
        return mods, uBAM_files

def align_uBAM(args, uBAM_file, file_acc, sampleID, reads_logger):
    '''Converts the uBAM to fastq keeping the MM/ML tags, aligns it with minimap2 into the
    SAM that SQANTI3 QC will reuse and generates a sorted and indexed BAM for the mods module.
    Returns the fastq path'''
    outdir = os.path.join(args.OUTPUT, file_acc)
    logdir = os.path.join(outdir, "logs")
    os.makedirs(logdir, exist_ok=True)

    # Same names SQANTI3 QC uses, so it finds and reuses our SAM
    corrGTF, corrSAM, corrFASTA, _, _ = get_corr_filenames(outdir, sampleID)
    fastq = os.path.join(outdir, f"{sampleID}_mods_reads.fastq")
    mods_bam = os.path.join(outdir, f"{sampleID}_mods.bam")

    if all(os.path.isfile(f) for f in (fastq, corrSAM, mods_bam, f"{mods_bam}.bai")):
        reads_logger.info(f"Aligned uBAM for {sampleID} already exists. Using it...")
        return fastq

    # -T so that tags go to the read header
    reads_logger.info(f"**** Converting uBAM to fastq for {sampleID}...")
    cmd = f"samtools fastq -T MM,ML,MN,Mm,Ml {uBAM_file} > {fastq}.tmp"
    run_command(cmd, reads_logger, os.path.join(logdir, "samtools_fastq.log"), description="converting uBAM to fastq")
    os.replace(f"{fastq}.tmp", fastq)

    # Write to a temporary file so a failed alignment never leaves a SAM that QC would reuse
    reads_logger.info(f"**** Aligning {sampleID} reads with Minimap2...")
    cmd = MINIMAP2_MODS_CMD.format(cpus=args.cpus, g=args.refFasta, i=fastq, o=f"{corrSAM}.tmp")
    run_command(cmd, reads_logger, os.path.join(logdir, "minimap2_mods_alignment.log"), description="aligning uBAM reads")
    os.replace(f"{corrSAM}.tmp", corrSAM)

    # A corrected fasta from a previous run makes QC skip the correction and ignore the new SAM
    for f in (corrFASTA, corrGTF):
        if os.path.isfile(f):
            reads_logger.warning(f"[WARNING] Removing {f} from a previous run, it does not match the new alignment")
            os.remove(f)

    reads_logger.info(f"**** Sorting and indexing {sampleID} alignment...")
    cmd = f"samtools sort -@ {args.cpus} -O bam -o {mods_bam}.tmp {corrSAM}"
    run_command(cmd, reads_logger, os.path.join(logdir, "samtools_sort.log"), description="sorting aligned uBAM reads")
    os.replace(f"{mods_bam}.tmp", mods_bam)
    cmd = f"samtools index {mods_bam}"
    run_command(cmd, reads_logger, os.path.join(logdir, "samtools_index.log"), description="indexing aligned uBAM reads")

    return fastq

def get_mods_info(args, df, reads_logger):
    '''
    To check in every sample if the uBAM is valid and has the modifications.
    It also saves the probabilities for each mod type in a dictionary.
    '''
    uBAM_files = {}
    types_found = set()
    # For each sample we check and add all modifications found
    for _, row in df.iterrows():
        mods, uBAM_files[row['file_acc']] = mod_files_check(args, row['file_acc'], reads_logger)
        types_found |= mods

    # If user does not give modtypes 
    if args.MODTYPES == []:
        # We use the mods found in all of the samples
        args.MODTYPES = sorted(types_found)
        
    # if a single probability is used for every modification type
    probs = args.MODPROB * len(args.MODTYPES) if len(args.MODPROB) == 1 else args.MODPROB
    args.MODTOTHRESH = dict(zip(args.MODTYPES, probs))

    reads_logger.info("Modifications to analyse (threshold): " +
                      ", ".join(f"{mod} ({prob})" for mod, prob in args.MODTOTHRESH.items()))
    return uBAM_files

def convert_prob(prob):
    '''
    We convert the prob of the mods (from 0 to 256) to a number between 0 and 1
    '''
    return (prob+0.5) / 256

def new_row(read_id, read, args):
    row = {"isoform": read_id, "read_length": read.query_length, 
           "start": read.reference_start + 1, "end": read.reference_end, # 1-based
           "chrom": read.reference_name, "strand": "-" if read.is_reverse else "+"}
    for mod in args.MODTYPES:
        row[f"{mod}_positions"] = ""
        row[f"{mod}_probs"] = ""
        row[f"{mod}_n_mod"] = 0
    return row

def process_read(read,args,reads_logger):
    '''
    Function to process a read and return the modifications that pass the probability threshold
    '''
    read_id = read.query_name
    mods = read.modified_bases
    # We create the row with the read id, start, end, ... and a list for each mod with pos and probs 
    row = new_row(read_id, read, args)
    # We will have to convert read pos to genomic 
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    read_to_genomic = {read_pos: ref_pos for read_pos, ref_pos in aligned_pairs}
    # For every mod (canonical base, strand, CHEBI code) with its [(pos,qual),...] list
    for tupla, probs_list in mods.items():
        valid_mods = [] # mods that pass the probability threshold
        base, strand, code = tupla
        code = str(code) # because CHEBI codes are interpreted as integers 
        if CHEBI_to_acc.get(code) not in args.MODTYPES:
            continue
        mod_name = CHEBI_to_acc.get(code)
        prob_lim = args.MODTOTHRESH.get(mod_name)
        # For each position we check that the probability is above the threshold
        for pos, prob in probs_list:
            pos = read_to_genomic.get(pos,None) # convert to genomic position
            if pos is not None:
                pos = pos + 1 # convert to 1-based
                prob = convert_prob(prob)
                if prob >= prob_lim:
                    valid_mods.append((pos, prob))
        valid_mods.sort(key=lambda x: x[0]) # we make sure list is sorted by position
        row[f"{mod_name}_probs"] = '_'.join([str(round(prob,3)) for _, prob in valid_mods])
        row[f"{mod_name}_positions"] = '_'.join([str(pos) for pos, _ in valid_mods])
        row[f"{mod_name}_n_mod"] = len(valid_mods)
    return row

def gen_mod_tsv(args,file_acc,sampleID, reads_logger):
    '''
    Function to extract mod information from the aligned BAM to a TSV
    '''
    path_to_bam = os.path.join(args.OUTPUT, file_acc, f"{sampleID}_mods.bam")
    rows = []
    with pysam.AlignmentFile(path_to_bam, "rb", threads=1) as bamfile:
        if not bamfile.check_index():
            reads_logger.error("Error: The BAM file does not have an index (.bai). Aborting.")
            sys.exit(-1)
        
        for read in bamfile.fetch():
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            row = process_read(read,args,reads_logger)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.OUTPUT, file_acc, f"{sampleID}_modifications.tsv"), sep="\t", index=False)
    return str(os.path.join(args.OUTPUT,file_acc, f"{sampleID}_modifications.tsv"))

def explode_mods(df,args,mod):
    df = df.copy()
    df = df[df[f"{mod}_n_mod"] > 0]
    df[f"{mod}_positions"] = df[f"{mod}_positions"].str.split(",")
    df[f"{mod}_probs"] = df[f"{mod}_probs"].str.split(",")
    df = df[["isoform", "jxnHash", "chrom", "strand", f"{mod}_positions", f"{mod}_probs"]]
    df_ex = df.explode([f"{mod}_positions" ,f"{mod}_probs"])
    df_ex = df_ex.rename(columns={f"{mod}_positions": "mod_pos", f"{mod}_probs": "mod_prob"})
    # We have to turn positions from str to numeric type to order them later
    df_ex["mod_pos"] = df_ex["mod_pos"].astype(int)
    df_ex["mod_prob"] = df_ex["mod_prob"].astype(float)
    df_ex["mod_type"] = mod
    
    return df_ex    

def ujc_mod_sites(args, file_acc, sampleID, reads_logger):
    '''
    Funtion to collapse modification information to UJC level and adding to classification table
    '''
    # The first step will be to add to the tsv the jxnHash 
    tsv_path = os.path.join(args.OUTPUT,file_acc, f"{sampleID}_modifications.tsv")
    clasif_path = os.path.join(args.OUTPUT, file_acc, f"{sampleID}_reads_classification.txt")
    mod_clasif_path = os.path.join(args.OUTPUT, file_acc, f"{sampleID}_reads_mods_classification.txt")
    read_mods = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    num_cols = ["start", "end", "read_length"] + [f"{mod}_n_mod" for mod in args.MODTYPES]
    read_mods[num_cols] = read_mods[num_cols].astype(int)
    clasif = pd.read_csv(clasif_path, sep="\t", usecols=["isoform", "jxnHash"])
    read_mods = read_mods.merge(clasif, on="isoform", how="left")
    # Spans of ALL the reads of each UJC (I think dropna is not needed but just in case)
    # We generate a df with the start and end positions of all reads that belong to the same jxnHash
    spans = read_mods.dropna(subset=["jxnHash"]).groupby("jxnHash").agg(
        starts=("start", lambda s: sorted(s)),
        ends=("end", lambda s: sorted(s)),
        n_reads_ujc=("isoform", "size"))
    # We can use this to count the total number of reads that cover a position in the UJC
    exploded_tables = []
    for i in args.MODTYPES:
        read_mods_i = explode_mods(read_mods, args, i)
        exploded_tables.append(read_mods_i)
    exploded_mods = pd.concat(exploded_tables, ignore_index=True)
    # And now we can group by jxnHash, mod_type and mod_pos to count the number of reads that have that modification at that position
    sites = (exploded_mods
         .groupby(["jxnHash", "chrom", "strand", "mod_type", "mod_pos"], as_index=False)
         # For each group of same hash, modtype and position we extract number of mods and mean of probability 
         .agg(n_mod=("isoform", "size"), mean_prob=("mod_prob", "mean")))
    # Now we merge with spans based on jxnHash
    sites = sites.merge(spans, on="jxnHash", how="left")
    # For a same UJC, number of reads that cover a position p is:
    # Nreads with start <= p - Nreads with end < p --> Nreads that start before the mod - Nreads that end before the mod
    sites["n_cov"] = [np.searchsorted(st, p, side="right") - np.searchsorted(en, p, side="left") 
                      for st, en, p in zip(sites["starts"], sites["ends"], sites["mod_pos"])]
    sites["occupancy"] = sites["n_mod"]/sites["n_cov"]
    sites = sites.drop(columns=["starts","ends"]).sort_values(["jxnHash", "mod_type", "mod_pos"])
    fmt = sites.assign(mod_pos=sites["mod_pos"].astype(str),
                   occupancy=sites["occupancy"].round(3).astype(str),
                   n_mod=sites["n_mod"].astype(str),
                   n_cov=sites["n_cov"].astype(str))
    # We unstack our columns based on mod_types
    ujc = (fmt.groupby(["jxnHash", "mod_type"])
          .agg(positions=("mod_pos", ",".join),
               occupancy=("occupancy", ",".join),
               n_mod=("n_mod", ",".join),
               n_cov=("n_cov", ",".join))
          .unstack("mod_type"))
    # Columns are now named like (positions,m6A) so we change it
    ujc.columns = [f"ujc_{mod}_{field}" for field, mod in ujc.columns]
    # Now we only have to fill nas of 
    cols = [f"ujc_{mod}_{field}" for mod in args.MODTYPES for field in ("positions", "occupancy", "n_mod", "n_cov")]
    ujc = ujc.reindex(columns=cols).fillna("").reset_index()
    # Each read gets the modified sites of its UJC and their occupancy
    read_mods = read_mods.merge(ujc, on="jxnHash", how="left")
    ujc_cols = [col for col in ujc.columns if col != "jxnHash"]
    read_mods[ujc_cols] = read_mods[ujc_cols].fillna("")
    # Reads table with the UJC info and long table with one row per UJC site
    out_prefix = os.path.join(args.OUTPUT, file_acc, sampleID)
    read_mods.to_csv(f"{out_prefix}_reads_modifications_ujc.tsv", sep="\t", index=False)
    sites.to_csv(f"{out_prefix}_ujc_mod_sites.tsv", sep="\t", index=False)
    # Modification columns are added to the classification
    clasif_completa = pd.read_csv(clasif_path, sep="\t", dtype=str, keep_default_na=False)
    cols_tsv = [col for mod in args.MODTYPES for col in read_mods.columns
                if col.startswith(f"{mod}_") or col.startswith(f"ujc_{mod}_")]
    # Columns from a previous run would be duplicated
    clasif_completa = clasif_completa.drop(columns=[col for col in cols_tsv if col in clasif_completa.columns])
    n_rows = len(clasif_completa)
    # Now we merge everything to the classification
    clasif_completa = clasif_completa.merge(read_mods[["isoform"] + cols_tsv].astype(str),
                                            on="isoform", how="left")
    if len(clasif_completa) != n_rows:
        reads_logger.error(f"[ERROR] Classification rows changed after adding modifications for {sampleID}")
        sys.exit(-1)
    clasif_completa.fillna("").to_csv(f"{mod_clasif_path}.tmp", sep="\t", index=False)
    os.replace(f"{mod_clasif_path}.tmp", mod_clasif_path)

    reads_logger.info(f"{sampleID}: {len(sites)} modified UJC sites written to {out_prefix}_ujc_mod_sites.tsv")
    return mod_clasif_path