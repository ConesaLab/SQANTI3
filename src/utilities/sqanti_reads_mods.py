import sys
import pysam
import os
import glob
from src.commands import MINIMAP2_MODS_CMD, run_command
from src.helpers import get_corr_filenames
ACCEPTED = ["m5C","m6A","inosine","pseU","2OmeC","2OmeA","2OmeG","2OmeU"]
mod_lookup = {
    # m5C
    "5-Methylcytosine": "m5C",
    "C+m": "m5C",
    "27551": "m5C",
    # m6A
    "N(6)-Methyladenosine": "m6A",
    "A+a": "m6A",
    "21891": "m6A",
    # inosine
    "Inosine": "inosine",
    "A+17596": "inosine",
    "17596": "inosine",
    # pseU
    "Pseudouridine": "pseU",
    "T+17802": "pseU",
    "17802": "pseU",
    # 2OmeC
    "2'-O-methylcytidine": "2OmeC",
    "C+19228": "2OmeC",
    "19228": "2OmeC",
    # 2OmeA
    "2'-O-methyladenosine": "2OmeA",
    "A+69426": "2OmeA",
    "69426": "2OmeA",
    # 2OmeG
    "2'-O-methylguanosine": "2OmeG",
    "G+19229": "2OmeG",
    "19229": "2OmeG",
    # 2OmeU
    "2'-O-methyluridine": "2OmeU",
    "T+19227": "2OmeU",
    "19227": "2OmeU"
    }

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
        if i < 0 or i > 1:
            reads_logger.error("[ERROR] The probabilities have to be a number between 0 and 1")
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
        reads_logger.error("[ERRROR] At least one of the given modifications is not accepted (accepted: m5C,m6A,inosine, pseU,2OmeC,2OmeA,2OmeG,2OmeU or their CHEBI code)")
        sys.exit(-1)

def check_uBAM(path, reads_logger,n=1000):
    '''Function to check if first n reads of a bam have MM/ML tags'''
    bam = pysam.AlignmentFile(path,"rb", check_sq=False)
    if bam.nreferences > 0:
        reads_logger.error("[ERROR] BAM seems to have been aligned. uBAM is required" )
        sys.exit(-1)
    reads_checked = 0
    reads_with_mods = 0
    for read in bam:
        if reads_checked >= n:
            break
        tags = dict(read.tags)
        if read.has_tag("MM") and read.has_tag("ML"):
            reads_with_mods +=1
        reads_checked += 1
    bam.close()
    return reads_with_mods > 0

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
        
        if not check_uBAM(uBAM_files,reads_logger):
            reads_logger.error(f"[ERROR] While checking first {1000} reads, no MM/ML tags where found")
            sys.exit(-1)
        return True, uBAM_files

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

    # -y so that tags go to the read header
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