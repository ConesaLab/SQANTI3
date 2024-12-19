import os
import sys
import shutil
import subprocess

from .parsers import STARcov_parser
from .utilities.short_reads import star

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)

GMAP_CMD = "gmap --cross-species -n 1 --max-intronlength-middle=2000000 --max-intronlength-ends=2000000 -L 3000000 -f samse -t {cpus} -D {dir} -d {name} -z {sense} {i} > {o}"
MINIMAP2_CMD = "minimap2 -ax splice --secondary=no -C5 -u{sense} -t {cpus} {g} {i} > {o}"
DESALT_CMD = "deSALT aln {dir} {i} -t {cpus} -x ccs -o {o}"
ULTRA_CMD = "uLTRA pipeline {g} {a} {i} {o_dir} --t {cpus} --prefix {prefix} --isoseq"

# To dynamically get the utilities path
def get_gmst_prog(utilities_path):
    return os.path.join(utilities_path, "gmst", "gmst.pl")


GMST_CMD = f"perl {get_gmst_prog(utilitiesPath)} -faa --strand direct --fnn --output {{o}} {{i}}"

# GTF
GTF2GENEPRED_PROG = os.path.join(utilitiesPath,"gtfToGenePred")
GFFREAD_PROG = "gffread"

if shutil.which(GTF2GENEPRED_PROG) is None:
    print("Cannot find executable {0}. Abort!".format(GTF2GENEPRED_PROG), file=sys.stderr)
    sys.exit(1)
if shutil.which(GFFREAD_PROG) is None:
    print("Cannot find executable {0}. Abort!".format(GFFREAD_PROG), file=sys.stderr)
    sys.exit(1)


# Rscript
RSCRIPTPATH = shutil.which('Rscript')

RSCRIPT_REPORT = '/report_qc/SQANTI3_report.R'

ISOANNOT_PROG =  os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")


def short_reads_mapping(args):
    """
    Runs STAR for mapping short reads to the transcriptome
    """
    if args.short_reads is not None:
        print("**** Running STAR for calculating Short-Read Coverage.", file=sys.stdout)
        star_out, star_index = star(args.genome, args.short_reads, args.dir, args.cpus)
        SJcovNames, SJcovInfo = STARcov_parser(star_out)
    else:
        star_out, star_index, SJcovNames, SJcovInfo = None, None, None, None
    return(star_out, star_index, SJcovNames, SJcovInfo)


def get_aligner_command(aligner_choice, genome, isoforms, annotation, 
                        outdir, corrSAM, n_cpu, gmap_index, sense):
    # Even though the speed does not change form the ifelse, this is cleaner
    match aligner_choice:
        case "gmap":
            print("****Aligning reads with GMAP...", file=sys.stdout)
            cmd = GMAP_CMD.format(
                cpus=n_cpu,
                dir=os.path.dirname(gmap_index),
                name=os.path.basename(gmap_index),
                sense=sense,
                i=isoforms,
                o=corrSAM,
            )
        case "minimap2":
            print("****Aligning reads with Minimap2...", file=sys.stdout)
            cmd = MINIMAP2_CMD.format(
                cpus=n_cpu,
                sense=sense,
                g=genome,
                i=isoforms,
                o=corrSAM,
            )
        case "deSALT":
            print("****Aligning reads with deSALT...", file=sys.stdout)
            cmd = DESALT_CMD.format(
                cpus=n_cpu,
                dir=gmap_index,
                i=isoforms,
                o=corrSAM,
            )
        case "uLTRA":
            print("****Aligning reads with uLTRA...", file=sys.stdout)
            cmd = ULTRA_CMD.format(
                cpus=n_cpu,
                prefix="../" + os.path.splitext(os.path.basename(corrSAM))[0],
                g=genome,
                a=annotation,
                i=isoforms,
                o_dir=outdir + "/uLTRA_out/",
            )
        case _:
            raise ValueError(f"Unsupported aligner choice: {aligner_choice}")
    return cmd

def run_gmst(corrFASTA,orf_input,gmst_pre):
    if orf_input is not None:
        print("Running ORF prediction of input on {0}...".format(orf_input))
        cmd = GMST_CMD.format(i=os.path.realpath(orf_input), o=gmst_pre)
    else:
        cmd = GMST_CMD.format(i=corrFASTA, o=gmst_pre)
    run_command(cmd, description="GMST ORF prediction")


def run_command(cmd, description="command execution"):
    """
    Executes a shell command and handles errors gracefully.
    
    :param cmd: The command to execute (string).
    :param description: A short description of the operation for better error messages (default: "command execution").
    :raises SystemExit: Exits the script if the command fails.
    """
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR during {description}: {cmd}", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)

