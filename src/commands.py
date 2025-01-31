import logging
import os
import sys
import shutil
import subprocess

from .logging_config import MY_LOGGING_CONFIG, qc_logger

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)

GMAP_CMD = "gmap --cross-species -n 1 --max-intronlength-middle=2000000 --max-intronlength-ends=2000000 -L 3000000 -f samse -t {cpus} -D {dir} -d {name} -z {sense} {i} > {o}"
MINIMAP2_CMD = "minimap2 -ax splice --secondary=no -C5 -u{sense} -t {cpus} {g} {i} > {o}"
DESALT_CMD = "deSALT aln {dir} {i} -t {cpus} -x ccs -o {o}"
ULTRA_CMD = "uLTRA pipeline {g} {a} {i} {o_dir} --t {cpus} --prefix {prefix} --isoseq"

# GTF
GTF2GENEPRED_PROG = os.path.join(utilitiesPath,"gtfToGenePred")
GFFREAD_PROG = "gffread"

if shutil.which(GTF2GENEPRED_PROG) is None:
    qc_logger.info(f"Cannot find executable {GTF2GENEPRED_PROG}. Abort!")
    sys.exit(1)
if shutil.which(GFFREAD_PROG) is None:
    qc_logger.error("Cannot find executable {GFFREAD_PROG}. Abort!")
    sys.exit(1)


# Rscript
RSCRIPTPATH = shutil.which('Rscript')

RSCRIPT_REPORT = '/report_qc/SQANTI3_report.R'

ISOANNOT_PROG =  os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")

def get_aligner_command(aligner_choice, genome, isoforms, annotation, 
                        outdir, corrSAM, n_cpu, gmap_index, sense):
    # Even though the speed does not change form the ifelse, this is cleaner
    match aligner_choice:
        case "gmap":
            qc_logger.info("****Aligning reads with GMAP...")
            cmd = GMAP_CMD.format(
                cpus=n_cpu,
                dir=os.path.dirname(gmap_index),
                name=os.path.basename(gmap_index),
                sense=sense,
                i=isoforms,
                o=corrSAM,
            )
        case "minimap2":
            qc_logger.info("****Aligning reads with Minimap2...")
            cmd = MINIMAP2_CMD.format(
                cpus=n_cpu,
                sense=sense,
                g=genome,
                i=isoforms,
                o=corrSAM,
            )
        case "deSALT":
            qc_logger.info("****Aligning reads with deSALT...")
            cmd = DESALT_CMD.format(
                cpus=n_cpu,
                dir=gmap_index,
                i=isoforms,
                o=corrSAM,
            )
        case "uLTRA":
            qc_logger.info("****Aligning reads with uLTRA...")
            cmd = ULTRA_CMD.format(
                cpus=n_cpu,
                prefix="../" + os.path.splitext(os.path.basename(corrSAM))[0],
                g=genome,
                a=annotation,
                i=isoforms,
                o_dir=outdir + "/uLTRA_out/",
            )
        case _:
            qc_logger.error(f"Unsupported aligner choice: {aligner_choice}")
            raise ValueError()
    return cmd

# To dynamically get the utilities path
def get_gmst_prog(utilities_path):
    return os.path.join(utilities_path, "gmst", "gmst.pl")

GMST_CMD = f"perl {get_gmst_prog(utilitiesPath)} -faa --strand direct --fnn --output {{o}} {{i}}"

def run_gmst(corrFASTA,orf_input,gmst_pre):
    if orf_input is not None:
        qc_logger.info(f"Running ORF prediction of input on {orf_input}...")
        cmd = GMST_CMD.format(i=os.path.realpath(orf_input), o=gmst_pre)
    else:
        cmd = GMST_CMD.format(i=corrFASTA, o=gmst_pre)
    cmd = f"cd {os.path.dirname(gmst_pre)}; {cmd}"
    run_command(cmd, description="GMST ORF prediction")

def GTF_to_genePred(corrGTF):
    """
    Converts a GTF file to a genePred file.
    """
    queryFile = os.path.splitext(corrGTF)[0] +".genePred"
    if os.path.exists(queryFile):
        qc_logger.info(f"{queryFile} already exists. Using it.")
    else:
        # gtf to genePred
        cmd = f"{GTF2GENEPRED_PROG} {corrGTF} {queryFile} -genePredExt -allErrors -ignoreGroupsWithoutExons"
        run_command(cmd, "GTF to genePred conversion")
    return queryFile
   

def run_command(cmd, out_file='log/program.log',description="command execution"):
    """
    Executes a shell command and handles errors gracefully.
    
    :param cmd: The command to execute (string).
    :param description: A short description of the operation for better error messages (default: "command execution").
    :raises SystemExit: Exits the script if the command fails.
    """
    try:
        qc_logger.debug(out_file)
        MY_LOGGING_CONFIG['handlers']['process_handler']['filename'] = out_file
        logging.config.dictConfig(MY_LOGGING_CONFIG)
        process_logger = logging.getLogger('process_logger')

        result = subprocess.run(cmd, shell=True,capture_output=True,
                                check=True,encoding="utf-8")
        
        qc_logger.debug(f"Process {cmd} had {result.returncode}")
        qc_logger.debug(f"Returncode {result.returncode} class is {type(result.returncode)}")
        process_logger.info(result.stdout)
        if result.returncode !=0:
            qc_logger.debug("EROOOOOOR")
            raise BrokenPipeError(result.stderr)
    except subprocess.CalledProcessError as e:
        process_logger.error(f"Details: {e}")
        process_logger.info(e.stdout)
        process_logger.error(e.stderr)
        qc_logger.error(f"Something went wrong during {description}")
        qc_logger.error(f"For more inflo, check {out_file}")
        sys.exit(1)
    except BrokenPipeError as e:

        sys.exit(1)
