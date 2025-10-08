import logging
import os
import sys
import shutil
import subprocess
import platform

from src.module_logging import qc_logger
from src.logging_config import MAIN_LOGGING_CONFIG
from src.config import utilitiesPath

sys.path.insert(0, utilitiesPath)

GMAP_CMD = "gmap --cross-species -n 1 --max-intronlength-middle=2000000 --max-intronlength-ends=2000000 -L 3000000 -f samse -t {cpus} -D {dir} -d {name} -z sense-force {i} > {o}"
MINIMAP2_CMD = "minimap2 -ax splice --secondary=no -C5 -uf -t {cpus} {g} {i} > {o}"
DESALT_CMD = "deSALT aln {dir} {i} -t {cpus} -x ccs -o {o}"
ULTRA_CMD = "uLTRA pipeline {g} {a} {i} {o_dir} --t {cpus} --prefix {prefix} --isoseq"

# GTF - Platform-specific binary selection
def _get_gtftogenepred_binary():
    """
    Select the correct gtfToGenePred binary based on the platform.
    Supports Linux and macOS (Darwin). Raises an error for unsupported platforms.
    """
    system = platform.system()
    if system == "Linux":
        binary_name = "gtfToGenePred-linux-x86_64"
    elif system == "Darwin":
        binary_name = "gtfToGenePred-darwin-x86_64"
    else:
        raise OSError(f"Unsupported platform: {system}. gtfToGenePred binary is only available for Linux and macOS.")

    binary_path = os.path.join(utilitiesPath, binary_name)
    if not os.path.exists(binary_path):
        raise FileNotFoundError(f"gtfToGenePred binary not found at {binary_path}")

    return binary_path

GTF2GENEPRED_PROG = _get_gtftogenepred_binary()
GFFREAD_PROG = "gffread"

# Rscript QC
RSCRIPTPATH = shutil.which('Rscript')
RSCRIPT_QC_REPORT = os.path.join(utilitiesPath,"report_qc","SQANTI3_report.R")
RSCRIPT_TUSCO_REPORT = os.path.join(utilitiesPath,"report_qc","TUSCO_report.R")

# Rscript filter
RSCRIPT_FILTER_REPORT = os.path.join(utilitiesPath,"report_filter","SQANTI3_filter_report.R")
RSCRIPT_ML = os.path.join(utilitiesPath,"filter","SQANTI3_MLfilter.R")


#Rscript rescue
RESCUE_AUTO_PATH = os.path.join(utilitiesPath,"rescue","automatic_rescue.R")
RSCRIPT_RESCUE_RULES = os.path.join(utilitiesPath, "rescue","rescue_by_mapping_rules.R")
RSCRIPT_RESCUE_ML = os.path.join(utilitiesPath, "rescue","rescue_by_mapping_ML.R")
RESCUE_RANDOM_FOREST = os.path.join(utilitiesPath, "rescue","run_randomforest_on_reference.R")

# IsoAnnot
ISOANNOT_PROG =  os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")

# PYTHONPATH
PYTHONPATH = shutil.which('python')

def get_aligner_command(aligner_choice, genome, isoforms, annotation, 
                        outdir, corrSAM, n_cpu, gmap_index,logger=qc_logger):
    # Even though the speed does not change form the ifelse, this is cleaner
    match aligner_choice:
        case "gmap":
            logger.info("****Aligning reads with GMAP...")
            cmd = GMAP_CMD.format(
                cpus=n_cpu,
                dir=os.path.dirname(gmap_index),
                name=os.path.basename(gmap_index),
                i=isoforms,
                o=corrSAM,
            )
        case "minimap2":
            logger.info("****Aligning reads with Minimap2...")
            cmd = MINIMAP2_CMD.format(
                cpus=n_cpu,
                g=genome,
                i=isoforms,
                o=corrSAM,
            )
        case "deSALT":
            logger.info("****Aligning reads with deSALT...")
            cmd = DESALT_CMD.format(
                cpus=n_cpu,
                dir=gmap_index,
                i=isoforms,
                o=corrSAM,
            )
        case "uLTRA":
            logger.info("****Aligning reads with uLTRA...")
            cmd = ULTRA_CMD.format(
                cpus=n_cpu,
                prefix="../" + os.path.splitext(os.path.basename(corrSAM))[0],
                g=genome,
                a=annotation,
                i=isoforms,
                o_dir=outdir + "/uLTRA_out/",
            )
        case _:
            logger.error(f"Unsupported aligner choice: {aligner_choice}")
            raise ValueError()
    return cmd

# TransDecoder2
def run_td2(corrFASTA, orf_input,threads):
    """
    Runs the TD2 ORF prediction tool on the corrected FASTA file.
    """
    # Initial setup
    sqanti_path = os.path.dirname(corrFASTA)
    if orf_input is None:
        orf_input = corrFASTA
    qc_logger.info(f"Running TD2 ORF search on {orf_input}...")
    td2_path = os.path.join(sqanti_path, "TD2")

    # First we run the ORF search
    search_cmd = f"TD2.LongOrfs -t {corrFASTA} -O {td2_path} -S --threads {threads}"
    logFile = f"{sqanti_path}/logs/TD2_LongOrfs.log"
    run_command(search_cmd, qc_logger, logFile, description="TD2 ORF search")

    # Now the actual prediction
    predict_cmd = f"cd {td2_path}; TD2.Predict -t {corrFASTA} -O ./ "
    logFile = f"{sqanti_path}/logs/TD2_Predict.log"
    run_command(predict_cmd, qc_logger, logFile, description="TD2 ORF prediction")
    
    # get the name of the output file
    orf_output = os.path.join(td2_path, f"{os.path.basename(orf_input)}.TD2.pep")
    return orf_output

def GTF_to_genePred(corrGTF):
    """
    Converts a GTF file to a genePred file.
    """
    queryFile = os.path.splitext(corrGTF)[0] +".genePred"
    logFile = f"{os.path.dirname(corrGTF)}/logs/GTF_to_genePred.log"
    if os.path.exists(queryFile):
        qc_logger.info(f"{queryFile} already exists. Using it.")
    else:
        # gtf to genePred
        cmd = f"{GTF2GENEPRED_PROG} {corrGTF} {queryFile} -genePredExt -allErrors -ignoreGroupsWithoutExons"
        run_command(cmd,qc_logger,logFile, "GTF to genePred conversion")
    return queryFile

def run_command(cmd,logger,out_file='log/program.log',description="command execution",silent=True, fail_ok=False):
    """
    Executes a shell command and handles errors gracefully.
    
    :param cmd: The command to execute (string).
    :param description: A short description of the operation for better error messages (default: "command execution").
    :raises SystemExit: Exits the script if the command fails.
    """
    try:
        logger.debug(out_file)
        os.makedirs(os.path.dirname(out_file), exist_ok=True) # Just in case

        MAIN_LOGGING_CONFIG['handlers']['process_handler']['filename'] = out_file
        logging.config.dictConfig(MAIN_LOGGING_CONFIG)
        process_logger = logging.getLogger('process_logger')
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=silent,
            check=True,
            encoding="utf-8",
        )
        # On success, log both STDOUT and STDERR so that tools like Rscript
        # (which emit informative messages to STDERR) are not lost.
        if result.stdout:
            process_logger.info(result.stdout)
        if result.stderr:
            process_logger.info(result.stderr)

        logger.debug(f"Process {cmd} had {result.returncode}")
        logger.debug(f"Returncode {result.returncode} class is {type(result.returncode)}")
        
        if result.returncode !=0:
            logger.debug("EROOOOOOR")
            raise BrokenPipeError(result.stderr)
    except subprocess.CalledProcessError as e:
        process_logger.error(f"Details: {e}")
        process_logger.info(e.stdout)
        process_logger.error(e.stderr)
        logger.error(f"Something went wrong during {description}")
        logger.error(f"For more inflo, check {out_file}")
        if fail_ok:
            logger.warning(f"Continuing despite failure in {description}.")
            return -1
        else:
            sys.exit(1)
    except BrokenPipeError as e:
        if fail_ok:
            logger.warning(f"Continuing despite failure in {description}: {e}")
            return -1
        else:
            sys.exit(1)
    
