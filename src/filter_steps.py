import os

from src.module_logging import filter_logger
from src.commands import (
    run_command, RSCRIPTPATH,RSCRIPT_ML,
    RSCRIPT_FILTER_REPORT,utilitiesPath
)
from src.filter_output import (
    filter_isoforms, filter_gtf, filter_sam, filter_gff3, filter_faa
)
from src.utilities.filter.sqanti3_rules_filter import rules_filter

def filter_files(isoforms,gtf,sam,faa,isoAnnotGFF3,
                 outdir,prefix, ids_to_keep, inclusion_f):
    # TODO: Update this prefix so it dinamically detects if it is a corrected file or not
    prefix = f"{outdir}/{prefix}" 
    # filter FASTA/FASTQ file
    if isoforms is not None:
        filter_isoforms(isoforms, prefix, ids_to_keep)

    # filter GTF
    if gtf is not None:
        filter_gtf(gtf, prefix, ids_to_keep)

    # filter SAM
    if sam is not None:
        filter_sam(sam, prefix, ids_to_keep)

    # filter FAA
    if faa is not None:
        filter_faa(faa, prefix, ids_to_keep)

    # filter isoAnnot GFF3
    if isoAnnotGFF3 is not None:
        filter_gff3(isoAnnotGFF3, prefix, inclusion_f)

RSCRIPT_RULES = os.path.join(utilitiesPath,"filter","SQANTI3_rules_filter.R")

def run_rules(args):

    prefix = os.path.join(args.dir, args.output)
    rules_filter(args.sqanti_class,args.json_filter,args.filter_mono_exonic,
                 prefix,filter_logger)
    if not args.skip_report:
      report_cmd = f"{RSCRIPTPATH} {RSCRIPT_FILTER_REPORT} -d {args.dir} -o {args.output} -u {utilitiesPath} -f rules"
      logFile = os.path.join(args.dir, 'logs', 'filter_report.log')
      run_command(report_cmd,filter_logger,logFile,"Rules filtering report")

    # After running Rules Filter code, an inclusion list will be generated. Those IDs must be passed to the filter files function
    inclusion_list = os.path.join(args.dir,f"{args.output}_inclusion-list.txt")
    seqs_to_keep = set(line.strip() for line in open(inclusion_list))
    return(seqs_to_keep, inclusion_list)


def prepare_ml_cmd(sq_class,prefix,outdir,percent_training,threshold,intrapriming,force_fsm_in,filter_mono_exonic,intermediate_files,max_class_size,TP=None,TN=None,remove_columns=None):
    cmd = f"{RSCRIPTPATH} {RSCRIPT_ML} -c {sq_class} -o {prefix} -d {outdir} -t {percent_training} \
        -j {threshold} -i {intrapriming} -f {force_fsm_in} -e {filter_mono_exonic} -m {intermediate_files} \
        -z {max_class_size}"

    if TP is not None:
        cmd += f" --TP {TP}"
    if TN is not None:
        cmd += f" --TN {TN}"
    if remove_columns is not None:
        cmd += f" --remove_columns {remove_columns}"
    return cmd

# TODO: These two functions are almost the same. There should be some way around to merge them ;)
def run_ML(args):
    cmd = prepare_ml_cmd(args.sqanti_class,args.output,args.dir,args.percent_training,args.threshold,args.intrapriming,
                        args.force_fsm_in,args.filter_mono_exonic,args.intermediate_files,args.max_class_size,args.TP,
                        args.TN,args.remove_columns)
    
    logFile = os.path.join(args.dir, 'logs', 'filter_ml.log')
    run_command(cmd,filter_logger,logFile,"Machine learning filtering", silent=False)

    if not args.skip_report:
      report_cmd=f"{RSCRIPTPATH} {RSCRIPT_FILTER_REPORT} -d {args.dir} -o {args.output} -u {utilitiesPath} -f ml "
      logFile = os.path.join(args.dir, 'logs', 'filter_report.log')
      run_command(report_cmd,filter_logger,logFile,"Machine learning filtering report")

    # After running ML R code, an inclusion list will be generated. Those IDs must be passed to the filter files function
    inclusion_list = os.path.join(args.dir,f"{args.output}_inclusion-list.txt")
    seqs_to_keep = set(line.strip() for line in open(inclusion_list))
    return(seqs_to_keep, inclusion_list)

