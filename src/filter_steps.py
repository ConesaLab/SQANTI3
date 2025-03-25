import os,sys
from Bio import SeqIO

from src.utilities.cupcake.io.BioReaders import GMAPSAMReader
from src.utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

from src.module_logging import filter_logger
from src.commands import run_command, RSCRIPTPATH,RSCRIPT_ML,RSCRIPT_RULES,RSCRIPT_FILTER_REPORT,utilitiesPath

def filter_isoforms(filename,prefix,ids_to_keep):
    fafq_type = 'fasta'
    with open(filename) as h:
            if h.readline().startswith('@'): fafq_type = 'fastq'
    fout=open(prefix + '.filtered.' + fafq_type, 'w')
    for r in SeqIO.parse(open(filename), fafq_type):
        if r.id in ids_to_keep:
            SeqIO.write(r, fout, fafq_type)
    fout.close()
    filter_logger.info(f"Output written to: {fout.name}")

def filter_gtf(filename,prefix,ids_to_keep):
    outputGTF = prefix + '.filtered.gtf'
    with open(outputGTF, 'w') as f:
        for r in collapseGFFReader(filename):
            if r.seqid in ids_to_keep:
                write_collapseGFF_format(f, r)
        filter_logger.info(f"Output written to: {f.name}")
    f.close()

def filter_sam(filename,prefix,ids_to_keep):
    outputSAM = prefix + '.filtered.sam'
    with open(outputSAM, 'w') as f:
        for r in GMAPSAMReader(filename):
            if r.qID in ids_to_keep:
                f.write(r)
        filter_logger.info(f"Output written to: {f.name}")
    f.close()

def filter_gff3(filename,prefix,inclusion_f):
    outputGFF3 = prefix + '.filtered.gff3'
    awk_cmd = """awk 'FNR==NR {{ a[$1]; next }} ($1 in a)' {l} {g} > {o}""".format(l=inclusion_f, g=filename, o=outputGFF3)
    logFile = os.path.join(os.path.dirname(prefix), 'logs', 'filter_gff3.log')
    run_command(awk_cmd,filter_logger,logFile,"Filtering GFF3 file")
    filter_logger.info(f"Output written to: {outputGFF3}")

def filter_faa(filename,prefix,ids_to_keep):
    outputFAA = prefix + '.filtered.faa'
    with open(outputFAA, 'w') as f:
        for r in SeqIO.parse(open(filename), 'fasta'):
            if r.id in ids_to_keep:
                f.write(">{0}\n{1}\n".format(r.description, r.seq))
        filter_logger.info(f"Output written to: {f.name}")
    f.close()

def filter_files(args,outdir,prefix, ids_to_keep, inclusion_f):
    # TODO: Update this prefix so it dinamically detects if it is a corrected file or not
    prefix = f"{outdir}/{prefix}" 
    # filter FASTA/FASTQ file
    if args.filter_isoforms is not None:
        filter_isoforms(args.filter_isoforms, prefix, ids_to_keep)

    # filter GTF
    if args.filter_gtf is not None:
        filter_gtf(args.filter_gtf, prefix, ids_to_keep)

    # filter SAM
    if args.filter_sam is not None:
        filter_sam(args.filter_sam, prefix, ids_to_keep)

    # filter FAA
    if args.filter_faa is not None:
        filter_faa(args.filter_faa, prefix, ids_to_keep)

    # filter isoAnnot GFF3
    if args.isoAnnotGFF3 is not None:
        filter_gff3(args.isoAnnotGFF3, prefix, inclusion_f)


def prepare_ml_cmd(sq_class,prefix,outdir,percent_training,threshold,intrapriming,force_fsm_in,filter_mono_exonic,intermediate_files,max_class_size,TP=None,TN=None,remove_columns=None):
    cmd = f"{RSCRIPTPATH} {RSCRIPT_ML} -c {sq_class} -o {prefix} -d {outdir} -t {percent_training} \
        -j {threshold} -i {intrapriming} -f {force_fsm_in} -e {filter_mono_exonic} -m {intermediate_files} \
        -z {max_class_size}"

    if TP is not None:
        cmd + f" -p {TP}"
    if TN is not None:
        cmd + f" -n {TN}"
    if remove_columns is not None:
        cmd + f" -r {remove_columns}"
    return cmd

# TODO: These two functions are almost the same. There should be some way around to merge them ;)
def run_ML(args):
    cmd = prepare_ml_cmd(args.sqanti_class,args.output,args.dir,args.percent_training,args.threshold,args.intrapriming,
                        args.force_fsm_in,args.filter_mono_exonic,args.intermediate_files,args.max_class_size,args.TP,
                        args.TN,args.remove_columns)
    
    logFile = os.path.join(args.dir, 'logs', 'filter_ml.log')
    run_command(cmd,filter_logger,logFile,"Machine learning filtering")

    if not args.skip_report:
      report_cmd=f"{RSCRIPTPATH} {RSCRIPT_FILTER_REPORT} -d {args.dir} -o {args.output} -u {utilitiesPath} -f ml "
      logFile = os.path.join(args.dir, 'logs', 'filter_report.log')
      run_command(report_cmd,filter_logger,logFile,"Machine learning filtering report")

    # After running ML R code, an inclusion list will be generated. Those IDs must be passed to the filter files function
    inclusion_list = os.path.join(args.dir,"{args.output}_inclusion-list.txt")
    seqs_to_keep = set(line.strip() for line in open(inclusion_list))
    return(seqs_to_keep, inclusion_list)

def run_rules(args):
    cmd = f"{RSCRIPTPATH} {RSCRIPT_RULES} -c {args.sqanti_class} -o {args.output} -d {args.dir} -j {args.json_filter} -u {utilitiesPath} -e {args.filter_mono_exonic}"

    report_cmd = f"{RSCRIPTPATH} {RSCRIPT_FILTER_REPORT} -d {args.dir} -o {args.output} -u {utilitiesPath} -f rules"

    logFile = os.path.join(args.dir, 'logs', 'filter_rules.log')
    run_command(cmd,filter_logger,logFile,"Rules filtering")
    if not args.skip_report:
      logFile = os.path.join(args.dir, 'logs', 'filter_report.log')
      run_command(report_cmd,filter_logger,logFile,"Rules filtering report")

    # After running Rules Filter code, an inclusion list will be generated. Those IDs must be passed to the filter files function
    inclusion_list = os.path.join(args.dir,f"{args.output}_inclusion-list.txt")
    seqs_to_keep = set(line.strip() for line in open(inclusion_list))
    return(seqs_to_keep, inclusion_list)
