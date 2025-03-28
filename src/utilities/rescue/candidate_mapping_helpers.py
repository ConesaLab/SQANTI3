import os, sys

from Bio import SeqIO

from src.module_logging import rescue_logger
from src.commands import run_command



def prepare_fasta_transcriptome(ref_gtf,ref_fasta,outdir):
    rescue_logger.info("Creating reference transcriptome FASTA from provided GTF (--refGTF).")

    # make FASTA file name
    pre, _ = os.path.splitext(os.path.basename(ref_gtf))
    ref_trans_Fasta = os.path.join(outdir,f"{pre}.fasta")

    # build gffread command
    ref_cmd = f"gffread -w {ref_trans_Fasta} -g {ref_fasta} {ref_gtf}"

  # run gffread
    logFile=os.path.join(outdir,"logs","rescue","create_reference_transcriptome.log")
    run_command(ref_cmd,rescue_logger,logFile,description="Converting reference transcriptome GTF to FASTA")
    rescue_logger.debug(f"File created in {ref_trans_Fasta}.")
    if os.path.isfile(ref_trans_Fasta):
        rescue_logger.info(f"Reference transcriptome FASTA was saved to {ref_trans_Fasta}")
        rescue_logger.info("gffread command used:")
        rescue_logger.info(ref_cmd)
    else:
        rescue_logger.error("Reference transcriptome FASTA was not created - file not found!")
        sys.exit(1)
    return ref_trans_Fasta

def filter_transcriptome(input_fasta, target_ids, output_fasta):
 
    # Filter and write sequences
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.id in target_ids:
                SeqIO.write(record, out_f, 'fasta')