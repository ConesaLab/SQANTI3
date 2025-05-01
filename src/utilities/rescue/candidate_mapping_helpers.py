import os, sys
import pysam
import pandas as pd

from Bio import SeqIO

from src.module_logging import rescue_logger
from src.commands import run_command



def prepare_fasta_transcriptome(ref_gtf,ref_fasta,outdir):
    """Creates a reference transcriptome from a genome and a gtf file

    Args:
        ref_gtf (str): Path to the GTF transcriptome
        ref_fasta (str): Path to the reference genome 
        outdir (str): Path to the output directory
    """
    rescue_logger.info("Creating reference transcriptome FASTA from provided GTF (--refGTF).")

    # make FASTA file name
    pre, _ = os.path.splitext(os.path.basename(ref_gtf))
    ref_trans_Fasta = os.path.join(outdir,f"{pre}.fasta")

    # build gffread command
    ref_cmd = f"gffread -w {ref_trans_Fasta} -g {ref_fasta} {ref_gtf}"

  # run gffread
    logFile=os.path.join(outdir,"logs","create_reference_transcriptome.log")
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

def filter_transcriptome(input_fasta, target_ids):
    """Reads a fasta file and select specific sequences

    Args:
        input_fasta (str): Path to the input FASTA file 
        target_ids (list): List of sequence IDs to be selected

    Returns:
        list: Selected sequences in SeqIO format
    """
    target_records = []
    # Filter and write sequences
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if record.id in target_ids:
            target_records.append(record)
    return target_records

def save_fasta(records, output_fasta):
    """Writes a SeqIO fasta object into a given file

    Args:
        records (SeqIO): SeqIO object containing the sequences to be written
        output_fasta (str): Path to the output FASTA file
    """
    with open(output_fasta, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

def process_sam_file(sam_file, output_dir, output_prefix):
    # Define file paths
    hits_file = f"{output_dir}/{output_prefix}_rescue_mapping_hits.tsv"

    # Open the SAM file and process it
    with pysam.AlignmentFile(sam_file, "r") as sam:
        # Extract candidate-target pairs and alignment type
        data = []
        for read in sam.fetch(until_eof=True):  # Skip header automatically
            data.append([read.query_name, read.reference_name, read.flag])

    # Convert to DataFrame and save as TSV
    hits_df = pd.DataFrame(data, columns=["candidate", "target", "alignment_type"])
    hits_df.to_csv(hits_file, sep="\t", index=False, header=False)

    rescue_logger.info(f"Mapping hit table was saved to {hits_file}")