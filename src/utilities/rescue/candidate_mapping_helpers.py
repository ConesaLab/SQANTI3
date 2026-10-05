import os, sys
import hashlib
import pysam
import pandas as pd

from Bio import SeqIO

from src.module_logging import rescue_logger
from src.commands import run_command

def prepare_fasta_transcriptome(ref_gtf,ref_fasta,outdir):
    rescue_logger.info("Creating reference transcriptome FASTA from provided GTF (--refGTF).")

    # make FASTA file name
    pre, _ = os.path.splitext(os.path.basename(ref_gtf))
    ref_trans_Fasta = os.path.join(outdir,f"{pre}.isoforms.fasta")

    # build gffread command
    ref_cmd = f"gffread -w {ref_trans_Fasta} -g {ref_fasta} {ref_gtf}"

  # run gffread
    logFile=os.path.join(outdir,"logs","create_reference_transcriptome.log")
    run_command(ref_cmd,rescue_logger,logFile,description="Converting reference transcriptome GTF to FASTA")
    rescue_logger.debug(f"File created in {ref_trans_Fasta}.")
    if os.path.isfile(ref_trans_Fasta):
        rescue_logger.info(f"Reference transcriptome FASTA was saved to {ref_trans_Fasta}")
    else:
        rescue_logger.error("Reference transcriptome FASTA was not created - file not found!")
        sys.exit(1)
    return ref_trans_Fasta

def filter_transcriptome(input_fasta, target_ids):
    target_records = []
    # Filter and write sequences
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if record.id in target_ids:
            target_records.append(record)
    return target_records

def merge_target_records(lr_records, ref_records):
    """
    Merge long read and reference target records, keeping only the first record per ID.
    Long read records take priority, so a reference transcript with the same ID (e.g. Bambu
    reports known transcripts with their reference ID) is dropped. Duplicates within each
    list are dropped too, as minimap2 would write them as duplicated @SQ header lines.
    Returns the merged records and the list of dropped IDs.
    """
    seen = set()
    merged, dropped = [], []
    for record in list(lr_records) + list(ref_records):
        if record.id in seen:
            dropped.append(record.id)
        else:
            seen.add(record.id)
            merged.append(record)
    return merged, dropped

def save_fasta(records, output_fasta):
    with open(output_fasta, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

def mapping_fingerprint(input_files, mapping_options):
    """
    Build a SHA-256 fingerprint of the mapping inputs: the content of each input file
    plus the mapping options. File paths are not included, so moving the output
    directory does not invalidate a previous mapping.
    """
    fingerprint = hashlib.sha256(mapping_options.encode())
    for input_file in input_files:
        file_hash = hashlib.sha256()
        with open(input_file, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                file_hash.update(chunk)
        fingerprint.update(file_hash.digest())
    return fingerprint.hexdigest()

def is_mapping_reusable(sam_file, fingerprint):
    """
    A previous SAM file is only reused if it has a fingerprint sidecar file
    matching the fingerprint of the current mapping inputs.
    """
    fingerprint_file = f"{sam_file}.sha256"
    if not (os.path.isfile(sam_file) and os.path.isfile(fingerprint_file)):
        return False
    with open(fingerprint_file) as fh:
        return fh.read().strip() == fingerprint

def save_mapping_fingerprint(sam_file, fingerprint):
    with open(f"{sam_file}.sha256", 'w') as fh:
        fh.write(f"{fingerprint}\n")

def tags_to_dict(tags):
    """Convert a list of (tag, value) tuples to a dictionary."""
    return {tag: value for tag, value in tags}

def process_sam_file(sam_file):

    # Open the SAM file and process it
    with pysam.AlignmentFile(sam_file, "r") as sam:
        # Extract candidate-target pairs and alignment type
        data = []
        for read in sam.fetch(until_eof=True):  # Skip header automatically
            try:
                data.append([read.query_name, read.reference_name, read.flag,tags_to_dict(read.tags)['AS']])
            except KeyError:
                data.append([read.query_name, read.reference_name, read.flag,0])

    return pd.DataFrame(data, columns=["rescue_candidate", "mapping_hit", "alignment_type","alignment_score"])