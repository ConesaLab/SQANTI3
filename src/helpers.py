import os
import gzip


from typing import Dict, Optional
from Bio import SeqIO #type: ignore

from src.utilities.cupcake.sequence.err_correct_w_genome import err_correct
from src.utilities.cupcake.sequence.sam_to_gff3 import convert_sam_to_gff3

from src.config import seqid_rex1, seqid_rex2, seqid_fusion
from src.commands import get_aligner_command, GFFREAD_PROG, run_command, run_td2
from src.parsers import parse_TD2, parse_corrORF
from src.module_logging import qc_logger

### Environment manipulation functions ###
def rename_isoform_seqids(input_fasta):
    """
    Rename input isoform fasta/fastq by extracting the first part of the sequence ID.
    
    Handles various ID formats by taking the content before '|' or space characters.
    For example:
    - "PB.1.1|chr1:10-100|xxxxxx" becomes "PB.1.1"
    - "transcript_name some_annotation" becomes "transcript_name"

    :param input_fasta: Could be either fasta or fastq, autodetect. Can be gzipped.
    :return: output fasta with the cleaned up sequence ID
    """
    type = 'fasta'
    # gzip.open and open have different default open modes:
    # gzip.open uses "rb" (read in binary format)
    # open uses "rt" (read in text format)
    # This can be solved by making explicit the read text mode (which is required
    # by SeqIO.parse)
    if input_fasta.endswith('.gz'):
        open_function = gzip.open
        in_file = os.path.splitext(input_fasta)[0]
        out_file = os.path.splitext(in_file)[0] + '.renamed.fasta'
    else:
        open_function = open
        out_file = os.path.splitext(input_fasta)[0] + '.renamed.fasta'
    with open_function(input_fasta, mode="rt") as h:
        if h.readline().startswith('@'): type = 'fastq'
        
    with open(out_file, mode='wt') as f:
        for r in SeqIO.parse(open_function(input_fasta, "rt"), type):
            # Extract the first part of the ID (before '|' or space)
            newid = r.id.split('|')[0].split()[0]
            f.write(">{0}\n{1}\n".format(newid, r.seq))
    return out_file

### Input/Output functions ###
def get_corr_filenames(outdir, prefix):
    corrPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    corrGTF = corrPathPrefix + "_corrected.gtf"
    corrSAM = corrPathPrefix + "_corrected.sam"
    corrFASTA = corrPathPrefix + "_corrected.fasta"
    corrORF = corrPathPrefix + "_corrected.faa"
    corrCDS_GTF_GFF = corrPathPrefix + "_corrected.cds.gff3"
    return corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF

def get_isoform_hits_name(outdir, prefix):
    corrPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    isoform_hits_name = corrPathPrefix + "_isoform_hits.txt"
    return isoform_hits_name

def get_class_junc_filenames(outdir, prefix):
    outputPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    outputClassPath = outputPathPrefix + "_classification.txt"
    outputJuncPath = outputPathPrefix + "_junctions.txt"
    return outputClassPath, outputJuncPath

def get_pickle_filename(outdir, prefix):
    pklPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    pklFilePath = pklPathPrefix + ".isoforms_info.pkl"
    return pklFilePath

def get_omitted_name(outdir, prefix):
    corrPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    omitted_name = corrPathPrefix + "_omitted_due_to_min_ref_len.txt"
    return omitted_name

def sequence_correction(
    outdir: str,
    output: str,
    cpus: int,
    chunks: int,
    fasta: bool,
    genome_dict: Dict[str, str],
    badstrandGTF: str,
    genome: str,
    isoforms: str,
    aligner_choice: str,
    gmap_index: Optional[str] = None,
    annotation: Optional[str] = None
    ) -> None:
    """
    Use the reference genome to correct the sequences (unless a pre-corrected GTF is given)
    """
    qc_logger.info("**** Correcting sequences")
    corrGTF, corrSAM, corrFASTA, _ , _ = get_corr_filenames(outdir, output)
    n_cpu = max(1, cpus // chunks)

    # Step 1. IF GFF or GTF is provided, make it into a genome-based fasta
    #         IF sequence is provided, align as SAM then correct with genome
    if os.path.exists(corrFASTA):
        qc_logger.info(f"Error corrected FASTA {corrFASTA} already exists. Using it...")
    else:
        qc_logger.info("Correcting fasta")
        if fasta:
            qc_logger.info("Cleaning up isoform IDs...")
            isoforms = rename_isoform_seqids(isoforms)
            qc_logger.info(f"Cleaned up isoform fasta file written to: {isoforms}")
    
            if os.path.exists(corrSAM):
                qc_logger.info(f"Aligned SAM {corrSAM} already exists. Using it...")
            else:
                logFile = f"{os.path.dirname(corrSAM)}/logs/{aligner_choice}_alignment.log"
                cmd = get_aligner_command(aligner_choice, genome, isoforms, annotation, 
                                          outdir,corrSAM, n_cpu, gmap_index)
                run_command(cmd,qc_logger, logFile,description="aligning reads")

            # error correct the genome (input: corrSAM, output: corrFASTA)
            err_correct(genome, corrSAM, corrFASTA, genome_dict=genome_dict)
            # convert SAM to GFF --> GTF
            convert_sam_to_gff3(corrSAM, f'{corrGTF}.tmp', source=os.path.basename(genome).split('.')[0])  # convert SAM to GFF3
        else:
            qc_logger.info("Skipping aligning of sequences because GTF file was provided.")
            filter_gtf(isoforms, f'{corrGTF}.tmp', badstrandGTF, genome_dict)
            if not os.path.exists(corrSAM):
                qc_logger.info("Indels will be not calculated since you ran SQANTI3 without alignment step (SQANTI3 with gtf format as transcriptome input).")

            # GTF to FASTA
            cmd = f"{GFFREAD_PROG} {corrGTF}.tmp -g {genome} -w {corrFASTA}"
            logFile = f"{outdir}/logs/gtf2fasta.log"
            run_command(cmd,qc_logger,logFile,description="Converting corrected GTF to FASTA")
        # Final step of converting the GFF3 to GTF or normalizing the GTF
        cmd = f"{GFFREAD_PROG} {corrGTF}.tmp -T -o {corrGTF}"
        logFile= f"{outdir}/logs/normalize_gtf.log"
        run_command(cmd,qc_logger,logFile, description="converting SAM to GTF")
        try:
            os.remove(f'{corrGTF}.tmp')
        except OSError as e:
            qc_logger.error(f"Error removing temporary file: {e}")
            raise

def filter_gtf(isoforms: str, corrGTF, badstrandGTF, genome_dict: Dict[str, str]) -> None:
    try:
        with open(corrGTF, 'w') as corrGTF_out, \
            open(isoforms, 'r') as isoforms_gtf, \
            open(badstrandGTF, 'w') as discard_gtf:
            for line in isoforms_gtf:
                process_gtf_line(line, genome_dict, corrGTF_out, discard_gtf)
    except IOError as e:
        qc_logger.error(f"Something went wrong processing GTF files: {e}")
        raise

def process_gtf_line(line: str, genome_dict: Dict[str, str], corrGTF_out: str, discard_gtf: str,logger=qc_logger):
    """
    Processes a single line from a GTF file, validating and categorizing it based on certain criteria.

    Args:
        line (str): A single line from a GTF file.
        genome_dict (Dict[str, str]): A dictionary containing genome reference data, where keys are chromosome names.
        corrGTF_out (str): Path to a file to write valid GTF lines with known strand information.
        discard_gtf (str): Path to a file to write GTF lines with unknown strand information.
    Raises:
        ValueError: If the chromosome in the GTF line is not found in the genome reference dictionary.

    Notes:
        - Lines starting with '#' are ignored.
        - Lines with fewer than 7 fields are considered malformed and skipped with a warning.
        - Lines with 'transcript' or 'exon' feature types are further processed:
            - If the strand is unknown ('-' or '+'), the line is written to the discard_gtf file with a warning.
            - Otherwise, the line is written to the corrGTF_out file.
    """
    if line.startswith("#"):
        return

    fields = line.strip().split("\t")
    if len(fields) < 7:
        logger.warning(f"Skipping malformed GTF line: {line.strip()}")
        return

    chrom, feature_type, strand = fields[0], fields[2], fields[6]

    if chrom not in genome_dict:
        logger.error(f"GTF chromosome {chrom} not found in genome reference file.")
        raise ValueError()

    if feature_type in ('transcript', 'exon'):
        if strand not in ['-', '+']:
            logger.warning(f"Discarding unknown strand transcript: {line.strip()}")
            discard_gtf.write(line)
        else:
            corrGTF_out.write(line)

def predictORF(outdir, include_ORF, orf_input, corrFASTA, corrORF, psauron_thr, threads):
    # ORF generation
    qc_logger.info("**** Predicting ORF sequences...")

    td2_dir = os.path.join(os.path.abspath(outdir), "TD2")
    if not os.path.exists(td2_dir):
        os.makedirs(td2_dir)

    # TD2 output --> myQueryProteins object
    cdsDict = {}
    if not include_ORF:
        qc_logger.warning("Skipping ORF prediction because user requested it. All isoforms will be non-coding!")
    elif os.path.exists(corrORF):
        qc_logger.info(f"ORF file {corrORF} already exists. Using it.")
        cdsDict = parse_corrORF(corrORF)
    else:
        td2_output = run_td2(corrFASTA, orf_input, psauron_thr, threads)  # threads is not used in TD2.Predict
        # Modifying ORF sequences by removing sequence before ATG
        cdsDict = parse_TD2(corrORF,td2_output)
    if len(cdsDict) == 0:
        qc_logger.warning("All input isoforms were predicted as non-coding")

    return(cdsDict)


def rename_novel_genes(isoform_info,novel_gene_prefix=None):
    """
    Rename novel genes to be "novel_X" where X is a number
    """
    novel_gene_index= 1
    for isoform_hit in isoform_info.values():
        if isoform_hit.structural_category in ("intergenic", "genic_intron"):
            # Liz: I don't find it necessary to cluster these novel genes. They should already be always non-overlapping.
            prefix = f'novelGene_{novel_gene_prefix}_' if novel_gene_prefix is not None else 'novelGene_'
            isoform_hit.genes = [f'{prefix}{novel_gene_index}']
            isoform_hit.transcripts = ['novel']
            novel_gene_index += 1
    return isoform_info