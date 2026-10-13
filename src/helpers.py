import os
import re
import gzip


from collections import defaultdict
from typing import Dict, Optional
from Bio import SeqIO #type: ignore
from Bio.SeqIO.FastaIO import SimpleFastaParser #type: ignore
from Bio.SeqIO.QualityIO import FastqGeneralIterator #type: ignore

from src.utilities.cupcake.sequence.err_correct_w_genome import err_correct
from src.utilities.cupcake.sequence.sam_to_gff3 import convert_sam_to_gff3

from src.config import seqid_rex1, seqid_rex2, seqid_fusion
from src.commands import get_aligner_command, GFFREAD_PROG, run_command, run_td2
from src.parsers import parse_TD2, parse_corrORF
from src.module_logging import qc_logger

### Environment manipulation functions ###
def clean_isoform_id(seqid):
    """
    Keep the part of a sequence ID before '|' or the first space.
    For example, "PB.1.1|chr1:10-100|xxxxxx" becomes "PB.1.1".
    """
    return seqid.split('|')[0].split()[0]

def is_fastq_file(input_fasta):
    """
    :param input_fasta: fasta or fastq. Can be gzipped.
    :return: True if the file is FASTQ (first line starts with '@'), False if FASTA
    """
    open_function = gzip.open if input_fasta.endswith('.gz') else open
    with open_function(input_fasta, mode="rt") as h:
        return h.readline().startswith('@')

def read_fasta_fastq(input_fasta):
    """
    Iterate over a FASTA or FASTQ file (autodetected, can be gzipped) without building SeqRecords,
    which is much faster for large read files. FASTQ qualities are not returned.

    :return: generator of (title, sequence); title is the full header line without '>' or '@'
    """
    open_function = gzip.open if input_fasta.endswith('.gz') else open
    is_fastq = is_fastq_file(input_fasta)
    with open_function(input_fasta, mode="rt") as h:
        if is_fastq:
            for title, seq, _ in FastqGeneralIterator(h):
                yield title, seq
        else:
            yield from SimpleFastaParser(h)

def find_isoform_id_collisions(input_fasta):
    """
    Find sequences that end up with the same ID after clean_isoform_id (or that already share an ID).
    Duplicated IDs give duplicated isoforms downstream (issue #216), so they must be detected before
    aligning.

    :param input_fasta: fasta or fastq, autodetected. Can be gzipped.
    :return: dict of (clean ID --> list of original IDs), only for IDs shared by more than one sequence
    """
    originals = defaultdict(list)
    for title, _ in read_fasta_fastq(input_fasta):
        originals[clean_isoform_id(title)].append(title.split()[0])
    return {new_id: ids for new_id, ids in originals.items() if len(ids) > 1}

def rename_isoform_seqids(input_fasta, out_dir):
    """
    Rename input isoform fasta/fastq by extracting the first part of the sequence ID.

    Handles various ID formats by taking the content before '|' or space characters.
    For example:
    - "PB.1.1|chr1:10-100|xxxxxx" becomes "PB.1.1"
    - "transcript_name some_annotation" becomes "transcript_name"

    :param input_fasta: Could be either fasta or fastq, autodetect. Can be gzipped.
    :param out_dir: directory for the output, so nothing is written next to the input (it may be read-only)
    :return: output fasta (<out_dir>/<input name>.renamed.fasta) with the cleaned up sequence IDs
    """
    in_name = os.path.basename(input_fasta)
    if in_name.endswith('.gz'):
        in_name = in_name[:-3]
    out_file = os.path.join(out_dir, os.path.splitext(in_name)[0] + '.renamed.fasta')

    with open(out_file, mode='wt') as f:
        for title, seq in read_fasta_fastq(input_fasta):
            f.write(f">{clean_isoform_id(title)}\n{seq}\n")
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

def get_supplementary_name(outdir, prefix):
    corrPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    return corrPathPrefix + "_supplementary_alignments.txt"

def report_non_primary_alignments(sam_file, out_file):
    """
    SQANTI3 only uses the primary alignment of each isoform. Isoforms that also have supplementary
    (split/chimeric) or secondary alignments are reported here, with all their alignments, so that
    the information is not lost (issue #216).

    :return: number of isoforms with non-primary alignments. out_file is only written if > 0.
    """
    def records():
        with open(sam_file) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                raw = line.split('\t', 6)
                if raw[2] == '*':
                    continue
                yield raw

    split_ids = {raw[0] for raw in records() if int(raw[1]) & (256 | 2048)}
    if not split_ids:
        return 0

    with open(out_file, 'w') as out:
        out.write("isoform\talignment\tchrom\tstrand\tstart\tend\tMAPQ\n")
        for qid, flag, chrom, pos, mapq, cigar, _ in records():
            if qid not in split_ids:
                continue
            flag = int(flag)
            kind = "supplementary" if flag & 2048 else "secondary" if flag & 256 else "primary"
            ref_len = sum(int(n) for n, op in re.findall(r'(\d+)([MDN=X])', cigar))
            start = int(pos)
            out.write(f"{qid}\t{kind}\t{chrom}\t{'-' if flag & 16 else '+'}\t{start}\t{start + ref_len - 1}\t{mapq}\n")
    return len(split_ids)

def get_unmapped_name(outdir, prefix):
    corrPathPrefix = os.path.abspath(os.path.join(outdir, prefix))
    return corrPathPrefix + "_unmapped.txt"

def report_unmapped_isoforms(isoforms_fasta, sam_file, out_file):
    """
    Isoforms without a primary alignment are dropped from all SQANTI3 outputs. They are listed here.
    The input FASTA is compared against the SAM instead of looking for unmapped records (flag 4),
    because not every aligner writes unmapped sequences to the SAM.

    :return: number of unmapped isoforms. out_file is only written if > 0.
    """
    mapped_ids = set()
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            qid, flag, chrom = line.split('\t', 3)[:3]
            if chrom != '*' and not int(flag) & (4 | 256 | 2048):
                mapped_ids.add(qid)

    unmapped = [(title.split()[0], len(seq)) for title, seq in read_fasta_fastq(isoforms_fasta)
                if title.split()[0] not in mapped_ids]
    if not unmapped:
        return 0

    with open(out_file, 'w') as out:
        out.write("isoform\tlength\n")
        for qid, length in unmapped:
            out.write(f"{qid}\t{length}\n")
    return len(unmapped)

def warn_supplementary(n_split, supplementary_file, in_chunk=False):
    """
    Log the isoforms with supplementary/secondary alignments. Same message with and without chunks:
    inside a chunk it is only a debug line, since the chunk file is deleted after combine_alignment_outputs
    merges it and logs the warning with the final path.
    """
    if n_split == 0:
        return
    if in_chunk:
        qc_logger.debug(f"Chunk: {n_split} isoforms have supplementary or secondary alignments.")
    else:
        qc_logger.warning(f"{n_split} isoforms have supplementary or secondary alignments (e.g. chimeric "
                          f"or split sequences). Only their primary alignment is used. "
                          f"All their alignments are listed in {supplementary_file}")

def warn_unmapped(n_unmapped, unmapped_file, in_chunk=False):
    """
    Log the isoforms that could not be aligned. See warn_supplementary for the chunk behaviour.
    """
    if n_unmapped == 0:
        return
    if in_chunk:
        qc_logger.debug(f"Chunk: {n_unmapped} isoforms could not be aligned to the genome.")
    else:
        qc_logger.warning(f"{n_unmapped} isoforms could not be aligned to the genome and are not included "
                          f"in the SQANTI3 output. They are listed in {unmapped_file}")

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
            isoforms = rename_isoform_seqids(isoforms, outdir)
            qc_logger.debug(f"Cleaned up isoform fasta file written to: {isoforms}")
    
            if os.path.exists(corrSAM):
                qc_logger.info(f"Aligned SAM {corrSAM} already exists. Using it...")
            else:
                logFile = f"{os.path.dirname(corrSAM)}/logs/{aligner_choice}_alignment.log"
                cmd = get_aligner_command(aligner_choice, genome, isoforms, annotation, 
                                          outdir,corrSAM, n_cpu, gmap_index)
                run_command(cmd,qc_logger, logFile,description="aligning reads")

            # Only primary alignments are used downstream; report the rest (issue #216)
            # Inside a chunk of a parallel run (chunks > 1), the warnings are logged after combining the chunks
            supplementary_file = get_supplementary_name(outdir, output)
            n_split = report_non_primary_alignments(corrSAM, supplementary_file)
            warn_supplementary(n_split, supplementary_file, in_chunk=chunks > 1)

            # Isoforms that did not align are dropped from all outputs; report them
            unmapped_file = get_unmapped_name(outdir, output)
            n_unmapped = report_unmapped_isoforms(isoforms, corrSAM, unmapped_file)
            warn_unmapped(n_unmapped, unmapped_file, in_chunk=chunks > 1)

            # The renamed copy of the input is only needed until here, and it is rebuilt on every run.
            # Remove it: for SQANTI-reads it is a full copy of the reads.
            os.remove(isoforms)

            # error correct the genome (input: corrSAM, output: corrFASTA)
            err_correct(genome, corrSAM, corrFASTA, genome_dict=genome_dict)
            qc_logger.debug(f"The corrected fasta file has been written to: {corrFASTA}")
            # convert SAM to GFF --> GTF
            n_skipped = convert_sam_to_gff3(corrSAM, f'{corrGTF}.tmp', source=os.path.basename(genome).split('.')[0])  # convert SAM to GFF3
            qc_logger.info(f"Skipped {n_skipped} unmapped SAM records when converting the alignments to GTF.")
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
        - Lines with 'transcript/exon/mRNA/match/cDNA_match/match_part/CDS' feature types are further processed:
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

    if feature_type in ('transcript', 'mRNA', 'exon', 'CDS', 'match', 'match_part', 'cDNA_match'):
        # Convert GFF3 alignment types and mRNA to standard GTF transcript/exon
        if feature_type in ('mRNA', 'match', 'cDNA_match'):
            fields[2] = 'transcript'
        elif feature_type == 'match_part':
            fields[2] = 'exon'
            
        line = "\t".join(fields) + "\n"
        
        if strand not in ['-', '+']:
            logger.warning(f"Discarding unknown strand feature: {line.strip()}")
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