import shutil
from Bio import SeqIO

# GTF output
# TODO: Perhaps fix the order of the - strand exons
def write_rescue_gtf(input_gtf, ref_gtf, inclusion_list, prefix):
    output_gtf = f"{prefix}_rescued.gtf"
    shutil.copy(input_gtf, output_gtf)
    with open(ref_gtf, 'r') as infile, open(output_gtf, 'a') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.strip().split('\t')
            feature_type = fields[2]
            if feature_type not in ['transcript', 'exon']:
                continue
            attributes = fields[8]
            attr_dict = {
                key: value.strip('"') 
                for key, value in (attr.strip().split() for attr in attributes.split(';') if attr.strip())
            }
            if 'transcript_id' in attr_dict and attr_dict['transcript_id'] in inclusion_list:
                outfile.write(line)

    outfile.close()

# FASTA output
def write_fasta_file(fasta_file, sequences):
    """Save sequences to a FASTA file."""
    with open(fasta_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")

def read_fasta_file(fasta_file, filter_ids=None):
    """Read sequences from a FASTA file, optionally filtering by IDs."""
    sequences_raw = SeqIO.parse(fasta_file, "fasta")
    if filter_ids is None:
        return list(sequences_raw)
    else:
        sequences = []
        for record in sequences_raw:
            if filter_ids is None or record.id in filter_ids:
                sequences.append(record)
        return sequences

def write_rescue_fasta(input_fasta, ref_fasta_file, good_transcripts, inclusion_list, prefix):
    """
    Merge the rescued transcripts with the good isoforms and write to a new file.
    """
    # Load reference FASTA and filter
    ref_fasta = read_fasta_file(ref_fasta_file, inclusion_list)
    tr_fasta = read_fasta_file(input_fasta, good_transcripts)

    write_fasta_file(f"{prefix}_rescued.fasta", tr_fasta + ref_fasta)

    