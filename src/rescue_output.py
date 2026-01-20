import shutil
from Bio import SeqIO

# GTF output
# TODO: Perhaps fix the order of the - strand exons
def write_rescue_gtf(input_gtf, ref_gtf, inclusion_list, prefix):
    output_gtf = f"{prefix}_rescued.gtf"
    shutil.copy(input_gtf, output_gtf)
    
    # Optimization 1: Convert list to set for O(1) lookups (huge speedup for large lists)
    inclusion_set = set(inclusion_list)

    with open(ref_gtf, 'r') as infile, open(output_gtf, 'a+') as outfile:
        outfile.seek(0, 2) # Move to the very end of the file
        if outfile.tell() > 0: # Check if file is not empty
            outfile.seek(outfile.tell() - 1, 0) # Move back one character
            last_char = outfile.read(1)
            if last_char != '\n':
                outfile.write('\n')
        
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            
            # Robustness: Skip malformed lines that don't have enough columns
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type not in ['transcript', 'exon']:
                continue
            
            attributes = fields[8]
            
            # Optimization 2: Simpler string search is faster than full parsing
            # We don't need to parse the whole dictionary just to check existence.
            # We check if the transcript_id string exists in the attributes.
            # NOTE: This is a heuristic. For 100% precision, use the parsing code below.
            
            # --- Robust Parsing (Fixes potential crash with spaces in values) ---
            found = False
            for attr in attributes.split(';'):
                attr = attr.strip()
                if not attr: 
                    continue
                
                # Split only on the first space to handle values like "Gene Name" safely
                parts = attr.split(' ', 1) 
                if len(parts) == 2:
                    key, value = parts
                    if key == 'transcript_id':
                        # Clean quotes
                        clean_tid = value.strip('"')
                        if clean_tid in inclusion_set:
                            found = True
                        break # Stop parsing attributes once we found the ID
            
            if found:
                outfile.write(line)
    return output_gtf

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
    
    # Use a set for O(1) lookups if filter_ids is a large list
    filter_ids_set = set(filter_ids)

    # List comprehension is faster and cleaner a the manual loop
    return [record for record in sequences_raw if record.id in filter_ids_set]

def write_rescue_fasta(input_fasta, ref_fasta_file, good_transcripts, inclusion_list, prefix):
    """
    Merge the rescued transcripts with the good isoforms and write to a new file.
    """
    # Load reference FASTA and filter
    ref_fasta = read_fasta_file(ref_fasta_file, inclusion_list)
    tr_fasta = read_fasta_file(input_fasta, good_transcripts)

    write_fasta_file(f"{prefix}_rescued.fasta", tr_fasta + ref_fasta)

    