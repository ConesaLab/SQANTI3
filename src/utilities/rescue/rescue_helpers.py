import pandas as pd
from gtfparse import read_gtf
from src.module_logging import rescue_logger

def read_classification(filename):
    return pd.read_csv(filename, sep="\t")

def get_rescue_gene_targets(df, rescue_candidates):
    """Get the genes with associated rescue candidates.
    
    Args:
        classification_file (str): Path to SQANTI3 classification file
        rescue_candidates (pd.Series): List of rescue candidates
        
    Returns:
        pd.Series: List of genes with associated rescue candidates
    """
    return df[df['isoform'].isin(rescue_candidates)]['associated_gene'].unique()
   
def get_rescue_reference_targets(ref_gtf, target_genes):
    """
    Extracts transcript IDs for specific target genes from a GTF file.
    """
    # 1. Load data
    # Assuming read_gtf returns an object with a .to_pandas() method
    df = read_gtf(ref_gtf).to_pandas()

    # 2. Optimize lookup
    # Converting list to set allows O(1) lookups, though pandas .isin handles this well.
    target_genes_set = set(target_genes)
    
    # 3. Vectorized Filtering
    # We apply both filters at once:
    # - Rows that are transcripts
    # - Rows where gene_id is in our target list
    mask = (df['feature'] == 'transcript') & (df['gene_id'].isin(target_genes_set))
    
    # 4. Direct Selection
    # Select only the column we need. 
    # .unique() ensures we don't return duplicates if the GTF is malformed.
    target_transcripts = df.loc[mask, 'transcript_id'].dropna().unique()
    
    # Logging
    # Note: We access target_genes_set for the log to show unique gene count requested
    rescue_logger.debug(f"Input: {len(target_genes_set)} target genes.")
    rescue_logger.debug(f"Output: Found {len(target_transcripts)} matching transcripts in GTF.")

    return pd.Series(target_transcripts, name='isoform')