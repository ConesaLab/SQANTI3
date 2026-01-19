import pandas as pd
from gtfparse import read_gtf
from src.module_logging import rescue_logger
from src.utilities.rescue.automatic_rescue import get_lost_reference_id

def read_classification(filename):
    return pd.read_csv(filename, sep="\t")


def identify_rescue_candidates(classif_df, monoexons='all'):
    """
    Pure Logic: Filters the classification dataframe to find rescue candidates.
    
    Args:
        classif_df (pd.DataFrame): The classification data.
        monoexons (str): 'all' to keep everything, otherwise filters single exons.
        
    Returns:
        pd.DataFrame: A subset of classif_df containing only the candidates.
    """
    # 1. Helper Logic: Identify transcripts that already have a match (FSM)
    # Using a set for O(1) lookups
    transcripts_with_fsm = set(
        classif_df.loc[classif_df['structural_category'] == 'full-splice_match', 'associated_transcript']
    )

    # 2. Select NIC and NNC artifacts
    # These are candidates regardless of FSM status
    nic_nnc_mask = (
        (classif_df['structural_category'].isin(['novel_in_catalog', 'novel_not_in_catalog'])) &
        (classif_df['filter_result'] == 'Artifact')
    )
    nic_nnc_candidates = classif_df[nic_nnc_mask]

    # 3. Select ISM candidates
    # Only select ISMs if their reference transcript was "lost" (not found in FSMs)
    lost_ids = get_lost_reference_id(classif_df)
    
    ism_mask = (
        (classif_df['structural_category'] == 'incomplete-splice_match') &
        (classif_df['associated_transcript'].isin(lost_ids)) & 
        (~classif_df['associated_transcript'].isin(transcripts_with_fsm))
    )
    ism_candidates = classif_df[ism_mask]

    # 4. Merge and Filter Monoexons
    final_candidates = pd.concat([nic_nnc_candidates, ism_candidates])

    if monoexons != 'all':
        final_candidates = final_candidates[final_candidates['exons'] > 1]
        
    return final_candidates

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