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
    target_genes = df[df['isoform'].isin(rescue_candidates)]['associated_gene'].unique()

    return target_genes

def parse_rescue_gtf(gtf_file):
    """Parse a GTF file to return a dictionary of gene_id and associated transcript_id.
    
    Args:
        gtf_filename (str): Path to GTF file
        
    Returns:
        Dictionary of gene_id and transcript_id
    """
    # Load the GTF file into a DataFrame
    df = read_gtf(gtf_file).to_pandas()
   
    # 1. Filter rows AND select only needed columns in one step to reduce memory footprint
    # 2. Use .copy() to avoid SettingWithCopy warnings later
    subset = df.loc[df['feature'] == 'transcript', ['gene_id', 'transcript_id']].copy()
    
    # 3. Drop rows with missing values (robustness)
    subset = subset.dropna()
    
    # 4. Group by gene_id. 
    # using .agg(list) is idiomatic and often slightly optimized over .apply(list)
    return subset.groupby('gene_id')['transcript_id'].agg(list).to_dict()

def get_rescue_reference_targets(ref_gtf, target_genes):

    gene_to_transcript = parse_rescue_gtf(ref_gtf)
    rescue_logger.debug(f"Found {len(gene_to_transcript)} genes in reference GTF.")
    rescue_logger.debug(f"There are {len(target_genes)} target genes.")
    target_transcripts = []
    for gene,transcripts in gene_to_transcript.items():
        if gene in target_genes:
            target_transcripts.extend(transcripts)
    rescue_logger.debug(f"Found {len(target_transcripts)} target transcripts.")
    return pd.Series(target_transcripts, name='isoform')
        
    # Get the genes with associated rescue candidates
def get_good_transcripts(class_file):
    """
    Get the transcripts that pass the SQANTI3 filter.
    """
    df = pd.read_csv(class_file, sep='\t', comment='#')
    return list(df[df['filter_result'] == 'Isoform']['isoform'])