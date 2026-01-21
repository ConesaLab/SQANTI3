import pandas as pd
import warnings
warnings.filterwarnings("ignore") #TODO: See why pandas raises the warning on copying

from collections import defaultdict

from src.utilities.rescue.requant_helpers import (
    build_artifact_table, export_counts, redistribute_counts, calculate_tpm
)

def run_requant(counts, rescue_df, classif_df, prefix):
    """
    Redistribute counts from artifacts to rescued isoforms.
    Handles multiple sample columns.
    
    Args:
        counts: DataFrame with 'transcript_id' column and one or more count columns
        rescue_df: DataFrame with rescue results (artifact -> assigned_transcript mappings)
        classif_df: DataFrame with SQANTI3 classification (Isoform vs Artifact)
        prefix: Output file prefix
    
    Returns:
        DataFrame with reassigned counts (artifacts removed, counts redistributed)
    """
    # Creation of a table with all artifacts (rescued and not rescued)
    artifacts_df = build_artifact_table(rescue_df, classif_df)
    
    # Redistribute counts based on the rescue results
    new_counts = redistribute_counts(artifacts_df, classif_df, counts)
    
    return export_counts(counts, new_counts, prefix)

def to_tpm(counts_df, class_df, prefix):

    class_df.rename(columns={'isoform': 'transcript_id'}, inplace=True)
    class_df = class_df[["transcript_id","length"]]
    
    counts_d = defaultdict(int)
    counts_df.apply(lambda x: counts_d.update({x.iloc[0] : x.iloc[1]}), axis = 1) 
    class_df['counts'] = class_df['transcript_id'].apply(lambda x: counts_d[x]) # Apparently this line triggers a SettingWithCopyWarning
    #remove zero values
    class_df = class_df[class_df['counts'] != 0]
    # Calculate TPM
    class_df['TPM'] = calculate_tpm(class_df['counts'], class_df['length'])
    final = class_df[['transcript_id', 'TPM']]
    final.to_csv(f"{prefix}_reassigned_counts_TPM.tsv", header = True, sep='\t', index=False)