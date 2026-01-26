import numpy as np
import warnings

from src.parsers import parse_counts
warnings.filterwarnings("ignore") #TODO: See why pandas raises the warning on copying

from collections import defaultdict

from src.module_logging import rescue_logger
from src.utilities.rescue.requant_helpers import (
    build_artifact_table, calculate_distribution_fractions, distribute_integer_counts, export_counts, prepare_count_matrices, calculate_tpm
)

def requantificaiton_pipeline(output_dir, output_prefix, counts_file, rescue_df, rescue_class):
    prefix = f"{output_dir}/{output_prefix}"
    #TODO: Make this take the variables from python directly
    counts_df = parse_counts(counts_file)
    rescue_logger.info("Counts file parsed.")
    requant_df = run_requant(counts_df, rescue_df, rescue_class, prefix)
    rescue_logger.info("Requantification of counts completed.")
    rescue_logger.info(f"New count table saved to {prefix}_reassigned_counts.tsv")
    # Doing this, we loose the counts assigned to multi_transcript and artifacts (they have no length, so TPM cannot be calculated)
    to_tpm(requant_df,rescue_class, prefix)
    rescue_logger.info("Requantification finished!")
    return

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
    artifacts_df = build_artifact_table(rescue_df[["artifact","assigned_transcript"]], classif_df)

    # Redistribute counts based on the rescue results
    new_counts = redistribute_counts_vectorized(artifacts_df, classif_df, counts)
    
    return export_counts(counts, new_counts, prefix)

def redistribute_counts_vectorized(rescue_df, classid_df, old_counts):
    """
    Main pipeline to reassign artifact counts.
    Now accepts the RAW (un-collapsed) rescue_df.
    """
    # 1. Setup
    if 'isoform' in old_counts.columns:
        temp_df = old_counts.set_index('isoform')
    else:
        temp_df = old_counts
    sample_cols = temp_df.select_dtypes(include=[np.number]).columns

    # 2. Prepare Matrices (Base vs Source)
    base_df, source_df = prepare_count_matrices(old_counts, classid_df)

    # 3. Calculate Weights
    fractions = calculate_distribution_fractions(rescue_df, base_df, sample_cols)
  
    fractions = fractions.fillna(0)
    # 5. Distribute & Conserve Integers
    final_additions = distribute_integer_counts(source_df, rescue_df, fractions, sample_cols)

    # 6. Merge Result
    final_counts = base_df.add(final_additions, fill_value=0).astype(int)

    return final_counts.reset_index().rename(columns={'index': 'isoform'})

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