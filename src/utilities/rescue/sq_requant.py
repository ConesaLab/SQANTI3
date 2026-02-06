import pandas as pd
import numpy as np
import warnings

from src.parsers import parse_counts
warnings.filterwarnings("ignore") #TODO: See why pandas raises the warning on copying

from collections import defaultdict

from src.module_logging import rescue_logger
from src.utilities.rescue.requant_helpers import (
    calculate_distribution_fractions, distribute_integer_counts, 
    export_counts, get_unrescued_artifacts, prepare_count_matrices, 
    calculate_tpm
)

def requantification_pipeline(output_dir, output_prefix, counts_file, rescue_df, original_class, rescue_class):
    """
    Executes the requantification pipeline for transcript abundance estimation.

    This function processes count data and rescue information to reassign transcript
    counts, calculate TPM values, and save the results to output files.

    Parameters
    ----------
    output_dir : str
        Directory path where output files will be saved.
    output_prefix : str
        Prefix string to use for naming output files.
    counts_file : str
        Path to the input counts file to be parsed.
    rescue_df : pandas.DataFrame
        DataFrame containing rescued transcript information.
    rescue_class : object
        Class or object containing rescue classification data and methods.

    Returns
    -------
    None
        The function writes results to files but does not return a value.

    Notes
    -----
    - Parses the counts file and logs the completion
    - Runs requantification on counts using rescue data
    - Saves reassigned counts to {output_prefix}_reassigned_counts.tsv
    - Converts counts to TPM (Transcripts Per Million) values
    - Multi-transcript and artifact counts are lost during TPM calculation
      as they lack length information required for TPM computation
    """
    prefix = f"{output_dir}/{output_prefix}"
    #TODO: Make this take the variables from python directly
    counts_df = parse_counts(counts_file)
    # Keep only counts present in the original classification
    counts_df = counts_df[counts_df['isoform'].isin(original_class['isoform'])]
    rescue_logger.info("Counts file parsed.")
    requant_df = requantify(counts_df, rescue_df, original_class, prefix)
    rescue_logger.info("Requantification of counts completed.")
    rescue_logger.info(f"New count table saved to {prefix}_reassigned_counts.tsv")
    update_classification(requant_df, rescue_class, prefix)
    # Doing this, we loose the counts assigned to multi_transcript and artifacts (they have no length, so TPM cannot be calculated)
    #to_tpm(requant_df,rescue_class, prefix)
    rescue_logger.info("Requantification finished!")
    return

def requantify(counts, rescue_df, classif_df, prefix):
    """
    Redistribute counts from artifacts to rescued isoforms.
    Handles multiple sample columns.
    
    Args:
        counts: DataFrame with 'isoform' column and one or more count columns
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

def build_artifact_table(rescue_df, classif_df):
    """Combine rescued and non-rescued artifacts into one table."""
    not_rescued = get_unrescued_artifacts(classif_df, rescue_df)
    return pd.concat([rescue_df, not_rescued], ignore_index=True).rename(columns={'artifact': 'isoform'})

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
    # 4. Distribute & Conserve Integers
    final_additions = distribute_integer_counts(source_df, rescue_df, fractions, sample_cols)

    # 5. Merge Result
    final_counts = base_df.add(final_additions, fill_value=0).astype(int)

    return final_counts.reset_index().rename(columns={'index': 'isoform'})

def update_classification(requant_df, rescue_class, prefix):
    """
    Update the FL count columns in the rescue classification based on reassigned counts.
    - Safe: Only updates isoforms present in BOTH files.
    - Prevents creation of 'ghost' rows (NaNs) if requant_df has extra IDs.
    """
    # 1. Set Indices
    req_df = requant_df.set_index('isoform')
    class_df = rescue_class.set_index('isoform').copy()
    
    # 2. Safety Filter: Intersection of indices
    # We only want to update rows that actually exist in the classification file.
    common_ids = req_df.index.intersection(class_df.index)
    
    # If no matches found, warn and return early
    if common_ids.empty:
        print("Warning: No matching isoforms found between requantification and classification files.")
        return class_df.reset_index()

    # Create a subset of the source data ensuring alignment
    req_subset = req_df.loc[common_ids]
    sample_cols = req_df.columns

    # 3. Update Individual Sample Columns
    for sample in sample_cols:
        target_col = f"FL.{sample}"
        
        # Initialize column if it doesn't exist (important to avoid NaNs for rows we don't update)
        if target_col not in class_df.columns:
            class_df[target_col] = 0
        
        # Update ONLY the common rows
        class_df.loc[common_ids, target_col] = req_subset[sample]

    # 4. Update Total 'FL' Column
    # Calculate sum only for the subset
    new_totals = req_subset.sum(axis=1)
    
    # Update 'FL', initializing it if needed
    if 'FL' not in class_df.columns:
        class_df['FL'] = 0
    class_df.loc[common_ids, 'FL'] = new_totals

    # 5. Save
    output_file = f"{prefix}_rescued_classification.txt"
    class_df.reset_index().to_csv(output_file, sep='\t', index=False)

def to_tpm(counts_df, class_df, prefix):

    class_df = class_df.rename(columns={'isoform': 'transcript_id'})
    class_df = class_df[["transcript_id","length"]].copy()
    
    counts_d = defaultdict(int)
    counts_df.apply(lambda x: counts_d.update({x.iloc[0] : x.iloc[1]}), axis = 1) 
    class_df['counts'] = class_df['transcript_id'].apply(lambda x: counts_d[x])
    #remove zero values
    class_df = class_df[class_df['counts'] != 0]
    # Calculate TPM
    class_df['TPM'] = calculate_tpm(class_df['counts'], class_df['length'])
    final = class_df[['transcript_id', 'TPM']]
    final.to_csv(f"{prefix}_reassigned_counts_TPM.tsv", header = True, sep='\t', index=False)