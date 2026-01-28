import pandas as pd
import numpy as np

def get_unrescued_artifacts(classif_df, rescue_df):
    """Find artifacts that remain unrescued or were not linked to any rescued transcript."""
    not_rescued = (
        classif_df[
            (classif_df["filter_result"] == "Artifact") &
            ~(classif_df["isoform"].isin(rescue_df["artifact"]))
        ][["isoform", "associated_gene"]]
        .rename(columns={"isoform": "artifact"})
    )
    # Assign transcript divergency for unrescued artifacts
    # If the associated gene is not novel, use gene_TD; else, use general_TD
    not_rescued["assigned_transcript"] = not_rescued["associated_gene"].apply(
        lambda gene: f"{gene}_TD" if not str(gene).startswith('novel') else "general_TD"
    )
    # remove the associated_gene column as it's no longer needed
    not_rescued.drop(columns=["associated_gene"], inplace=True)
 
    return not_rescued

def prepare_count_matrices(old_counts, classid_df):
    """
    Splits the original count table into two matrices:
    1. base_df: Valid isoforms (Destinations)
    2. source_df: Artifacts (Sources to be redistributed)
    """
    # Standardize index
    counts_df = old_counts.set_index('isoform') if 'isoform' in old_counts.columns else old_counts
    
    # Identify IDs
    valid_ids = classid_df[classid_df['filter_result'] == 'Isoform']['isoform']
    artifact_ids = classid_df[classid_df['filter_result'] == 'Artifact']['isoform']
    
    # Create aligned matrices (Fill NaN with 0 for safe math)
    base_df = counts_df.reindex(valid_ids).fillna(0).astype(int)
    source_df = counts_df.reindex(artifact_ids).fillna(0).astype(int)
    
    return base_df, source_df

def calculate_distribution_fractions(exploded_map, base_df, sample_cols):
    """
    Calculates the proportion of the artifact count that should go to each target.
    Logic: 
      - Proportional to target's abundance in base_df.
      - If all targets have 0 counts, split evenly (Uniform).
    """

    # 1. Align Target Counts to the Map
    # Join the counts of the 'assigned_transcript' onto the exploded rows
    target_weights = exploded_map.merge(
        base_df, 
        left_on='assigned_transcript', 
        right_index=True, 
        how='left'
    )[sample_cols].fillna(0)

    # Ensure index alignment with exploded_map for grouping
    target_weights.index = exploded_map.index

    # 2. Calculate Totals per Artifact Group
    grouped = target_weights.groupby(exploded_map['isoform'])
    total_weights = grouped.transform('sum')

    # 3. Handle Edge Cases (Total Weight = 0)
    zero_mask = (total_weights == 0)
    group_sizes = grouped.transform('count')

    # 4. Compute Fractions
    fractions = pd.DataFrame(0.0, index=target_weights.index, columns=sample_cols)
    
    # Case A: Proportional Split
    fractions[~zero_mask] = target_weights[~zero_mask] / total_weights[~zero_mask]
    
    # Case B: Uniform Split (1/N)
    fractions[zero_mask] = 1.0 / group_sizes[zero_mask]
    
    return fractions

def distribute_integer_counts(source_df, exploded_map, fractions, sample_cols):
    """
    Applies fractions to source counts and enforces integer conservation.
    Strategy: Floor(Value) + Remainder added to the first target.
    """
    # 1. Align Source (Artifact) counts to the exploded map
    # We look up the count for the artifact (isoform) corresponding to each row
    artifact_counts_aligned = source_df.loc[exploded_map['isoform']].reset_index(drop=True)[sample_cols]

    # 2. Basic Distribution (Floored to Integer)
    # e.g., 10 * 0.33 = 3.3 -> 3
    counts_to_add = np.floor(artifact_counts_aligned * fractions).astype(int)

    # 3. Calculate Remainder (Integer Conservation)
    # Check how many counts we lost due to flooring
    # e.g., Original: 10. Allocated Sum: 3+3+3=9. Remainder: 1.
    allocated_sum = counts_to_add.groupby(exploded_map['isoform']).transform('sum')
    remainders = (artifact_counts_aligned - allocated_sum).astype(int)

    # 4. Distribute Remainder
    # Add the remainder to the FIRST target of each artifact group to preserve the sum.
    # This avoids creating decimals or losing counts.
    is_first = (exploded_map.groupby('isoform').cumcount() == 0)
    
    # We only update rows where 'is_first' is True
    # Pandas aligns indices automatically, so this safely adds the remainder vector
    counts_to_add.loc[is_first] += remainders.loc[is_first]

    # 5. Assign Target IDs for final aggregation
    counts_to_add['isoform'] = exploded_map['assigned_transcript']
    counts_to_add.fillna(0, inplace=True)
    # 6. Sum up all contributions per Target
    final_additions = counts_to_add.groupby('isoform')[sample_cols].sum()
    
    return final_additions

def export_counts(old_counts_df, new_counts_df, prefix):
    """Save two count tables:
    1. Reassigned counts (new counts only, with 'isoform' column)
    2. Extended table for ALL isoforms (with old and new counts for each sample)
       - Good isoforms (kept their counts)
       - Artifacts (lost their counts: old > 0, new = 0)
       - Rescued/TD isoforms (gained counts: old = 0, new > 0)
    
    Handles multiple sample columns.
    """
   # Standardize column names if needed to ensure 'isoform' key exists
    if 'isoform' not in old_counts_df.columns:
        old_counts_df = old_counts_df.rename(columns={old_counts_df.columns[0]: 'isoform'})
    if 'isoform' not in new_counts_df.columns:
        new_counts_df = new_counts_df.rename(columns={'transcript_id': 'isoform'})

    # Get sample columns
    sample_cols = [col for col in old_counts_df.columns if col != 'isoform']
    
    # Create first dataframe: new counts with 'isoform' column + sample columns
    final_df = new_counts_df[['isoform'] + sample_cols].copy()
    final_df = final_df.sort_values('isoform').reset_index(drop=True)
    
    # Create second dataframe: extended table for ALL isoforms
    # Optimized: Use Outer Join to align Old and New counts automatically (avoids manual looping)
    extended_df = pd.merge(
        old_counts_df, 
        new_counts_df, 
        on='isoform', 
        how='outer', 
        suffixes=('_old', '_new')
    )
    
    # Fill missing values (0 if not present in one of the tables)
    # This covers cases where artifacts are missing in new, or rescued are missing in old
    extended_df = extended_df.fillna(0)
    
    # Rename columns to match desired format (old_Sample, new_Sample)
    # and reorder them for readability
    ordered_cols = ['isoform']
    rename_map = {}
    
    for sample in sample_cols:
        col_old = f"{sample}_old"
        col_new = f"{sample}_new"
        
        # Verify columns exist before renaming (robustness check)
        if col_old in extended_df.columns and col_new in extended_df.columns:
            rename_map[col_old] = f"old_{sample}"
            rename_map[col_new] = f"new_{sample}"
            ordered_cols.extend([f"old_{sample}", f"new_{sample}"])

    extended_df = extended_df.rename(columns=rename_map)
    
    # Sort and save both dataframes
    if not extended_df.empty:
        # Select specific columns to clean up any extra metadata and apply sorting
        extended_df = extended_df[ordered_cols].sort_values('isoform').reset_index(drop=True)
        
        # Ensure integer counts (Merge converts to float due to NaNs)
        numeric_cols = extended_df.columns.drop('isoform')
        extended_df[numeric_cols] = extended_df[numeric_cols].astype(int)
    
    # Save files
    final_df.to_csv(f"{prefix}_reassigned_counts.tsv", sep='\t', index=False)
    extended_df.to_csv(f"{prefix}_reassigned_counts_extended.tsv", sep='\t', index=False)
    
    return final_df

def calculate_tpm(counts, lengths):
    # Convert lengths to kilobases
    lengths_kb = lengths / 1000
    # Calculate RPK
    rpk = counts / lengths_kb
    # Calculate Total RPK
    total_rpk = rpk.sum()
    # Calculate TPM
    tpm = (rpk / total_rpk) * 1e6
    return tpm