import pandas as pd
import numpy as np

from collections import defaultdict

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

def build_artifact_table(rescue_df, classif_df):
    """Combine rescued and non-rescued artifacts into one table."""
    not_rescued = get_unrescued_artifacts(classif_df, rescue_df)
    return pd.concat([rescue_df, not_rescued], ignore_index=True).rename(columns={'artifact': 'isoform'})


def redistribute_counts(artifacts_df, classif_df, old_counts):
    """Reassign counts from artifacts to their rescued isoforms.
    - Multi-transcript artifacts are grouped into a single placeholder entry.
    - Ensures count conservation while keeping multi-mapped cases separate.
    """
    new_counts = defaultdict(int)
    true_isoforms = classif_df[classif_df['filter_result'] == 'Isoform']['isoform'].tolist()

    # Retain valid isoforms’ counts
    for iso in true_isoforms:
        new_counts[iso] = old_counts.get(iso, 0)

    new_counts["multi_transcript_artifact"] = 0

    for _, row in artifacts_df.iterrows():
        iso = row['artifact']
        assigned = row['assigned_transcript']
        artifact_count = old_counts.get(iso, 0)

        if isinstance(assigned, list): # TODO: How to better handle this casuistic?
            new_counts["multi_transcript_artifact"] += artifact_count
        else: # Single isoform match
            new_counts[assigned] += artifact_count

    return new_counts

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
    # DEBUG
    # import pickle
    # with open("trial_issues.pkl", "wb") as f:
    #     pickle.dump((exploded_map, base_df, sample_cols), f)
    # f.close()
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

def export_counts(old_counts_df, new_counts_dict, prefix):
    """Save full, changed, and summarized count tables.
    Handles multiple sample columns and excludes artifacts from final output."""
    # Get sample columns
    sample_cols = [col for col in old_counts_df.columns if col != 'isoform']
    
    # Convert old counts to dict for comparison
    old_counts_dict = old_counts_df.set_index('isoform')[sample_cols].to_dict('index')
    
    # Get all transcript IDs (union of old and new)
    all_transcripts = list(set(old_counts_dict.keys()) | set(new_counts_dict.keys()))
    
    # Build extended comparison table
    extended_data = []
    changed_transcripts = []
    
    for transcript_id in all_transcripts:
        row = {'transcript_id': transcript_id}
        
        # Add old counts
        old_vals = old_counts_dict.get(transcript_id, {col: 0 for col in sample_cols})
        for sample in sample_cols:
            row[f'old_{sample}'] = old_vals.get(sample, 0)
        
        # Add new counts
        new_vals = new_counts_dict.get(transcript_id, {col: 0 for col in sample_cols})
        for sample in sample_cols:
            row[f'new_{sample}'] = new_vals.get(sample, 0)
        
        # Check if any sample changed
        if any(old_vals.get(s, 0) != new_vals.get(s, 0) for s in sample_cols):
            changed_transcripts.append(transcript_id)
        
        extended_data.append(row)
    
    # Create DataFrames
    extended_df = pd.DataFrame(extended_data)
    changed_df = pd.DataFrame({'changed_transcript': changed_transcripts})
    
    # Create final reassigned counts (only new counts, excluding artifacts)
    final_data = []
    for transcript_id, counts in new_counts_dict.items():
        row = {'transcript_id': transcript_id}
        row.update(counts)
        final_data.append(row)
    
    final_df = pd.DataFrame(final_data)
    
    # Sort by transcript_id for consistency
    extended_df = extended_df.sort_values('transcript_id').reset_index(drop=True)
    final_df = final_df.sort_values('transcript_id').reset_index(drop=True)
    
    # Save files
    extended_df.to_csv(f"{prefix}_reassigned_counts_extended.tsv", sep='\t', index=False)
    changed_df.to_csv(f"{prefix}_changed_counts.tsv", sep='\t', index=False)
    final_df.to_csv(f"{prefix}_reassigned_counts.tsv", sep='\t', index=False)
    
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