import pandas as pd

from collections import defaultdict

def reassign_counts(isoform, old_counts, new_counts, list_of_changed):
    if isoform in new_counts.keys():
        count = new_counts[isoform]
        list_of_changed.append(isoform)
    else:
        try:
            count = old_counts[isoform]
        except KeyError:
            count = 0
    return(count)

def fill_old_counts(isoform, old_counts):
    try:
        count = old_counts[isoform]
    except KeyError:
        count = 0
    return(count)

def collapse_rescued_transcripts(rescue_df):
    """Group artifacts by assigned transcripts.
    Using a list only when multiple transcripts exist avoids nested lists later."""
    collapsed = (
        rescue_df.groupby("artifact", as_index=False)
        .agg({"assigned_transcript": lambda x: x.tolist() if len(x) > 1 else x.iloc[0]})
    )
    collapsed["num_assigned_transcripts"] = collapsed["assigned_transcript"].apply(
        lambda x: len(x) if isinstance(x, list) else 1
    )
    return collapsed


def get_unrescued_artifacts(classif_df, collapsed_df):
    """Find artifacts that remain unrescued or were not linked to any rescued transcript."""
    not_rescued = (
        classif_df[
            (classif_df["filter_result"] == "Artifact") &
            ~(classif_df["isoform"].isin(collapsed_df["artifact"]))
        ][["isoform", "assigned_gene"]]
        .rename(columns={"isoform": "artifact"})
    )
    # Assign transcript divergency for unrescued artifacts
    # If there is an associated gene, assign it to '<gene_name>_TD; else, label as 'general_TD'
    not_rescued["assigned_transcript"] = not_rescued["assigned_gene"].apply(
        lambda gene: f"{gene}_TD" if pd.notna(gene) else "general_TD"
    )
    # remove the assigned_gene column as it's no longer needed
    not_rescued.drop(columns=["assigned_gene"], inplace=True)
    not_rescued["num_assigned_transcripts"] = 1
    
    return not_rescued


def build_artifact_table(rescue_df, classif_df):
    """Combine rescued and non-rescued artifacts into one table."""
    collapsed = collapse_rescued_transcripts(rescue_df)
    not_rescued = get_unrescued_artifacts(classif_df, collapsed)
    return pd.concat([collapsed, not_rescued], ignore_index=True)


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


def export_counts(old_counts_df, new_counts_dict, prefix):
    """Save full, changed, and summarized count tables.
    Handles multiple sample columns and excludes artifacts from final output."""
    # Get sample columns
    sample_cols = [col for col in old_counts_df.columns if col != 'transcript_id']
    
    # Convert old counts to dict for comparison
    old_counts_dict = old_counts_df.set_index('transcript_id')[sample_cols].to_dict('index')
    
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