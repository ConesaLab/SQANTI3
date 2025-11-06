import pandas as pd
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

def parse_files(count_file):
    # Load counts file
    counts = pd.read_csv(count_file, sep = '\t', comment = '#')
    counts.columns = ['transcript_id', 'count']

    return(counts)

def select_hit(isoform, rescued, new_counts, old_counts):
    group = rescued[rescued["mapping_hit"] == isoform]
    group['counts'] = group.apply(lambda x: old_counts.get(x['rescue_candidate'], 0), axis = 1)
    new_counts[isoform] = group['counts'].sum()
    return new_counts

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

def run_requant(counts, rescue_df, classif_df, prefix):

    # Modify input table to have the artifact, the transcripts that are associated and how many
    collapsed_df = (
        rescue_df.groupby("artifact", as_index=False)
        .agg({
            "assigned_transcript": lambda x: x.tolist() if len(x) > 1 else x.iloc[0],
        })
    )

    collapsed_df["num_assigned_transcripts"] = collapsed_df["assigned_transcript"].apply(
        lambda x: len(x) if isinstance(x, list) else 1
    )

    # Select only isoforms that were not rescued or had a transcript that was already rescued
    
    not_rescued = (
        classif_df[
            (classif_df["filter_result"] == "Artifact") &
            ~(classif_df["isoform"].isin(collapsed_df["artifact"]))
        ][["isoform"]]
        .rename(columns={"isoform": "artifact"})
    )

    not_rescued["assigned_transcript"] = "artifact"
    not_rescued["num_assigned_transcripts"] = 1

    artifacts_df = pd.concat([collapsed_df, not_rescued], ignore_index=True)

    #create dictionaries of old and new counts
    old_counts = counts.set_index('transcript_id')['count'].to_dict()
    new_counts = defaultdict(int)

    # Add good counts
    true_isoforms = classif_df[classif_df['filter_result'] == 'Isoform']['isoform'].tolist()
    for isoform in true_isoforms:
        new_counts[isoform] = old_counts.get(isoform, 0)
    
    # Add artifact counts
    new_counts["multi_transcript_artifact"] = 0
    for isoform in artifacts_df['artifact']:
        new_counts[isoform] = 0
        assigned_tr = artifacts_df.loc[artifacts_df['artifact'] == isoform, 'assigned_transcript'].values[0]
        if isinstance(assigned_tr, list):  # How to handle this?
            new_counts["multi_transcript_artifact"] += old_counts.get(isoform, 0)
        else: # Single isoform match
            if assigned_tr in new_counts:
                new_counts[assigned_tr] += old_counts.get(isoform, 0)
            else:
                new_counts[assigned_tr] = old_counts.get(isoform, 0)
    
    
    #combine old and new counts
    counts_df = pd.DataFrame({'transcript_id': list(set(old_counts.keys()) | set(new_counts.keys()))})
    counts_df['old_count'] = counts_df['transcript_id'].apply(lambda x: fill_old_counts(x, old_counts))
    list_of_changed = []
    counts_df['new_count'] = counts_df['transcript_id'].apply(lambda x: reassign_counts(x, old_counts, new_counts, list_of_changed))
    changed = pd.DataFrame()
    changed['changed_count'] = list_of_changed
    counts_df.to_csv(f"{prefix}_reassigned_counts_extended.tsv", header = True, index = False, sep = '\t')
    changed.to_csv(f"{prefix}_changed_counts.tsv", header = True, index = False, sep = '\t')
    counts_df_short = counts_df[['transcript_id', 'new_count']]
    counts_df_short.to_csv(f"{prefix}_reassigned_counts.tsv", header = True, index = False, sep = '\t')
    return counts_df_short

def to_tpm(counts_df, class_df, prefix):
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

    class_df.rename(columns={'isoform': 'transcript_id'}, inplace=True)
    class_df = class_df[["transcript_id","length"]]
    
    counts_d = defaultdict(int)
    counts_df.apply(lambda x: counts_d.update({x.iloc[0] : x.iloc[1]}), axis = 1)
    class_df['counts'] = class_df['transcript_id'].apply(lambda x: counts_d[x])
    #remove zero values
    class_df = class_df[class_df['counts'] != 0]
    # Calculate TPM
    class_df['TPM'] = calculate_tpm(class_df['counts'], class_df['length'])
    final = class_df[['transcript_id', 'TPM']]
    final.to_csv(f"{prefix}_reassigned_counts_TPM.tsv", header = True, sep='\t', index=False)