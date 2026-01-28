import pandas as pd
import numpy as np

def merge_classifications(reference_rules_path, classif, strategy, thr = 0.7):
    """Load input data files."""
    classif = classif.copy()
    columns = ["isoform", "filter_result"]
    if strategy == "ml":
        class_ref = pd.read_csv(reference_rules_path, sep="\t")
        # Create a filter_result column based on the threshold with Isoform/Artifact respectively
        class_ref["filter_result"] = np.where(
            class_ref["POS_MLprob"] >= thr, "Isoform", "Artifact"
        )
        
        columns += ["POS_MLprob"]
    else: # rules
        class_ref = pd.read_csv(reference_rules_path, sep="\t", usecols=columns)
    # Add the reference classification to the transcript classification
    class_ref["origin"] = "reference"
    classif.loc[:,"origin"] = "lr_defined"
    combined_class = pd.concat([class_ref, classif[columns + ["origin"]]])

    return  combined_class

def add_filter_results(mapping_hits, combined_class, classif):
    """Add the filter results to the mapping hits."""
    
    # Add filter result of mapping hits
    mapping_hits = mapping_hits.merge(
        combined_class.rename(columns={"isoform": "mapping_hit"}), 
        on="mapping_hit", 
        how="left"
    ).rename(columns={"filter_result": "hit_filter_result",
                      "POS_MLprob": "hit_POS_MLprob",
                      "origin": "hit_origin"})

    # Add structural categories of candidates
    mapping_hits = mapping_hits.merge(
        classif[["isoform", "structural_category","associated_transcript"]], 
        left_on="rescue_candidate", 
        right_on="isoform", 
        how="left"
    ).rename(columns={"structural_category": "candidate_structural_category"})
    mapping_hits.drop(columns=["isoform"], inplace=True)
    return mapping_hits
    
def select_best_hits(df):
    # Keep rows with max score per group
    df = df[df['score'] == df['score'].max()]
    # Prefer lr_defined over reference
    if (df['hit_origin'] == 'lr_defined').any():
        df = df[df['hit_origin'] == 'lr_defined']
    # For now, we are keeping all best hits (in case of ties) 
    return df

def filter_mapping_hits(mapping_hits, strategy):
    """Filter mapping hits and remove redundant references."""
    # 1. Filter the mapping hits reference based on the results of filtering
    mapping_hits_filt = mapping_hits[mapping_hits["hit_filter_result"] == "Isoform"].copy()
    if strategy == "rules":
        # Apply rules-based filtering
        mapping_hits_filt['score'] = mapping_hits_filt['alignment_score']

    elif strategy == "ml": #TODO: make this weighted?
        mapping_hits_filt['score'] = mapping_hits_filt['alignment_score'] * mapping_hits_filt['hit_POS_MLprob']

    # Select the best hit(s) per rescue candidate
    best_mapping_hits = mapping_hits_filt.groupby(
    'rescue_candidate', 
    group_keys=False
    )[mapping_hits_filt.columns].apply(select_best_hits).copy()

    return best_mapping_hits

def find_reintroduced_transcripts(mapping_hits_filt, automatic_rescue_array):
    """Find reintroduced transcripts after filtering."""
    full_rescue_array = mapping_hits_filt.loc[
        (mapping_hits_filt["hit_origin"] == "reference"), "mapping_hit"
    ].unique()
    # We decided to reintroduce a reference transcript even though there is an FSM that covers it
    if automatic_rescue_array.empty:
        combined = pd.Series(full_rescue_array)
    else:
        combined = pd.Series(np.concatenate([automatic_rescue_array, full_rescue_array]))
    return combined.drop_duplicates().reset_index(drop=True)

def merge_rescue_modes(auto_rescue_df, mapping_rescue_df, full_inclusion_list, strategy):
    """
    Merge rescue modes into a single DataFrame efficiently.
    """
    # 1. Convert lookup lists to Sets for O(1) speed
    # Checking "x in list" is slow (O(N)). Checking "x in set" is constant (O(1)).
    valid_transcripts = set(full_inclusion_list)
    existing_transcripts = set(auto_rescue_df["assigned_transcript"])

    # 2. Rename columns
    # We create a copy implicitly here, preventing side-effects on the original DF
    mapping_clean = mapping_rescue_df.rename(columns={
        "rescue_candidate": "artifact",
        "mapping_hit": "assigned_transcript",
        "hit_origin": "origin"
    })
    
    mapping_clean["rescue_mode"] = f"{strategy}_mapping"

    # 3. Calculate "Reintroduced" using Vectorized Boolean Masks
    # A. Is it in the valid list?
    is_valid = mapping_clean["assigned_transcript"].isin(valid_transcripts)
    # B. Is it NEW? (Not present in the auto_rescue dataframe)
    is_new = ~mapping_clean["assigned_transcript"].isin(existing_transcripts)
    # C. Is it the FIRST occurrence in this file?
    # (Ensures we only flag the transcript as 'reintroduced' once, avoiding double counting)
    is_first_occurrence = ~mapping_clean.duplicated(subset="assigned_transcript", keep="first")

    # Combine conditions
    mapping_clean["reintroduced"] = (is_valid & is_new & is_first_occurrence).map({True: "yes", False: "no"})

    # 4. Filter and Concatenate
    cols = ["artifact", "assigned_transcript", "rescue_mode", "origin", "reintroduced"]
    rescue_df = pd.concat([auto_rescue_df, mapping_clean[cols]], ignore_index=True)
    
    return rescue_df

def rescue_by_mapping(mapping_hits, reference_rules_path, classif, automatic_rescue_list,
                      rescue_df, strategy, thr=0.7):
    # Load data
    combined_class = merge_classifications(reference_rules_path,classif,strategy,thr)

    # Join rules
    mapping_hits = add_filter_results(mapping_hits, combined_class, classif)

    mapping_rescue_df = filter_mapping_hits(mapping_hits, strategy)
    
    full_inclusion_list = find_reintroduced_transcripts(mapping_rescue_df, automatic_rescue_list)

    rescue_df = merge_rescue_modes(rescue_df, mapping_rescue_df, full_inclusion_list, strategy)

    return full_inclusion_list, rescue_df

