import pandas as pd
import numpy as np
from src.module_logging import rescue_logger

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
    # If there is a tie between two reference hits, keep one 
    # if (df['hit_origin'] == 'reference').sum() > 1:
    #     df = df[df['hit_origin'] == 'reference'].iloc[[0]]


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

def merge_rescue_modes(auto_rescue_df, full_rescue_df,full_inclusion_list, strategy):
    """Merge rescue modes into a single DataFrame."""
    # Fix full rescue dataframe
    full_rescue_df = full_rescue_df.rename(columns={
        "rescue_candidate": "artifact",
        "mapping_hit": "assigned_transcript",
        "hit_origin":"origin"
    })
    full_rescue_df["rescue_mode"] = f"{strategy}_mapping"
    full_rescue_df["reintroduced"] = (
        full_rescue_df["assigned_transcript"].isin(full_inclusion_list)
    ).map({True: "yes", False: "no"})
    
    full_rescue_df = full_rescue_df[[
        "artifact", "assigned_transcript", "rescue_mode","origin", "reintroduced"]]
    # Combine automatic and full rescue dataframes
    rescue_df = pd.concat([auto_rescue_df, full_rescue_df])
    return rescue_df

def find_reintroduced_transcripts(mapping_hits_filt, automatic_rescue_array):
    """Find reintroduced transcripts after filtering."""
    full_rescue_array = mapping_hits_filt.loc[
        (mapping_hits_filt["hit_origin"] == "reference"), "mapping_hit"
    ].unique()
    # TODO: Should we not reintroduce a reference transcript if it is represented by an FSM?
    if automatic_rescue_array.empty:
        combined = pd.Series(full_rescue_array)
    else:
        combined = pd.Series(np.concatenate([automatic_rescue_array, full_rescue_array]))
    return combined.drop_duplicates().reset_index(drop=True)

def old_reintroduction(mapping_hits_filt, class_ref, classif, automatic_rescue_df):

    #DEBUG: mapping_hits_filt.to_csv("debug_mapping_hits_filt.tsv", sep="\t", index=False)
    # 2. Select only reference rescued transcripts
    rescued_ref = mapping_hits_filt[mapping_hits_filt["origin"] == "reference"]
 
    # 3. Remove reference transcripts already in the transcriptome
    ## 3.1 Retrieve all reference transcripts already represented by isoforms
    isoform_assoc_tr = classif[
        (classif["filter_result"] == "Isoform") & 
        (classif["associated_transcript"] != "novel")
    ][["associated_transcript"]]
    
    ## 3.2 Include those retrieved in automatic rescue
    isoform_assoc_tr = pd.concat([isoform_assoc_tr, automatic_rescue_df["associated_transcript"]]).drop_duplicates()
    # Find truly rescued references
    rescued_mapping_final = rescued_ref[
        ~rescued_ref["mapping_hit"].isin(isoform_assoc_tr["associated_transcript"])
    ][["mapping_hit"]].rename(columns={"mapping_hit": "ref_transcript"}).drop_duplicates()
    
    # Generate final list of rescued transcripts
    automatic_ref_rescued = automatic_rescue_df.rename(columns={"associated_transcript": "ref_transcript"})["ref_transcript"]
    reintroduced_final = pd.concat([automatic_ref_rescued, rescued_mapping_final]).drop_duplicates()
    return reintroduced_final, rescued_mapping_final, mapping_hits_filt, rescued_ref, isoform_assoc_tr

def write_rescue_inclusion_list(rescued_final,output_prefix):
    """Write the rescue inclusion list to a file."""
    output_path = f"{output_prefix}_full_inclusion_list.tsv"
    rescued_final.to_csv(output_path, sep="\t", index=False, header=False)
    rescue_logger.debug(f"Wrote rescue inclusion list to {output_path}")

def write_rescue_full_table(best_mapping_hits, output_prefix):
    """Write the full rescue mapping hits table to a file."""
    output_path = f"{output_prefix}_full_rescue_mapping_hits.tsv"
    best_mapping_hits.to_csv(output_path, sep="\t", index=False)
    rescue_logger.debug(f"Wrote full rescue mapping hits to {output_path}")

def process_automatic_rescue(classif, automatic_ref_rescued, class_ref):
    """Process automatic rescue results."""
    if automatic_ref_rescued.empty:
        automatic_fsm = pd.DataFrame(columns=["isoform","associated_transcript","structural_category","associated_gene"])
    else:
        automatic_fsm = classif[(classif["filter_result"] == "Artifact") &
                                    (classif["associated_transcript"].isin(automatic_ref_rescued["associated_transcript"]))
                                ][["isoform","associated_transcript","structural_category","associated_gene"]]

    automatic_fsm = automatic_fsm.merge(
            class_ref,
            left_on="associated_transcript",
            right_on="isoform",
            how="left"
        )
    
    automatic_fsm.drop(columns=["isoform_y"], inplace=True)
    automatic_fsm = automatic_fsm.rename(columns={"isoform_x":"rescue_candidate",
                                  "associated_transcript":"mapping_hit",
                                  "POS_MLprob":"hit_POS_MLprob",
                                  "filter_result":"hit_filter_result",
                                  "structural_category":"candidate_structural_category"})
    automatic_fsm["sam_flag"] = None
    automatic_fsm["rescue_result"] = "rescued_automatic"
    automatic_fsm["exclusion_reason"] = None
    return automatic_fsm

def get_exclusion_reason(row,mapping_hits_max,probs_ref,rescued_ref,isoform_assoc_tr,strategy):
    """Get the reason for exclusion of a mapping hit."""
    if row['mapping_hit'] not in mapping_hits_max['mapping_hit'].values:
        if strategy == "rules":
            return 'artifact_by_rules'
        elif strategy == "ml":
            return 'ML_probability'
    elif (row['mapping_hit'] in mapping_hits_max['mapping_hit'].values and 
            row['mapping_hit'] not in probs_ref['isoform'].values):
        return 'long_read_transcript'
    elif (row['mapping_hit'] in rescued_ref['mapping_hit'].values and
            row['mapping_hit'] in isoform_assoc_tr['associated_transcript'].values):
        return 'reference_already_present'
    return None

def assign_best_match(group):
    if ((group['exclusion_reason'] == "reference_already_present").any() or
        (group['rescue_result'] == "rescued_mapping").any() or
        (group['rescue_result'] == "rescued_automatic").any()):
        group['best_match_for_candidate'] = "reference_transcript"
        
    elif ((group['rescue_result'] == "not_rescued").all() and
          (group['exclusion_reason'] != "reference_already_present").all() and
          (group['exclusion_reason'] == "long_read_transcript").any()):
        group['best_match_for_candidate'] = "long_read_transcript"
        
    elif (group['exclusion_reason'] == "ML_probability").all():
        group['best_match_for_candidate'] = "unknown"

    elif (group['exclusion_reason'] == "artifact_by_rules").all():
        group['best_match_for_candidate'] = "unknown"
        
    else:
        group['best_match_for_candidate'] = np.nan  # fallback if none match
    
    return group

def assign_best_match_id(rescue_table,strategy,threshold):
    # 1. Filter by ML probability threshold and keep only max ML prob within each group
    if strategy == "rules":
        rescue_table_max = (
            rescue_table[rescue_table["hit_filter_result"] == "Isoform"]
            .loc[lambda df: df.groupby("rescue_candidate")["hit_filter_result"].transform("max") == df["hit_filter_result"]]
        )
    elif strategy == "ml":

        rescue_table_max = (
            rescue_table[rescue_table["hit_POS_MLprob"] >= threshold]
            .loc[lambda df: df.groupby("rescue_candidate")["hit_POS_MLprob"].transform("max") == df["hit_POS_MLprob"]]
        )

# 2. Unique best match per candidate (only 1 row in group)
    match_unique_ids = (
        rescue_table_max.groupby("rescue_candidate")
        .filter(lambda g: len(g) == 1)[["rescue_candidate", "mapping_hit"]]
        .rename(columns={"mapping_hit": "best_match_id"})
    )

    # 3. Ambiguous (tie) cases
    rescue_table_ties = (
        rescue_table_max.groupby("rescue_candidate")
        .filter(lambda g: len(g) > 1)
    )

    # 4. Split ties into primary and non-primary alignment groups
    rescue_table_ties_prim = rescue_table_ties[rescue_table_ties["sam_flag"] == 0]

    # Candidates without primary alignments
    nonprim_candidates = rescue_table_ties.loc[
        ~rescue_table_ties["rescue_candidate"].isin(rescue_table_ties_prim["rescue_candidate"])
    ]

    # 5. Get best match IDs for tie groups
    match_tie_ids_prim = (
        rescue_table_ties_prim.groupby("rescue_candidate")["mapping_hit"]
        .apply(lambda x: ",".join(x))
        .reset_index(name="best_match_id")
    )

    match_tie_ids_nonprim = (
        nonprim_candidates.groupby("rescue_candidate")["mapping_hit"]
        .apply(lambda x: ",".join(x))
        .reset_index(name="best_match_id")
    )

    # Combine all match IDs
    match_ids = pd.concat([match_tie_ids_prim, match_tie_ids_nonprim, match_unique_ids], ignore_index=True)

    # 6. Join back to rescue_table
    rescue_table = rescue_table.merge(match_ids, on="rescue_candidate", how="left")

    # 7. Fill in best_match_id for special cases
    rescue_table["best_match_id"] = np.where(
        rescue_table["best_match_for_candidate"] == "unknown",
        "unknown",
        rescue_table["best_match_id"]
    )

    rescue_table["best_match_id"] = np.where(
        rescue_table["rescue_result"] == "rescued_automatic",
        rescue_table["mapping_hit"],
        rescue_table["best_match_id"]
    )
    return(rescue_table)

def create_rescue_table(mapping_hits, best_mapping_hits, rescued_mapping_final, class_ref,
                        rescued_ref, isoform_assoc_tr, automatic_fsm, strategy, threshold):
    """Create the rescue table.
    
    Input:
        mapping_hits: DataFrame with mapping hits
        best_mapping_hits: DataFrame with the best mapping hits for each candidate
        rescued_mapping_final: DataFrame with the final rescued mapping
        class_ref: DataFrame with the results of filtering the reference
        isoform_assoc_tr: DataFrame with the isoform association
        automatic_fsm: DataFrame with the automatic rescue results
    """
    rescue_table = mapping_hits.copy()
    # Add result of rescue
    rescue_table["rescue_result"] = rescue_table["mapping_hit"].apply(
        lambda x: "rescued_mapping" if x in rescued_mapping_final["ref_transcript"].values else "not_rescued"
    )
    # Add the reason of why a transcript was excluded
    rescue_table['exclusion_reason'] = rescue_table.apply(get_exclusion_reason, axis=1, args=(best_mapping_hits, class_ref,
                                                                                              rescued_ref, isoform_assoc_tr, strategy))

    rescue_table = pd.concat([rescue_table, automatic_fsm], ignore_index=True)

    # Add the best match for each candidate
    rescue_table = rescue_table.groupby('rescue_candidate', group_keys=False).apply(
        assign_best_match
    )

    rescue_table = assign_best_match_id(rescue_table,strategy,threshold)

    return rescue_table

def write_rescue_table(rescue_table, output_prefix):
    """Write the rescue table to a file."""
    # Reorder columns
    cols = [col for col in rescue_table.columns if col != "associated_gene"] + ["associated_gene"]
    rescue_table = rescue_table[cols]
    output_path = f"{output_prefix}_full_rescue_table.tsv"
    rescue_table.to_csv(output_path, sep="\t", index=False)

def rescue_by_mapping(mapping_hits, reference_rules_path, classif, automatic_rescue_list,
                      rescue_df, strategy, thr=0.7):
    # Load data
    combined_class = merge_classifications(reference_rules_path,classif,strategy,thr)

    # Join rules
    mapping_hits = add_filter_results(mapping_hits, combined_class, classif)

    full_rescue_df = filter_mapping_hits(mapping_hits, strategy)
    
    full_inclusion_list = find_reintroduced_transcripts(full_rescue_df, automatic_rescue_list)

    rescue_df = merge_rescue_modes(rescue_df, full_rescue_df, full_inclusion_list, strategy)

    # write results
    # write_rescue_inclusion_list(full_inclusion_list, output_prefix)
    # write_rescue_full_table(rescue_df, output_prefix)

    return full_inclusion_list, rescue_df

    ### OLD PIPELINE (to generate the annoying table with useless information) ###
    # Filter mapping hits
    # reintroduced_final, rescued_mapping_final,\
    #       best_mapping_hits, rescued_ref, isoform_assoc_tr = filter_mapping_hits(mapping_hits, class_ref, classif,
    #                                                            automatic_rescued,strategy)

    # Process automatic rescue
    # automatic_fsm = process_automatic_rescue(classif, automatic_rescued, class_ref)

    # # Create rescue table
    # rescue_table = create_rescue_table(mapping_hits, best_mapping_hits, rescued_mapping_final,
    #                                    class_ref, rescued_ref, isoform_assoc_tr, automatic_fsm,
    #                                    strategy, threshold)
    # # Write rescue table
    # write_rescue_table(rescue_table, output_prefix)

