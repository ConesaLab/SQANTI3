import pandas as pd
import numpy as np

def load_data(mapping_hits_path, reference_rules_path, sqanti_rules_classif_path, 
              automatic_rescue_path,strategy):
    """Load input data files."""
    mapping_hits = pd.read_csv(mapping_hits_path, sep="\t")
    
    if strategy == "rules":
        columns = ["isoform", "filter_result"]
    elif strategy == "ml":
        columns = ["isoform", "POS_MLprob"]
    class_ref = pd.read_csv(reference_rules_path, sep="\t", usecols=columns)
    classif = pd.read_csv(sqanti_rules_classif_path, sep="\t")
    # Add the reference classification to the transcript classification
    filter_combined = pd.concat([class_ref,classif[columns]])
    # Load automatic rescue results
    automatic_ref_rescued = pd.read_csv(automatic_rescue_path, sep="\t", header=0, names=["isoform", "associated_transcript","structural_category"])
    if automatic_ref_rescued.empty:
        # add a column with values "none" to avoid errors
        automatic_ref_rescued["associated_transcript"] = "none"
        
    return mapping_hits, class_ref, classif, filter_combined, automatic_ref_rescued

def add_filter_results(mapping_hits, filter_combined, classif):
    """Add the filter results to the mapping hits."""
    
    # Add filter result of mapping hits
    mapping_hits = mapping_hits.merge(
        filter_combined.rename(columns={"isoform": "mapping_hit"}), 
        on="mapping_hit", 
        how="left"
    ).rename(columns={"filter_result": "hit_filter_result",
                      "POS_MLprob": "hit_POS_MLprob"})

    # Add structural categories of candidates
    mapping_hits = mapping_hits.merge(
        classif[["isoform", "structural_category","associated_gene"]], 
        left_on="rescue_candidate", 
        right_on="isoform", 
        how="left"
    ).rename(columns={"structural_category": "candidate_structural_category"})
    mapping_hits.drop(columns=["isoform"], inplace=True)
    return mapping_hits

def filter_mapping_hits(mapping_hits, class_ref, classif, automatic_rescue_df,
                        strategy,threshold):
    """Filter mapping hits and remove redundant references."""
    # 1. Filter the mapping hits reference based on the results of filtering
    if strategy == "rules":
        # Apply rules-based filtering
        mapping_hits_filt = mapping_hits[mapping_hits["hit_filter_result"] == "Isoform"]
        mapping_hits_filt['score'] = mapping_hits_filt['alignment_score']
        
    elif strategy == "ml":
        mapping_hits_filt = mapping_hits[mapping_hits['hit_POS_MLprob'] >= threshold]
        mapping_hits_filt['score'] = mapping_hits_filt['alignment_score'] * mapping_hits_filt['hit_POS_MLprob']

    # Select the rows with the best score for each rescue candidate
    best_score = mapping_hits_filt.groupby('rescue_candidate')['score'].transform('max')
    mapping_hits_filt = mapping_hits_filt[mapping_hits_filt['score'] == best_score]

    mapping_hits_filt.to_csv("debug_mapping_hits_filt.tsv", sep="\t", index=False)
    input()

    # 2. Select only reference rescued transcripts
    rescued_ref = mapping_hits_filt[mapping_hits_filt["mapping_hit"].isin(class_ref["isoform"])]
 
    # 3. Remove reference transcripts already in the transcriptome

    # Retrieve all reference transcripts already represented by isoforms
    isoform_assoc_tr = classif[
        (classif["filter_result"] == "Isoform") & 
        (classif["associated_transcript"] != "novel")
    ][["associated_transcript"]]
    
    # Include those retrieved in automatic rescue
    isoform_assoc_tr = pd.concat([isoform_assoc_tr, automatic_rescue_df["associated_transcript"]]).drop_duplicates()
    # Find truly rescued references
    rescued_mapping_final = rescued_ref[
        ~rescued_ref["mapping_hit"].isin(isoform_assoc_tr["associated_transcript"])
    ][["mapping_hit"]].rename(columns={"mapping_hit": "ref_transcript"}).drop_duplicates()
    
    # Generate final list of rescued transcripts
    automatic_ref_rescued = automatic_rescue_df.rename(columns={"associated_transcript": "ref_transcript"})["ref_transcript"]
    rescued_final = pd.concat([automatic_ref_rescued, rescued_mapping_final]).drop_duplicates()
    return rescued_final, rescued_mapping_final, mapping_hits_filt, rescued_ref, isoform_assoc_tr

def write_rescue_inclusion_list(rescued_final,output_prefix):
    """Write the rescue inclusion list to a file."""
    output_path = f"{output_prefix}_full_inclusion_list.tsv"
    rescued_final.to_csv(output_path, sep="\t", index=False, header=False)

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

def rescue_by_mapping(mapping_hits_path, reference_rules_path, sqanti_rules_classif_path,
                      automatic_rescue_path, output_prefix,strategy, threshold=0.7):
    # Load data
    mapping_hits, class_ref, classif, \
        filter_combined, automatic_rescued = load_data(mapping_hits_path, reference_rules_path,
                                                       sqanti_rules_classif_path, automatic_rescue_path,strategy)

    # Join rules
    mapping_hits = add_filter_results(mapping_hits, filter_combined, classif)

    # Filter mapping hits
    rescued_final, rescued_mapping_final,\
          best_mapping_hits, rescued_ref, isoform_assoc_tr = filter_mapping_hits(mapping_hits, class_ref, classif,
                                                               automatic_rescued,strategy,threshold)

    # Write rescue inclusion list
    write_rescue_inclusion_list(rescued_final, output_prefix)

    # Process automatic rescue
    automatic_fsm = process_automatic_rescue(classif, automatic_rescued, class_ref)

    # Create rescue table
    rescue_table = create_rescue_table(mapping_hits, best_mapping_hits, rescued_mapping_final,
                                       class_ref, rescued_ref, isoform_assoc_tr, automatic_fsm,
                                       strategy, threshold)
    # Write rescue table
    write_rescue_table(rescue_table, output_prefix)

