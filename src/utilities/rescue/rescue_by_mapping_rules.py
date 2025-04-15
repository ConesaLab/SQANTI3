import pandas as pd

def load_data(mapping_hits_path, reference_rules_path, sqanti_rules_classif_path):
    """Load input data files."""
    mapping_hits = pd.read_csv(mapping_hits_path, sep="\t", names=["rescue_candidate", "mapping_hit", "sam_flag"])
    rules_ref = pd.read_csv(reference_rules_path, sep="\t", usecols=["isoform", "filter_result"])
    classif = pd.read_csv(sqanti_rules_classif_path, sep="\t")
    return mapping_hits, rules_ref, classif

def join_rules(mapping_hits, rules_ref, classif):
    """Join reference and long-read rules results."""
    rules_lr = classif[["isoform", "filter_result"]]
    rules = pd.concat([rules_ref, rules_lr], ignore_index=True)

    # Add filter result of mapping hits
    mapping_hits = mapping_hits.merge(
        rules.rename(columns={"isoform": "mapping_hit"}), 
        on="mapping_hit", 
        how="left"
    ).rename(columns={"filter_result": "hit_filter_result"})

    # Add structural categories of candidates
    mapping_hits = mapping_hits.merge(
        classif[["isoform", "structural_category"]], 
        left_on="rescue_candidate", 
        right_on="isoform", 
        how="left"
    ).rename(columns={"structural_category": "candidate_structural_category"})
    return mapping_hits

def filter_mapping_hits(mapping_hits, rules_ref, classif, automatic_rescue_path):
    """Filter mapping hits and remove redundant references."""
    # Filter mapping hits that passed rules
    mapping_hits_iso = mapping_hits[mapping_hits["hit_filter_result"] == "Isoform"]

    # Select only reference rescued transcripts
    rescued_ref = mapping_hits_iso[mapping_hits_iso["mapping_hit"].isin(rules_ref["isoform"])]

    # Retrieve all reference transcripts already represented by isoforms
    isoform_assoc_tr = classif[
        (classif["filter_result"] == "Isoform") & 
        (classif["associated_transcript"] != "novel")
    ][["associated_transcript"]]

    # Include those retrieved in automatic rescue
    automatic_ref_rescued = pd.read_csv(automatic_rescue_path, sep="\t", names=["associated_transcript"])
    isoform_assoc_tr = pd.concat([isoform_assoc_tr, automatic_ref_rescued]).drop_duplicates()

    # Find truly rescued references
    rescued_mapping_final = rescued_ref[
        ~rescued_ref["mapping_hit"].isin(isoform_assoc_tr["associated_transcript"])
    ][["mapping_hit"]].rename(columns={"mapping_hit": "ref_transcript"}).drop_duplicates()

    # Generate final list of rescued transcripts
    automatic_ref_rescued = automatic_ref_rescued.rename(columns={"associated_transcript": "ref_transcript"})
    rescued_final = pd.concat([automatic_ref_rescued, rescued_mapping_final]).drop_duplicates()
    return rescued_final, mapping_hits_iso, rescued_mapping_final

def write_rescue_inclusion_list(rescued_final,output_prefix):
    """Write the rescue inclusion list to a file."""
    output_path = f"{output_prefix}_rescue_inclusion-list.tsv"
    rescued_final.to_csv(output_path, sep="\t", index=False, header=False)

def process_automatic_rescue(classif, automatic_ref_rescued, rules_ref):
    """Process automatic rescue results."""
    automatic_fsm = classif.merge(
        automatic_ref_rescued.rename(columns={"ref_transcript": "associated_transcript"}), 
        on="associated_transcript", 
        how="right"
    )
    automatic_fsm = automatic_fsm.merge(
        rules_ref.rename(columns={"isoform": "associated_transcript"}), 
        on="associated_transcript", 
        how="left"
    )
    automatic_fsm["sam_flag"] = None
    automatic_fsm["rescue_result"] = "rescued_automatic"
    automatic_fsm["exclusion_reason"] = None
    return automatic_fsm

def create_rescue_table(mapping_hits, rescued_mapping_final, isoform_assoc_tr, automatic_fsm):
    """Create the rescue table."""
    rescue_table = mapping_hits.copy()
    rescue_table["rescue_result"] = rescue_table["mapping_hit"].apply(
        lambda x: "rescued_mapping" if x in rescued_mapping_final["ref_transcript"].values else "not_rescued"
    )
    rescue_table["exclusion_reason"] = rescue_table.apply(
        lambda row: "artifact_by_rules" if row["mapping_hit"] not in mapping_hits["mapping_hit"].values else
                    "long_read_transcript" if row["mapping_hit"] not in rescued_mapping_final["ref_transcript"].values else
                    "reference_already_present" if row["mapping_hit"] in isoform_assoc_tr["associated_transcript"].values else None,
        axis=1
    )
    rescue_table = pd.concat([rescue_table, automatic_fsm], ignore_index=True)
    return rescue_table

def write_rescue_table(rescue_table, output_prefix):
    """Write the rescue table to a file."""
    output_path = f"{output_prefix}_rescue_table.tsv"
    rescue_table.to_csv(output_path, sep="\t", index=False)

def rescue_rules(mapping_hits_path, reference_rules_path, sqanti_rules_classif_path, automatic_rescue_path, output_prefix):
    # Load data
    mapping_hits, rules_ref, classif = load_data(mapping_hits_path, reference_rules_path, sqanti_rules_classif_path)

    # Join rules
    mapping_hits = join_rules(mapping_hits, rules_ref, classif)

    # Filter mapping hits
    rescued_final, mapping_hits_iso, rescued_mapping_final = filter_mapping_hits(
        mapping_hits, rules_ref, classif, automatic_rescue_path
    )

    # Write rescue inclusion list
    write_rescue_inclusion_list(rescued_final, output_prefix)

    # Process automatic rescue
    automatic_fsm = process_automatic_rescue(classif, rescued_final, rules_ref)

    # Create rescue table
    rescue_table = create_rescue_table(mapping_hits, rescued_mapping_final, classif, automatic_fsm)

    # Write rescue table
    write_rescue_table(rescue_table, output_prefix)

