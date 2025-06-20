import pandas as pd
import numpy as np

from src.module_logging import rescue_logger

def read_classification(filename):
    return pd.read_csv(filename, sep="\t")

def rescue_fsm_monoexons(df):
    df_mono = df[
        (df['structural_category'] == 'full-splice_match') &
        (df['exons'] == 1) &
        (df['filter_result'] == 'Artifact')
    ]
    # Convert to rescued
    rescue_mono = df_mono[['associated_transcript']].drop_duplicates()
    rescue_mono = rescue_mono.rename(columns={'associated_transcript': 'isoform'})
    return rescue_mono

def add_ism_monoexons(df_filt,df):
    classif_ism_mono = df[
        (df['structural_category'] == 'incomplete-splice_match') &
        (df['exons'] == 1)
    ]
    return pd.concat([df_filt, classif_ism_mono])

def rescue_monoexons(df,me):
    df_mono = df[
        (df['exons'] == 1) &
        (df['filter_result'] == 'Artifact')
    ]
    if me == 'fsm':
        df_mono = df_mono[df_mono['structural_category'] == 'full-splice_match']
    return df_mono

def get_lost_reference_id(df):
    # Find all reference IDs in associated_transcript column
    all_ref = df['associated_transcript'].unique()
    # Check reference IDs not represented by isoforms (lost in filtering)
    isoform_ref = df[df['filter_result'] == 'Isoform']['associated_transcript'].unique().tolist() 
    rescue_logger.debug(f"Found {len(isoform_ref)} references represented")
    # Find lost references
    lost_ref = all_ref[~np.isin(all_ref, isoform_ref)]
    return lost_ref

def rescue_lost_reference(ref_id, classif):
    # Filter classification
    classif_ref = classif[classif['associated_transcript'] == ref_id]
    
    # Check for FSM
    ref_check = any(classif_ref['structural_category'] == "full-splice_match")
    
    # If there is an FSM associated to the lost reference, return lost reference
    if ref_check:
        ref_df = pd.DataFrame({'isoform': [ref_id]})
        return ref_df
    
    # If there is no FSM associated, return ISM artifacts
    else:
        ism_df = classif_ref[['isoform']]
        return ism_df
    
def save_automatic_rescue(rescue_df,class_df,mode,prefix):
    # First write operation: rescue_auto without headers
    rescue_df.to_csv(
        f"{prefix}_automatic_inclusion_list.tsv",
        sep='\t',
        header=False,
        index=False
    )

    if mode == 'automatic':
        if rescue_df.iloc[0,0] == "none":
            rescue_logger.info("No FSM mono-exonic artifacts found for automatic rescue.")
        else:
            # Second write operation: rescue_table with headers 
            rescue_table = class_df[
                class_df['associated_transcript'].isin(rescue_df['isoform'])
            ][['isoform', 'associated_transcript', 'structural_category']].rename(
                columns={
                    'isoform': 'artifact',
                    'associated_transcript': 'rescued_transcript'
                }
            )

            rescue_table.to_csv(
                f"{prefix}_automatic_rescue_table.tsv",
                sep='\t',
                index=False
            )

