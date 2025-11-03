import pandas as pd
import numpy as np

from src.module_logging import rescue_logger

def rescue_fsm_monoexons(df):
    df_mono = df[
        (df['structural_category'] == 'full-splice_match') &
        (df['exons'] == 1)
    ]
    rescue_mono = pd.DataFrame()
    lost_mono = get_lost_reference_id(df_mono)
    # Rescue mono
    for ref_id in lost_mono:
        rescue_mono = pd.concat([rescue_mono,rescue_lost_reference(ref_id, df_mono)])
    return rescue_mono

def get_lost_reference_id(df):
    all_ref = df['associated_transcript'].to_numpy()
    # Find all references were one of their transcripts is classified as "Isoforms"  
    isoform_ref = df.loc[df['filter_result'].eq('Isoform'), 'associated_transcript'].unique()
    rescue_logger.debug(f"Found {isoform_ref.size} references represented")
    return np.setdiff1d(np.unique(all_ref), isoform_ref)

def rescue_lost_reference(ref_id, classif):
    # Filter classification
    classif_ref = classif[classif['associated_transcript'] == ref_id]
    
    # Check for FSM
    ref_check = any(classif_ref['structural_category'] == "full-splice_match")
    
    # If there is an FSM associated to the lost reference, return lost reference
    if ref_check:
        ref_df = pd.DataFrame({'isoform': [ref_id]})
        return ref_df
    
def save_automatic_rescue(rescue_df,class_df,prefix):
    # First write operation: rescue_auto without headers
    rescue_df.to_csv(
        f"{prefix}_automatic_inclusion_list.tsv",
        sep='\t',
        header=False,
        index=False
    )
    if rescue_df.iloc[0,0] == "none":
        rescue_logger.info("No FSM mono-exonic artifacts found for automatic rescue.")
        rescue_df["isoform"] = "none"
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
    if rescue_df.iloc[0,0] == "none":
        rescue_logger.info("No FSM mono-exonic artifacts found for automatic rescue.")

