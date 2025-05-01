import pandas as pd
import numpy as np

from src.module_logging import rescue_logger

def rescue_fsm_monoexons(df):
    """Rescue monoexons that are classified as full-splice_match and are artifacts

    Args:
        df (DataFrame): DataFrame with SQANTI3 classification filtered to contain only artifacts

    Returns:
        DataFrame: DataFrame containing rescued monoexons reference transcripts
    """
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
    """
    Given a SQANTI3 classification file, returns the reference transcipts for which none of their 
    long-reads defined isoforms passed the SQANTI3 filter strategy applied
    Args:
        df (pd.DataFrame): SQANTI3 filter classification DataFrame 

    Returns:
        list: Reference transcripts that were lost during filtering
    """
    # Find all reference IDs in associated_transcript column
    all_ref = df['associated_transcript'].unique()
    # Check reference IDs not represented by isoforms (lost in filtering)
    isoform_ref = df[df['filter_result'] == 'Isoform']['associated_transcript'].unique().tolist() 
    rescue_logger.debug(f"Found {len(isoform_ref)} references represented")
    # Find lost references
    lost_ref = all_ref[~np.isin(all_ref, isoform_ref)]
    return lost_ref

def rescue_lost_reference(ref_id, classif):
    """
    Rescues a transcipt if there is any FSM associated with it.

    Args:
        ref_id (str): reference transcript candidate to automatic rescue
        classif (DataFrame): Artifact only classification DataFrame

    Returns:
        DataFrame: Single column dataframe with the rescued transcript ID
    """
    # Filter classification to get all the isoforms associated with a transcript
    classif_ref = classif[classif['associated_transcript'] == ref_id]
    
    # Check for FSM
    ref_check = any(classif_ref['structural_category'] == "full-splice_match")
    
    # If there is an FSM associated to the lost reference, return lost reference
    if ref_check:
        ref_df = pd.DataFrame({'isoform': [ref_id]})
        return ref_df
    

def save_automatic_rescue(rescue_df,class_df,mode,prefix):
    """Saves the automatic rescue results to files.

    Args:
        rescue_df (DataFrame): DataFrame containing rescued transcripts.
        class_df (DataFrame): DataFrame containing classification information.
        mode (str): Rescue mode of operation for saving results.
        prefix (str): Prefix for output file names.
    """
    # First write operation: rescue_auto without headers
    rescue_df.to_csv(
        f"{prefix}_automatic_rescued_list.tsv",
        sep='\t',
        header=False,
        index=False
    )

    if mode == 'automatic':
        
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

