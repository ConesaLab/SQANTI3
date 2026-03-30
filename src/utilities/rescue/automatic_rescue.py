import pandas as pd
import numpy as np

from src.module_logging import rescue_logger

def get_lost_reference_id(df):
    """
    Get the reference IDs that have lost all their isoforms during filtering.
    
    :param df: SQANTI classification dataframe in Pandas
    """
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
    
def generate_automatic_table(inclusion_df,class_df):
    """
    Generate the rescue table with the results of automatic rescue

    :param inclusion_df: Pandas Series with the isoforms to be reintrouced
    :param class_df: SQANTI3 classification DataFrame (after filtering)
    """
    # Handle case with nothing to be rescued in the automatic mode
    if inclusion_df.iloc[0,0] == "none":
        rescue_logger.info("No FSM mono-exonic artifacts found for automatic rescue.")
        inclusion_df["isoform"] = "none"

    rescue_table = class_df[
        class_df['associated_transcript'].isin(inclusion_df['isoform']) & class_df['filter_result'].eq('Artifact')
    ][['isoform', 'associated_transcript']].rename(
        columns={'isoform': 'artifact', 'associated_transcript': 'assigned_transcript'}
    )
    # Adding the important columns
    rescue_table["rescue_mode"] = "automatic"
    rescue_table["origin"] = "reference"
    rescue_table["reintroduced"] = (
        ~rescue_table.duplicated("assigned_transcript")
    ).map({True: "yes", False: "no"})

    return rescue_table
