import pandas as pd

from src.module_logging import rescue_logger

def get_lost_reference_id(df):
    """
    Get the reference IDs that have lost all their isoforms during filtering.
    
    :param df: SQANTI classification dataframe in Pandas
    """
    all_ref = df['associated_transcript'].unique()
    # Find all references were one of their transcripts is classified as "Isoforms"  
    isoform_ref = df.loc[df['filter_result'].eq('Isoform'), 'associated_transcript'].unique()
    rescue_logger.debug(f"Found {isoform_ref.size} references represented")
    # Hash-based difference, sorted like np.setdiff1d; on string IDs numpy
    # compares every pair of IDs
    return pd.Index(all_ref).difference(isoform_ref).to_numpy()

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
