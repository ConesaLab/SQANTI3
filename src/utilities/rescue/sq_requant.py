import pandas as pd
import warnings
warnings.filterwarnings("ignore") # TODO: See why pandas raises the warning on copying

from collections import defaultdict

from src.module_logging import rescue_logger
from src.utilities.rescue.requant_helpers import (
    build_artifact_table, export_counts, redistribute_counts, calculate_tpm
)
def parse_files(count_file):
    # Load counts file
    counts = pd.read_csv(count_file, sep = '\t', comment = '#')
    counts.columns = ['transcript_id', 'count']
    return(counts)


def run_requant(counts, rescue_df, classif_df, prefix):
    # Creation of a table with all artifacts (rescued and not rescued)
    artifacts_df = build_artifact_table(rescue_df, classif_df)
    old_counts = counts.set_index('transcript_id')['count'].to_dict()
    # Redistribute counts based on the rescue results
    new_counts = redistribute_counts(artifacts_df, classif_df, old_counts)
    
    return export_counts(old_counts, new_counts, prefix)

def to_tpm(counts_df, class_df, prefix):

    class_df.rename(columns={'isoform': 'transcript_id'}, inplace=True)
    class_df = class_df[["transcript_id","length"]]
    
    counts_d = defaultdict(int)
    counts_df.apply(lambda x: counts_d.update({x.iloc[0] : x.iloc[1]}), axis = 1) 
    class_df['counts'] = class_df['transcript_id'].apply(lambda x: counts_d[x]) # Apparently this line triggers a SettingWithCopyWarning
    #remove zero values
    class_df = class_df[class_df['counts'] != 0]
    # Calculate TPM
    class_df['TPM'] = calculate_tpm(class_df['counts'], class_df['length'])
    final = class_df[['transcript_id', 'TPM']]
    final.to_csv(f"{prefix}_reassigned_counts_TPM.tsv", header = True, sep='\t', index=False)