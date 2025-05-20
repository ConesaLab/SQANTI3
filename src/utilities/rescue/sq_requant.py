import pandas as pd
import os
from collections import defaultdict
import sys
import warnings
warnings.filterwarnings("ignore")

def parse_files(gtf_path,inclusion_file,count_file,prefix,mode):

    col_names = [
        "Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"
    ]
    rescue_gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        names=col_names,
        dtype={"Chromosome": str, "Start": int, "End": int},
        low_memory=False
    )
    rescue_gtf["transcript_id"] = rescue_gtf["Attribute"].str.extract(r'transcript_id "([^"]+)"')
    
    # Inclusion list
    inclusion_df = pd.read_csv(inclusion_file, names = ['transcript_id'])

    # Load counts file
    counts = pd.read_csv(count_file, sep = '\t', comment = '#')
    counts.columns = ['transcript_id', 'count']

    #Rescued table
    rescued_path=f"{prefix}_{mode}_rescue_table.tsv"
    rescued_table = pd.read_csv(rescued_path, sep = '\t')


    return(rescue_gtf, inclusion_df, counts, rescued_table)

def select_hit(isoform, rescued, new_counts, old_counts):
    group = rescued[rescued["mapping_hit"] == isoform]
    group['counts'] = group.apply(lambda x: old_counts.get(x['rescue_candidate'], 0), axis = 1)
    new_counts[isoform] = group['counts'].sum()
    return new_counts

def reassign_counts(isoform, old_counts, new_counts, list_of_changed):
    if isoform in new_counts.keys():
        count = new_counts[isoform]
        list_of_changed.append(isoform)
    else:
        try:
            count = old_counts[isoform]
        except KeyError:
            count = 0
    return(count)

def fill_old_counts(isoform, old_counts):
    try:
        count = old_counts[isoform]
    except KeyError:
        count = 0
    return(count)
    
def run_requant(rescue_gtf, inclusion_df, counts, 
                rescue_df, prefix):


    #select list of isoforms that were rescued by SQANTI3 rescue
    rescued_iso = rescue_df[rescue_df['rescue_result'] != 'not_rescued']
    additional_isoforms = rescue_df[rescue_df['exclusion_reason'] == 'reference_already_present']
    rescued = pd.concat([rescued_iso, additional_isoforms])

    # Selec only isoforms that were not rescued or had a transcript that was already rescued
    summary_df = (
    rescue_df.groupby('rescue_candidate', as_index=False
                     ).apply(lambda g: pd.Series({
        'rescued': ((g['rescue_result'] != 'not_rescued') | 
                    (g['exclusion_reason'].fillna('')  == 'reference_already_present')).any(),
        'associated_gene': g['associated_gene'].unique()[0]  # assumes unique
    }))
    .reset_index(drop=True)
    )
    not_rescued = summary_df[summary_df['rescued'] == False]
    not_rescued['mapping_hit'] = not_rescued['associated_gene'] + '_TD'
    not_rescued = not_rescued.drop_duplicates(['rescue_candidate', 'mapping_hit'])
    

    # Mixing both we can include the counts for the transcript divergence of each gene
    rescued = pd.concat([rescued, not_rescued])
    rescued = rescued.drop_duplicates(['rescue_candidate', 'mapping_hit'])

    not_rescued = not_rescued[['rescue_candidate', 'mapping_hit']]
    not_rescued.columns = ['rescue_candidate', 'transcript_id']

    #create dictionaries of old and new counts
    old_counts = counts.set_index('transcript_id')['count'].to_dict()
    new_counts = defaultdict(int)
    
    #reassign counts to surviving isoforms
    rescued = rescued[rescued['rescue_candidate'].isin(old_counts.keys())]
    inclusion_df = pd.concat([inclusion_df, not_rescued['transcript_id']])
    inclusion_df.apply(lambda x: select_hit(x.iloc[0], rescued, new_counts, old_counts), axis = 1)
    #combine old and new counts
    counts_df = pd.DataFrame(rescue_gtf['transcript_id'].unique(), columns = ['transcript_id'])
    counts_df = pd.concat([counts_df,not_rescued['transcript_id']])
    counts_df['old_count'] = counts_df['transcript_id'].apply(lambda x: fill_old_counts(x, old_counts))
    list_of_changed = []
    counts_df['new_count'] = counts_df['transcript_id'].apply(lambda x: reassign_counts(x, old_counts, new_counts, list_of_changed))
    changed = pd.DataFrame()
    changed['changed_count'] = list_of_changed

    counts_df.to_csv(f"{prefix}_reassigned_counts_extended.tsv", header = True, index = False, sep = '\t')
    changed.to_csv(f"{prefix}_changed_counts.tsv", header = True, index = False, sep = '\t')

    counts_df_short = counts_df[['transcript_id', 'new_count']]
    counts_df_short.to_csv(f"{prefix}_reassigned_counts.tsv", header = True, index = False, sep = '\t')

    # Save transcript divergency relations
    not_rescued.to_csv(f"{prefix}_transcript_divergency_relations.tsv", header = True, index = False, sep = '\t')

def to_tpm(rescue_gtf, prefix):
    def calculate_tpm(counts, lengths):
        # Convert lengths to kilobases
        lengths_kb = lengths / 1000
        # Calculate RPK
        rpk = counts / lengths_kb
        # Calculate Total RPK
        total_rpk = rpk.sum()
        # Calculate TPM
        tpm = (rpk / total_rpk) * 1e6
        return tpm
    #check if requantification file was generated
    counts_df_path=f"{prefix}_reassigned_counts.tsv"
    if not os.path.isfile(counts_df_path):
        print(f"ERROR: {counts_df_path} doesn't exist. Abort!", file=sys.stderr)
    else:
        counts_df = pd.read_csv(counts_df_path, sep = '\t')

    #calculate transcripts length
    exons = rescue_gtf[rescue_gtf['Feature'] == 'exon']
    # Calculate exon lengths
    exons['exon_length'] = exons['End'] - exons['Start'] + 1
    # Sum exon lengths for each transcript to get the total transcript length
    transcript_lengths = exons.groupby('transcript_id')['exon_length'].sum().reset_index()
    transcript_lengths.columns = ['Transcript_ID', 'Length']
    lengths_d = transcript_lengths.set_index('Transcript_ID')['Length'].to_dict()
    
    rescue_gtf = rescue_gtf[rescue_gtf.Feature == 'transcript']
    rescue_gtf['t_length'] = rescue_gtf['transcript_id'].apply(lambda x: lengths_d[x])
    
    counts_d = defaultdict(int)
    counts_df.apply(lambda x: counts_d.update({x.iloc[0] : x.iloc[1]}), axis = 1)
    rescue_gtf['counts'] = rescue_gtf['transcript_id'].apply(lambda x: counts_d[x])
    #remove zero values
    rescue_gtf = rescue_gtf[rescue_gtf['counts'] != 0]
    # Calculate TPM
    rescue_gtf['TPM'] = calculate_tpm(rescue_gtf['counts'], rescue_gtf['t_length'])
    final = rescue_gtf[['transcript_id', 'TPM']]
    final.to_csv(f"{prefix}_reassigned_counts_TPM.tsv", header = True, sep='\t', index=False)