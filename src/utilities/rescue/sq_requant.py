import pandas as pd
import argparse
import os
from collections import defaultdict
import sys
import warnings
warnings.filterwarnings("ignore")

def parse_files(args):
    #Check if the rescued GTF file was generated
    gtf_path="{0}/{1}_rescued.gtf".format(args.dir, args.output)
    if not os.path.isfile(gtf_path):
        print("ERROR: {0} doesn't exist. Abort!".format(gtf_path), file=sys.stderr)
    else:
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

    #check if inclusion list was generated
    inclusion_list_path = args.dir + "/" + args.output + "_rescue_inclusion-list.tsv"
    if not os.path.isfile(inclusion_list_path):
        print("ERROR: {0} doesn't exist. Abort!".format(inclusion_list_path), file=sys.stderr)
        sys.exit(-1)
    else:
        inclusion_list = pd.read_csv(inclusion_list_path, sep = '\t')

    #cheeck if counts file exists
    if not os.path.isfile(args.counts):
        print("ERROR: {0} doesn't exist. Abort!".format(args.counts), file=sys.stderr)
    else:
        counts = pd.read_csv(args.counts, sep = '\t', comment = '#')
        counts.columns = ['transcript_id', 'count']

    #check if rescued table exists
    rescued_path="{0}/{1}_rescue_table.tsv".format(args.dir, args.output)
    if not os.path.isfile(rescued_path):
        print("ERROR: {0} doesn't exist. Abort!".format(rescued_path), file=sys.stderr)
    else:
        rescued = pd.read_csv(rescued_path, sep = '\t')
    return(rescue_gtf, inclusion_list, counts, rescued)

def run_requant(rescue_gtf, inclusion_list, counts, rescued, out, output):
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
    #select list of isoforms that were rescued by SQANTI3 rescue
    rescued_iso = rescued[rescued['rescue_result'] != 'not_rescued']
    additional_isoforms = rescued[rescued['exclusion_reason'] == 'reference_already_present']
    rescued = pd.concat([rescued_iso, additional_isoforms])
    #rescued = rescued.drop_duplicates(['rescue_candidate', 'mapping_hit'])
    rescued = rescued.drop_duplicates(['rescue_candidate', 'mapping_hit'])      
    #create dictionaries of old and new counts
    old_counts = counts.set_index('transcript_id')['count'].to_dict()
    new_counts = defaultdict(int)
    #reassign counts to surviving isoforms
    rescued = rescued[rescued['rescue_candidate'].isin(old_counts.keys())]
    inclusion_list.apply(lambda x: select_hit(x.iloc[0], rescued, new_counts, old_counts), axis = 1)
    #combine old and new counts
    counts_df = pd.DataFrame(rescue_gtf['transcript_id'].unique(), columns = ['transcript_id'])
    counts_df['old_count'] = counts_df['transcript_id'].apply(lambda x: fill_old_counts(x, old_counts))
    list_of_changed = []
    counts_df['new_count'] = counts_df['transcript_id'].apply(lambda x: reassign_counts(x, old_counts, new_counts, list_of_changed))
    changed = pd.DataFrame()
    changed['changed_count'] = list_of_changed

    counts_df.to_csv("{}/{}_reassigned_counts_extended.tsv".format(out, output), header = True, index = False, sep = '\t')
    changed.to_csv("{}/{}_changed_counts.tsv".format(out, output), header = True, index = False, sep = '\t')

    counts_df_short = counts_df[['transcript_id', 'new_count']]
    counts_df_short.to_csv("{}/{}_reassigned_counts.tsv".format(out, output), header = True, index = False, sep = '\t')

def to_tpm(rescue_gtf, out, output):
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
    counts_df_path="{}/{}_reassigned_counts.tsv".format(out, output)
    if not os.path.isfile(counts_df_path):
        print("ERROR: {0} doesn't exist. Abort!".format(counts_df_path), file=sys.stderr)
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
    final.to_csv("{}/{}_reassigned_counts_TPM.tsv".format(out, output), header = True, sep='\t', index=False)