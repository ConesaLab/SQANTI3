
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:32:18 2024

@author: nkeil
"""


import pandas as pd
import argparse
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.ticker import FixedLocator

## Update options
def getOptions():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Make sqanti reads QC summary tables and plots')
    
    parser.add_argument('-r','--ref', type=str, dest="inREF" ,required=True, help='Path to the reference gtf. This should be the same reference used for sqanti3 QC')
    parser.add_argument('-d', '--design', type=str, dest="inDESIGN" ,required=True, help='Path to design file, must have sampleID column, classification_file column and junction_file column')
    parser.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='Optional: This is the column name that plots are to be faceted by')
    parser.add_argument('-o','--out', type=str,  dest="OUT", required=True, help='Output directory for saving plots and tables')
    parser.add_argument('-p','--prefix', type=str, dest="PREFIX", required=True, help='Output filename prefix')
    parser.add_argument('-gp', '--generate-plots', action='store_true', dest="GENERATE_PLOTS", help='Generate plots if set (default: True)', default=True)
    parser.add_argument('-ge','--gene-expression', type=int, dest="ANNOTEXP", required=False, help='Expression cut off level for determining underannotated genes', default = 10)
    parser.add_argument('-je','--jxn-expression', type=int, dest="JXNEXP", required=False, help='Expression cut off level for junction cv plots', default = 3)
    #parser.add_argument('-ps','--perc-samples', type=int, dest="PERCSAMPLE", required=False, help='Percent samples for determining underannotated genes', default = 80)
    parser.add_argument('-pn','--perc-novel', type=int, dest="PERCNOVEL", required=False, help='Percent NIC+NNC for determining underannotated genes', default = 90)
    
    args = parser.parse_args()
    return args

def load_sqanti_file(file,col_Lst, dtype_Dct):
    return pd.read_csv(file, sep="\t", usecols=col_Lst, dtype=dtype_Dct, low_memory=True)


def merge_dfs(df1, df2, column1, column2, how='outer'):
    """
    Merge two pandas DataFrames on specified columns with different IDs.
    Returns:
    - Merged DataFrame.
    """
    # Rename the column in df2 to match df1 for the merge
    df2_renamed = df2.rename(columns={column2: column1})
    
    # Perform the merge
    merged_df = pd.merge(df1, df2_renamed, on=column1, how=how)
    
    return merged_df

def merge_three_dfs(df1, df2, df3, col1, col2, col3):
    """
    Merge three DataFrames using an outer join based on one specified column from each DataFrame,
    ensuring these columns contain the same unique set of values.
    
    Parameters:
    - df1, df2, df3: the DataFrames to merge.
    - col1: the column name in df1 to merge on.
    - col2: the column name in df2 to merge on.
    - col3: the column name in df3 to merge on.
    
    Returns:
    - Merged DataFrame with all rows and columns from df1, df2, and df3.
    """
    
    # Check if the merge columns are unique in each DataFrame
    if not df1[col1].is_unique or not df2[col2].is_unique or not df3[col3].is_unique:
        raise ValueError("Merge columns must be unique in each DataFrame.")
    
    # Merge df1 and df2 on their specified columns, then rename the column in df2 to match df1 for simplicity
    merged_df1_df2 = pd.merge(df1, df2.rename(columns={col2: col1}), how='outer', on=col1)
    
    # Merge the result with df3, renaming the column in df3 to match df1
    final_merged_df = pd.merge(merged_df1_df2, df3.rename(columns={col3: col1}), how='outer', on=col1)
    
    return final_merged_df

def categorize_by_readcount(read_count):
    if read_count == 1:
        return '1 read'
    elif read_count <= 10:
        return '2-10 reads'
    elif read_count <= 50:
        return '11-50 reads'
    elif read_count <= 100:
        return '50-100 reads'
    else:
        return '100+ reads'
    
def cv(x):
    if np.mean(x) == 0:
        return np.nan
    else:
        return np.std(x) / np.mean(x)
    
def calc_jxn_cv(jxnDF, classDF, refDF, dropFlag=True):
  
    """
    Calculate cv of reference donors and acceptors from sqanti jxn file
    
    Parameters:
    - jxnDF, dataframe of sqanti junction file
    - classDF: dataframe of sqanti classification file
    - refDF: dataframe with column 'gene_id' that contains all reference genes
    - dropFlag: flag to drop jxns in novel genes or jxns associated with multiple genes
    
    Returns:
    - Dataframe with cv of each reference donor and acceptor 
    """
    ##Drop all junctions where diff to nearest ref junction cannot be calculated
    filter_jxnDF=jxnDF.dropna(subset=['diff_to_Ref_start_site','diff_to_Ref_end_site'])
    
    ##Merge in gene id for all reads
    filter_jxnDF=pd.merge(filter_jxnDF, classDF[['isoform', 'associated_gene']], on='isoform', how='left')
    
    ##Calculate reference junction start and end
    filter_jxnDF['ref_junction_start']=(filter_jxnDF['genomic_start_coord']+filter_jxnDF['diff_to_Ref_start_site']).astype(int)
    filter_jxnDF['ref_junction_end']=(filter_jxnDF['genomic_end_coord']+filter_jxnDF['diff_to_Ref_end_site']).astype(int)
    filter_jxnDF['abs_diff_to_start']=abs(filter_jxnDF['diff_to_Ref_start_site'])
    filter_jxnDF['abs_diff_to_end']=abs(filter_jxnDF['diff_to_Ref_end_site'])
    
    ##Calculate mean, std and cv - Starts
    cv_startDF =filter_jxnDF.groupby(['chrom', 'strand', 'ref_junction_start']).agg({'abs_diff_to_start': ['mean', 'std', cv, 'size'], 'associated_gene': [lambda x: '|'.join(set(x)), lambda x: len(set(x)) ]}).reset_index()
    cv_startDF.columns = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff','cv', 'count','associated_gene','gene_count']
    
    ##Create flags
    cv_startDF['flag_multi_gene'] = cv_startDF['gene_count'].apply(lambda x: 1 if x > 1 else 0)
    cv_startDF['flag_annotated_gene'] = np.where(cv_startDF['associated_gene'].isin(refDF['gene_id'].unique()), 1, 0) 
    
    cv_startDF['flag_donor'] = cv_startDF['strand'].apply(lambda x: 1 if x == '+' else 0)
    cv_startDF['flag_acceptor'] = cv_startDF['strand'].apply(lambda x: 1 if x == '-' else 0)
    cv_startDF['flag_single'] = cv_startDF['count'].apply(lambda x: 1 if x == 1 else 0)
    cv_startDF['flag_mean_0'] = cv_startDF['mean_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    cv_startDF['flag_std_0'] = cv_startDF['std_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    cv_startDF['flag_cv_0'] = cv_startDF['cv'].apply(lambda x: 1 if x == 0 else 0)
    
    ##Calculate mean, std and cv - Ends
    cv_endDF =filter_jxnDF.groupby(['chrom', 'strand', 'ref_junction_end']).agg({'abs_diff_to_end': ['mean', 'std', cv, 'size'], 'associated_gene': [lambda x: '|'.join(set(x)), lambda x: len(set(x)) ]}).reset_index()
    cv_endDF.columns = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff','cv', 'count','associated_gene','gene_count']
    
    ##Create flags
    cv_endDF['flag_multi_gene'] = cv_endDF['gene_count'].apply(lambda x: 1 if x > 1 else 0)
    cv_endDF['flag_annotated_gene'] =  np.where(cv_endDF['associated_gene'].isin(refDF['gene_id'].unique()), 1, 0)
    
    cv_endDF['flag_donor'] = cv_endDF['strand'].apply(lambda x: 1 if x == '-' else 0)
    cv_endDF['flag_acceptor'] = cv_endDF['strand'].apply(lambda x: 1 if x == '+' else 0)
    cv_endDF['flag_single'] = cv_endDF['count'].apply(lambda x: 1 if x == 1 else 0)
    cv_endDF['flag_mean_0'] = cv_endDF['mean_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    cv_endDF['flag_std_0'] = cv_endDF['std_abs_diff'].apply(lambda x: 1 if x == 0 else 0)
    cv_endDF['flag_cv_0'] = cv_endDF['cv'].apply(lambda x: 1 if x == 0 else 0)
    
    cvDF=pd.concat([cv_startDF,cv_endDF])

    if dropFlag:
         
        #count=((cvDF['flag_annotated_gene'] != 1) | (cvDF['flag_multi_gene'] == 1)).sum()
        #print(str(count)+"donors and acceptors assigned to multiple and/or novel genes")
        
        drop=(cvDF['flag_annotated_gene'] != 1) | (cvDF['flag_multi_gene'] == 1)
        cvDF_filtered=cvDF[~drop]
        
        return cvDF_filtered
    else :
        return cvDF
    
def flag_ref_monoexon(inRef):
    
    """
    Creates a dataframe with all reference genes and identifies genes with monoexon transcri[td]
    
    Parameters:
    - inRef, path to reference gtf
    
    Returns:
    - Dataframe with all reference gene ids and flag if it has a monoexon transcript or not
    """
    # Get exon counts from reference GTF
    
    refgtf = pd.read_csv(inRef,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    refexon = refgtf[refgtf["feature"]=="exon"].copy()
#    print("Total lines in reference GTF= "+str(len(gtf)))

    # Get gene_id and transcript_id values from attributes column
    for i in refexon.index:
        raw_attrs = refexon.at[i, 'attribute']
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = np.nan, np.nan
        for item in g_t_attrs:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
            elif 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
        if transcript_id == np.nan and refexon.at[i, 'feature'] != "gene" :
            print("WARNING: transcript_id not found in {}".format(refexon[i]))
        refexon.at[i, 'gene_id'] = str(gene_id)
        refexon.at[i, 'transcript_id'] = str(transcript_id)
#    print("Total transcripts in original GTF= "+str(refexon["transcript_id"].nunique()))
#    print("Total genes in original GTF= "+str(refexon["gene_id"].nunique()))

    # Get min exons per transcript for each gene to flag genes with at least one monoexon transcript
    refxcrpt = refexon.groupby(["gene_id", "transcript_id"])["feature"].count().reset_index().rename(columns={"feature": "num_exon"})
    refgene = refxcrpt.groupby("gene_id")["num_exon"].min().reset_index().rename(columns={"num_exon": "min_exon"})
    refgene["flag_ref_monoexon"] = np.where(
        refgene["min_exon"] == 1,
        1,
        0
    )
    return refgene

def proc_samples(design_file, ref):
    # Read design file
    design_DF = pd.read_csv(design_file, sep=",")
    

    
    # Make dictionaries to store each file type
    gene_count_dfs = {}
    ujc_count_dfs = {}
    length_dfs = {}
    err_dfs = {}
    cv_dfs = {}
    fsm_dfs = {}
    ism_dfs = {}
    nic_nnc_dfs = {}
    nov_can_dfs = {}
    length_Dct ={}
    
    
    if args.inFACTOR == None:
        exp_factor = 'temp_factor'
    else:
        exp_factor = args.inFACTOR
    
    
    ##Flag ref genes with at least one mono exonic transcript
    ref_DF = flag_ref_monoexon(ref)
    
    ## CREATE SUMMARY FILES TO MAKE PLOTS FROM
    jxn_cols = ['isoform','chrom','strand','junction_number','genomic_start_coord','genomic_end_coord','junction_category',
                'diff_to_Ref_start_site','diff_to_Ref_end_site','canonical','junction_category']
    
    jxn_dtypes = {'isoform':'string', 'chrom': 'string', 'strand': 'string', 'junction_number': 'string',
                  'genomic_start_coord':'Int64', 'genomic_end_coord': 'Int64', 'junction_category':'string',
                  'diff_to_Ref_start_site': 'Int64', 'diff_to_Ref_end_site': 'Int64', 'canonical': 'string'}
    
    class_cols = ['isoform','chrom','strand','exons','associated_gene','associated_transcript','structural_category','subcategory',
                    'length', 'RTS_stage','perc_A_downstream_TTS','ref_length','ref_exons','all_canonical']
    
    ##!!! predicted NMD dtype? - Not su
    class_dtypes = {'isoform': 'string', 'chrom': 'string', 'strand': 'string', 'exons': 'Int64', 'associated_gene': 'string','associated_transcript': 'string', 
                    'structural_category': 'string', 'subcategory': 'string','length': 'Int64', 'RTS_stage': 'boolean', 'perc_A_downstream_TTS': float, 
                    'ref_length': 'Int64','ref_exons': 'Int64', 'all_canonical': 'string'}
    
    for index, row in design_DF.iterrows():
        # gtf = row['gtf_file']
        sampleID = row['sampleID']
        class_file = row['classification_file']
        jxn_file = row['junction_file']
        
        if exp_factor == 'temp_factor' :
            exp_factor_val = 0
        else:
            exp_factor_val = row[exp_factor]
        
        print("Loading junction file: "+ sampleID)
        jxn_DF = load_sqanti_file(jxn_file, jxn_cols, jxn_dtypes)
        #jxn_DF = load_sqanti_file(jxn_file, jxn_cols)
        
        print("Loading classification file: "+ sampleID)
        class_DF = load_sqanti_file(class_file, class_cols, class_dtypes)
        #class_DF = load_sqanti_file(class_file, class_cols)
    
        ##Merge in ref DF
        class_DF = merge_dfs(class_DF, ref_DF, 'associated_gene', 'gene_id', 'left')
    
        ##Get number of novel, known canonical and non-canonical junctions
        count_jxns_DF = jxn_DF.pivot_table(index='isoform',
                                           columns=['junction_category', 'canonical'],
                                           aggfunc='size',
                                           fill_value=0).reset_index()
        count_jxns_DF.columns = ['isoform', 'known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
    
    
        #Get number of canonical/nocel junctions by sample
        print("Making nov_can_DF: " + sampleID)
        nov_can_DF = count_jxns_DF.drop('isoform', axis=1).sum().to_frame().transpose()
        nov_can_DF['sampleID'] = sampleID
        nov_can_DF[exp_factor] = exp_factor_val
    
        # Get UJC ID (jxn hash) for each read. First subset the index file to be just the

        ##Merge classification DF, count jxns DF and jxn hash into one
        class_DF = merge_dfs(class_DF, count_jxns_DF, 'isoform', 'isoform')
    
        # Gene count DF 
        # counting the number of reads in each structural category, in each gene
        gene_category_count_DF = class_DF.pivot_table(index='associated_gene',
                                                      columns='structural_category',
                                                      aggfunc='count',
                                                      fill_value=0)['isoform'].reset_index()
        gene_category_count_DF['total_read_count'] = gene_category_count_DF.iloc[:, 1:].sum(axis=1)
    
        categories = list(class_DF['structural_category'].unique())
    
        ##ISM DF
        if 'incomplete-splice_match' in categories:
            ISM_DF = class_DF[class_DF['structural_category'] == 'incomplete-splice_match'].copy()
            ##New subcategory mono-in-multi
#            ISM_DF.loc[:, 'subcategory'] = np.where(
#                (ISM_DF['flag_ref_monoexon'] == 0) & (ISM_DF['exons'] == 1),
#                'mono_in_multi',
#                ISM_DF['subcategory']
#            )
            ISM_DF['sampleID'] = sampleID
            ISM_DF = ISM_DF.pivot_table(index='sampleID',
                                        columns='subcategory',
                                        aggfunc='count',
                                        fill_value=0)['isoform'].reset_index()
            ISM_DF[exp_factor] = exp_factor_val
    
        ##FSM DF
        if 'full-splice_match' in categories:
            FSM_DF = class_DF[class_DF['structural_category'] == 'full-splice_match'].copy()
            FSM_DF['sampleID'] = sampleID
            FSM_DF = FSM_DF.pivot_table(index='sampleID',
                                        columns='subcategory',
                                        aggfunc='count',
                                        fill_value=0)['isoform'].reset_index()
            FSM_DF[exp_factor] = exp_factor_val
    
        ##NIC/NNC DF
        if 'novel_in_catalog' in categories or 'novel_not_in_catalog' in categories:
            NIC_NNC_DF = class_DF[class_DF['structural_category'].isin(['novel_in_catalog', 'novel_not_in_catalog'])].copy()
            NIC_NNC_DF['sampleID'] = sampleID
            NIC_NNC_DF = NIC_NNC_DF.pivot_table(index='sampleID',
                                                columns='subcategory',
                                                aggfunc='count',
                                                fill_value=0)['isoform'].reset_index()
            NIC_NNC_DF[exp_factor] = exp_factor_val
        
    
        # counting the number of UJCs in each read
        gene_UJC_count_DF = class_DF.groupby('associated_gene')['jxnHash'].nunique().reset_index(name='unique_jxnHash_counts')
        gene_count_DF = pd.merge(gene_category_count_DF, gene_UJC_count_DF, how='outer', on='associated_gene')
        gene_count_DF['sampleID'] = sampleID
        gene_count_DF[exp_factor] = exp_factor_val
    
        ##UJC DF
        ujc_count_DF = class_DF.groupby('jxnHash').agg({
            'isoform': 'nunique',  # Count unique isoforms for each jxnHash
            'associated_gene': 'first',
            'known_canonical': 'first', 
            'known_non_canonical': 'first',
            'novel_canonical': 'first',
            'novel_non_canonical': 'first',
            'structural_category': 'first' }).reset_index()
        ujc_count_DF = ujc_count_DF.fillna(0)
        ujc_count_DF.rename(columns={'isoform': 'read_count'}, inplace=True)
        ujc_count_DF['sampleID'] = sampleID
        ujc_count_DF[exp_factor] = exp_factor_val
    
        # Length DF
        # Calculating the length stats
        total_reads = len(class_DF)
        reads_gt_1kb = (class_DF['length'] > 1000).sum()
        reads_gt_2kb = (class_DF['length'] > 2000).sum()
        reads_gt_3kb = (class_DF['length'] > 3000).sum()
        average_length = class_DF['length'].mean()
        median_length = class_DF['length'].median()
        min_length = class_DF['length'].min()
        max_length = class_DF['length'].max()
        q25_length = class_DF['length'].quantile(0.25)
        q75_length = class_DF['length'].quantile(0.75)
        perc_reads_gt_1kb = (reads_gt_1kb / total_reads) * 100
        perc_reads_gt_2kb = (reads_gt_2kb / total_reads) * 100
        perc_reads_gt_3kb = (reads_gt_3kb / total_reads) * 100
    
        # Creating a new dataframe with length stats
        length_summary_DF = pd.DataFrame({
            'total_reads': [total_reads],
            'reads_gt_1kb': [reads_gt_1kb],
            'reads_gt_2kb': [reads_gt_2kb],
            'reads_gt_3kb': [reads_gt_3kb],
            'perc_reads_gt_1kb': [perc_reads_gt_1kb],
            'perc_reads_gt_2kb': [perc_reads_gt_2kb],
            'perc_reads_gt_3kb': [perc_reads_gt_3kb],
            'average_length': [average_length],
            'median_length': [median_length],
            'min_length': [min_length],
            'max_length': [max_length],
            'q25_length': [q25_length],
            'q75_length': [q75_length],
            'sampleID': [sampleID]
        })
        length_summary_DF[exp_factor] = exp_factor_val
        
        ##Length arrays for violin plot
        length_Dct[sampleID] = np.array(class_DF['length'])
    
        ##Err DFs
        # Count RTS, intrapriming and reads with noncan jxns
        # Calculate counts
        num_reads_RTS = (class_DF['RTS_stage'] == True).sum()
        num_reads_intrapriming = (class_DF['perc_A_downstream_TTS'] > 60).sum()
        num_reads_non_can = (class_DF['all_canonical'] == 'non_canonical').sum()
    
        # Calculate percentages
        total_reads = len(class_DF)
        perc_reads_RTS = (num_reads_RTS / total_reads) * 100
        perc_reads_intrapriming = (num_reads_intrapriming / total_reads) * 100
        perc_reads_non_can = (num_reads_non_can / total_reads) * 100
    
        # Create err DataFrame
        err_DF = pd.DataFrame({
            'num_reads_RTS': [num_reads_RTS],
            'perc_reads_RTS': [perc_reads_RTS],
            'num_reads_intrapriming': [num_reads_intrapriming],
            'perc_reads_intrapriming': [perc_reads_intrapriming],
            'num_reads_non-canonical': [num_reads_non_can],
            'perc_reads_non-canonical': [perc_reads_non_can]
            
        })
        err_DF['sampleID'] = sampleID
        err_DF[exp_factor] = exp_factor_val
    
        ##Calculate junction cv for each of the samples
        cv_DF = calc_jxn_cv(jxn_DF, class_DF, ref_DF, dropFlag=True)
        cv_DF['sampleID'] = sampleID
        cv_DF[exp_factor] = exp_factor_val
    
        ##Store all summary dataframes in dictionaries
        gene_count_dfs[sampleID] = gene_count_DF
        ujc_count_dfs[sampleID] = ujc_count_DF
        length_dfs[sampleID] = length_summary_DF
        cv_dfs[sampleID] = cv_DF
        err_dfs[sampleID] = err_DF
        fsm_dfs[sampleID] = FSM_DF
        ism_dfs[sampleID] = ISM_DF
        nic_nnc_dfs[sampleID] = NIC_NNC_DF
        nov_can_dfs[sampleID] = nov_can_DF
    
        print(sampleID + " done processing")
        del(cv_DF)
        del(class_DF)
        del(jxn_DF)
        print("Memory cleared for next sample")  
        
    return(ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct )

def prep_tables(ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct ):
    
    if args.inFACTOR == None:
        exp_factor = 'temp_factor'
    else:
        exp_factor = args.inFACTOR
    
    
    #Combine datframes for all samples
    #Cat together all gene count dfs
    gene_count_DF = pd.concat(gene_count_dfs.values(), sort=False)
    
    # Add flags to gene count DF
  
    gene_count_DF['flag_annotated_gene'] = np.where(gene_count_DF['associated_gene'].isin(ref_DF['gene_id'].unique()), 1, 0)
    
    #Cat together all ujc dfs
    ujc_count_DF = pd.concat(ujc_count_dfs.values())
    
    ujc_count_DF['flag_annotated_gene'] = np.where(ujc_count_DF['associated_gene'].isin(ref_DF['gene_id'].unique()), 1, 0)
    
    ##Cat together all length DFs 
    length_DF = pd.concat(length_dfs.values(), sort=False)
    cols = ['sampleID', exp_factor] + [col for col in length_DF.columns if col not in ['sampleID', exp_factor]]
    length_DF = length_DF[cols]
    
    #Cat together all err DFs
    
    err_DF = pd.concat(err_dfs.values(), sort=False)
    cols = ['sampleID', exp_factor] + [col for col in err_DF.columns if col not in ['sampleID', exp_factor]]
    err_DF = err_DF[cols]
    ##Cat together all cv DFs
    
    cv_DF = pd.concat(cv_dfs.values(), sort=False)
    
    #Cat subcategory DFs
    
    FSM_DF = pd.concat(fsm_dfs.values(), sort=False)
    FSM_DF.fillna(0, inplace=True)
    
    ISM_DF = pd.concat(ism_dfs.values(), sort=False)
    ISM_DF.fillna(0, inplace=True)
    
    NIC_NNC_DF = pd.concat(nic_nnc_dfs.values(), sort=False)
    NIC_NNC_DF.fillna(0, inplace=True)
    
    #Cat nov_can_dfs
    nov_can_DF = pd.concat(nov_can_dfs.values(), sort=False)
    nov_can_DF.fillna(0, inplace=True)
        
    ##Export tables
    ##!!! Only exporting gene count DF, ujc count DF, length DF and cv_DF for now...
    ##!!! Others can also be exported
    
    if args.inFACTOR == None:
        
        gene_count_DF_drop = gene_count_DF.drop(columns=[exp_factor])
        gene_count_DF_drop.to_csv(os.path.join(args.OUT, args.PREFIX + '_gene_counts.csv'), index=False)
        
        ujc_count_DF_drop = ujc_count_DF.drop(columns=[exp_factor])
        ujc_count_DF_drop.to_csv(os.path.join(args.OUT, args.PREFIX + '_ujc_counts.csv'), index=False)
        
        length_DF_drop = length_DF.drop(columns=[exp_factor])
        length_DF_drop.to_csv(os.path.join(args.OUT, args.PREFIX + '_length_summary.csv'), index=False)
        
        cv_DF_drop = cv_DF.drop(columns=[exp_factor])
        cv_DF_drop.to_csv(os.path.join(args.OUT, args.PREFIX + '_cv.csv'), index=False)
    else:    
        gene_count_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_gene_counts.csv'), index=False)
        ujc_count_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_ujc_counts.csv'), index=False)
        length_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_length_summary.csv'), index=False)
        cv_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_cv.csv'), index=False)
    
    return (gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct)


def identify_cand_underannot(gene_count_DF):

    #Subset gene count DF to only annotated genes
    annot_novel_DF  = gene_count_DF.loc[gene_count_DF['flag_annotated_gene'] == 1, 
                                      ['associated_gene', 'novel_in_catalog', 'novel_not_in_catalog', 
                                       'total_read_count', 'sampleID']].copy()
    
    #Calculate percent novel (NIC+NNC)
    annot_novel_DF['perc_novel'] = ((annot_novel_DF['novel_in_catalog']+annot_novel_DF['novel_not_in_catalog'])/annot_novel_DF['total_read_count'])*100
    
    num_samples = gene_count_DF['sampleID'].nunique()
    perc_novel = args.PERCNOVEL
    annot_exp = args.ANNOTEXP
    #perc_novel = 90
    #perc_sample = 70
    #annot_exp = 50
    
    #Calculate percent samples NIC/NNC > perc_novel and expression above annot_exp
    novel_DF = annot_novel_DF[(annot_novel_DF['perc_novel'] >perc_novel) &
                              (annot_novel_DF['total_read_count'] > annot_exp)].groupby('associated_gene')['sampleID'].nunique().reset_index(name='num_samps_novel_perc_gt_' + str(perc_novel))
    
    novel_DF['perc_samps_novel_perc_gt_' + str(perc_novel)] = (novel_DF['num_samps_novel_perc_gt_' + str(perc_novel)]/num_samples)*100
    
    ##Output underannotated gene table
    novel_DF.to_csv(os.path.join(args.OUT, args.PREFIX + '_cand_underannotated_gene_table.csv'), index=False)
    

def prep_data_4_plots(gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct):
    
    #Map structureal categories to their abbreviations
    abbr_mapping = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
    "antisense": "AS",
    "fusion": "FUS",
    "genic": "GENIC",
    "genic_intron": "GI",
    "intergenic": "INTER"
    }
    
    #Set exp_factor
    if args.inFACTOR == None:
        exp_factor = 'temp_factor'
    else:
        exp_factor= args.inFACTOR
    
    #Make DF with sampleID and correspondign exp_factor values
    exp_factor_DF = gene_count_DF[['sampleID', exp_factor]].drop_duplicates()
    
    ##Prep data for plots 0 - 3 - Vertical dot plot and barplot  %structural categories (all genes and annotated genes)
    total_reads_per_sample = gene_count_DF.groupby('sampleID')['total_read_count'].sum()
    all_gene_percs = {}

    all_categories = [col for col in gene_count_DF.columns if col in abbr_mapping.keys()]
    
    # Loop through each category
    for category in all_categories:
        
        # Sum the category values for annotated genes grouped by sampleID
        all_category_sum_per_sample = gene_count_DF.groupby('sampleID')[category].sum()
        
        # Calculate the percentage of total reads per sample
        all_gene_percs[f'percent_{category}_annotated'] = (all_category_sum_per_sample / total_reads_per_sample) * 100
    
    # Convert the percentages dictionary to a DataFrame 
    all_gene_percs_DF = pd.DataFrame(all_gene_percs)
    all_gene_percs_DF.reset_index(inplace=True)
    
    all_gene_percs_DF.columns = ['sampleID'] + [abbr_mapping[cat] for cat in all_categories]
    all_gene_percs_long_DF = all_gene_percs_DF.melt(id_vars='sampleID', var_name='Category', value_name='Percentage')
    all_gene_percs_long_DF = merge_dfs(all_gene_percs_long_DF,exp_factor_DF, 'sampleID', 'sampleID')
    all_gene_percs_pivot_DF = all_gene_percs_long_DF.pivot_table(index=['sampleID', exp_factor], columns='Category', values='Percentage', fill_value=0).reset_index()

    # Identify annotated genes
    annotated_gene_DF = gene_count_DF[gene_count_DF['flag_annotated_gene'] == 1].copy()
    
    # Aggregate total reads per sample
    total_reads_in_annotated_gene_DF =  annotated_gene_DF.groupby('sampleID')['total_read_count'].sum()
    
    #Initialize gene perc dict
    annot_gene_percs = {}
    
    # Categories to calculate percentages for, in annotated genes
    annot_categories = ['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog','genic','genic_intron']
    
    # Loop through each category
    for category in annot_categories:
    
        if category in annotated_gene_DF.columns:
            # Sum the category values for annotated genes grouped by sampleID
            annot_category_sum_per_sample = annotated_gene_DF.groupby('sampleID')[category].sum()
            
            # Calculate the percentage of total reads per sample
            annot_gene_percs[f'{category}'] = (annot_category_sum_per_sample / total_reads_in_annotated_gene_DF) * 100
        
    # Convert the percentages dictionary to a DataFrame 
    annot_gene_percs_DF = pd.DataFrame(annot_gene_percs)
    annot_gene_percs_DF.reset_index(inplace=True)
    
    #Set column names
    annot_gene_percs_DF.columns = ['sampleID'] + [abbr_mapping[cat] for cat in annot_categories if category in annotated_gene_DF.columns] 
        
    annot_gene_percs_long_DF = annot_gene_percs_DF.melt(id_vars='sampleID', var_name='Category', value_name='Percentage')
    annot_gene_percs_long_DF = merge_dfs(annot_gene_percs_long_DF,exp_factor_DF, 'sampleID', 'sampleID')
    annot_gene_percs_pivot_DF = annot_gene_percs_long_DF.pivot_table(index=['sampleID', exp_factor], columns='Category', values='Percentage', fill_value=0).reset_index()

    #Prep data for plots 10 and 11 - Barplot - number of genes per sample, colured by number of reads assigned to gene
    annotated_gene_DF['read_category'] = annotated_gene_DF['total_read_count'].apply(categorize_by_readcount)
    
    gene_agg_DF = annotated_gene_DF.groupby(['sampleID', 'read_category'])['associated_gene'].nunique().unstack(fill_value=0)
    gene_agg_DF = merge_dfs(gene_agg_DF,exp_factor_DF, 'sampleID', 'sampleID')
    
    # Calculate total counts by sampleID to use for percentage calculation
    gene_counts_by_sample = annotated_gene_DF.groupby(['sampleID', 'read_category'])['associated_gene'].nunique()
    total_gene_counts = annotated_gene_DF.groupby('sampleID')['associated_gene'].nunique()
    
    # Calculate percentages
    gene_percs_read_cat_DF = gene_counts_by_sample.div(total_gene_counts, level='sampleID') * 100  # level='sampleID' ensures division is done within each sampleID group
    # Unstack for plotting
    gene_percs_unstacked = gene_percs_read_cat_DF.unstack(fill_value=0)
    gene_percs_unstacked = merge_dfs(gene_percs_unstacked,exp_factor_DF, 'sampleID', 'sampleID')
    
    ##Prep subcategory dataframes
    
    ## FSM plots 4 & 5
    FSM_DF = FSM_DF.copy()
    
    FSM_perc_DF = FSM_DF.copy()

    # Calculate the sum for the subcategories (excluding 'sampleID' and exp_factor)
    categories = [col for col in FSM_perc_DF.columns if col not in ['sampleID', exp_factor]]
    FSM_perc_DF['total'] = FSM_perc_DF[categories].sum(axis=1)
    
    # Calculate the percentage for each subcategory
    for category in categories:
        FSM_perc_DF[category] = (FSM_perc_DF[category] / FSM_perc_DF['total']) * 100
    FSM_perc_DF.drop('total', axis=1, inplace=True)
    FSM_perc_DF['sampleID'] = FSM_DF['sampleID']
    FSM_perc_DF[exp_factor] = FSM_DF[exp_factor]
    
    ##ISM plots 6 & 7
    ISM_DF = ISM_DF.copy()
    ISM_perc_DF = ISM_DF.copy()

    # Calculate the sum for the subcategories (excluding 'sampleID' and exp_factor)
    categories = [col for col in ISM_perc_DF.columns if col not in ['sampleID', exp_factor]]
    ISM_perc_DF['total'] = ISM_perc_DF[categories].sum(axis=1)

    # Calculate the percentage for each subcategory
    for category in categories:
        ISM_perc_DF[category] = (ISM_perc_DF[category] / ISM_perc_DF['total']) * 100
    ISM_perc_DF.drop('total', axis=1, inplace=True)
    ISM_perc_DF['sampleID'] = ISM_DF['sampleID']
    ISM_perc_DF[exp_factor] = ISM_DF[exp_factor]

    ##NIC_NNC plots 8 & 9
    NIC_NNC_DF = NIC_NNC_DF.copy()
    
    NIC_NNC_perc_DF = NIC_NNC_DF.copy()

    # Calculate the sum for the subcategories (excluding 'sampleID' and exp_factor)
    categories = [col for col in NIC_NNC_perc_DF.columns if col not in ['sampleID', exp_factor]]
    NIC_NNC_perc_DF['total'] = NIC_NNC_perc_DF[categories].sum(axis=1)
    
    # Calculate the percentage for each subcategory
    for category in categories:
        NIC_NNC_perc_DF[category] = (NIC_NNC_perc_DF[category] / NIC_NNC_perc_DF['total']) * 100
    NIC_NNC_perc_DF.drop('total', axis=1, inplace=True)
    NIC_NNC_perc_DF['sampleID'] = NIC_NNC_DF['sampleID']
    NIC_NNC_perc_DF[exp_factor] = NIC_NNC_DF[exp_factor]
    
    
    
    ##Prep data for plots 12-15 -  Boxplots Distribution of % structural category (FSM ISM NIC NNC) 
    #Only on annotated genes
    for category in ['full-splice_match','incomplete-splice_match','novel_in_catalog', 'novel_not_in_catalog']:
        annotated_gene_DF[f'percent_{category}'] = (annotated_gene_DF[category] / annotated_gene_DF['total_read_count']) * 100

    melted_annotated_gene_DF=  annotated_gene_DF.melt(id_vars=['sampleID'], 
                                         value_vars=['percent_full-splice_match','percent_incomplete-splice_match','percent_novel_in_catalog', 'percent_novel_not_in_catalog'],
                        var_name='category', value_name='percentage')
    
    melted_annotated_gene_DF = merge_dfs(melted_annotated_gene_DF,exp_factor_DF, 'sampleID', 'sampleID')
    
    ###UJCs
    
    ##Prep data for  plots 16 - 19 - #Barplot - numberand % UJCs per sample, colured by stackby
     # Add flags to ujc count DF
    
    ujc_count_DF['read_category'] = ujc_count_DF['read_count'].apply(categorize_by_readcount)
    annot_ujc_count_DF = ujc_count_DF[ujc_count_DF['flag_annotated_gene'] == 1]
    
    ujc_cnts_dct = {}
    ujc_percs_dct = {}
    
    for stack_by in ['read_category', 'structural_category']:
        # Group and count unique 'jxnHash' within each 'sampleID' and 'stack_by' category
        ujc_agg_DF = annot_ujc_count_DF.groupby(['sampleID', stack_by])['jxnHash'].nunique().unstack(fill_value=0).reset_index()
        if stack_by == 'structural_category':
           ujc_agg_DF.columns = ['sampleID'] + [abbr_mapping[cat] for cat in ujc_agg_DF.columns[1:] if cat in abbr_mapping]
        
        # Merge with exp_factor_DF 
        ujc_agg_DF = ujc_agg_DF.merge(exp_factor_DF, on='sampleID', how='left')

        # Calculate total unique 'jxnHash' counts for each 'sampleID'
        total_ujc_counts_DF = annot_ujc_count_DF.groupby('sampleID')['jxnHash'].nunique().reset_index()

        # Create an empty DataFrame for percentages
        ujc_percs_DF = pd.DataFrame()
        
       
    # Iterate over columns to calculate percentages, excluding 'sampleID' and any columns from 'exp_factor_DF'
        for column in ujc_agg_DF.columns:
            if column not in ['sampleID'] + list(exp_factor_DF.columns):
                # Align and divide by total counts using a mapping
                total_counts_map = total_ujc_counts_DF.set_index('sampleID')['jxnHash'].to_dict()
                ujc_agg_DF['total_counts'] = ujc_agg_DF['sampleID'].map(total_counts_map)
                ujc_percs_DF[column] = (ujc_agg_DF[column] / ujc_agg_DF['total_counts']) * 100
                
        ujc_percs_DF['sampleID'] = ujc_agg_DF['sampleID']
        
        
        ujc_total_DF = ujc_agg_DF[['total_counts','sampleID']]
        ujc_agg_DF.drop('total_counts', axis=1, inplace=True)
        ujc_percs_DF = merge_dfs(ujc_percs_DF,exp_factor_DF, 'sampleID', 'sampleID')
       
        # Store results in dictionaries
        ujc_cnts_dct[stack_by] = ujc_agg_DF
        ujc_percs_dct[stack_by] = ujc_percs_DF

    ##Lengths
    ##Prep data for plot 21 - - Barplot - number of reads caoloured by read count category
    # Calculate counts for each category
    length_DF['reads_lt_1kb'] = length_DF['total_reads'] - length_DF['reads_gt_1kb']
    length_DF['reads_1kb_to_2kb'] = length_DF['reads_gt_1kb'] - length_DF['reads_gt_2kb']
    length_DF['reads_2kb_to_3kb'] = length_DF['reads_gt_2kb'] - length_DF['reads_gt_3kb']
    
    # Calculate percentages for each category
    length_cols = ['reads_lt_1kb', 'reads_1kb_to_2kb', 'reads_2kb_to_3kb', 'reads_gt_3kb']
    for category in length_cols:
       length_DF[f'{category}_perc'] = (length_DF[category] / length_DF['total_reads']) * 100
    # Aggregating the counts into a new DataFrame suitable for plotting
    length_cnts_agg = length_DF.set_index('sampleID')[length_cols]
    length_cnts_agg = merge_dfs(length_cnts_agg,exp_factor_DF, 'sampleID', 'sampleID')
    
    ##Prep data for plot 22 - Barplot % reads by read count category
    # Calculating percentages for plotting
    percent_length_cols = [f'{col}_perc' for col in length_cols]
    length_percs_agg = length_DF.set_index('sampleID')[percent_length_cols]
    length_percs_agg = merge_dfs(length_percs_agg,exp_factor_DF, 'sampleID', 'sampleID')
    
    #Prep data for plot  and 24-29 - % structural category vs %reads greater than 1kb
    length_DF= merge_dfs(length_DF,all_gene_percs_DF, 'sampleID', 'sampleID')
    
    #Prep data for plot 23
    ##Length for violin plots
    length_DF2 = pd.DataFrame({k: pd.Series(v) for k, v in length_Dct.items()})
    length_DF2 = length_DF2.melt(var_name='sampleID', value_name='length')
    length_DF2 = merge_dfs(length_DF2,exp_factor_DF, 'sampleID', 'sampleID')
    
    
    ##Prep data for PCA - plot 32 + 33
    sample_data_DF = merge_three_dfs(length_DF, err_DF.drop(exp_factor, axis=1), ujc_total_DF, 'sampleID', 'sampleID', 'sampleID')
    pca_features = sample_data_DF.columns.difference(['sampleID', exp_factor])
    scaler = StandardScaler()    
    scaled_feature_DF = scaler.fit_transform(sample_data_DF[pca_features])
    pca = PCA()
    pca_results = pca.fit_transform(scaled_feature_DF)
    pca_DF = pd.DataFrame(data = pca_results)
    pca_DF['sampleID'] = sample_data_DF['sampleID']
    pca_DF[exp_factor] = sample_data_DF[exp_factor]
    
    #Get PCA loadings
    loadings = pca.components_.T
    loadings_DF = pd.DataFrame(data=loadings, index=pca_features)
    
    #Get variance ratios
    variance_ratio = pca.explained_variance_ratio_ 
     
    ##Prep for cv plots 36 & 37
    #jxn_exp=3
    jxn_exp = args.JXNEXP
    cv_don_DF = cv_DF[(cv_DF['count'] >= jxn_exp) & (cv_DF['flag_donor'] == 1)]
    cv_acc_DF = cv_DF[(cv_DF['count'] >= jxn_exp) & (cv_DF['flag_acceptor'] == 1)]
    
    if exp_factor == 'temp_factor':
        cv_don_DF = cv_don_DF.copy()
        cv_acc_DF = cv_acc_DF.copy()
        cv_don_DF.loc[:, exp_factor] = cv_don_DF[exp_factor].fillna(0)
        cv_acc_DF.loc[:, exp_factor] = cv_acc_DF[exp_factor].fillna(0)
        
  
    cv_don_summary = cv_don_DF.groupby(['sampleID',exp_factor]).apply(
            lambda x: pd.Series({
                    'ref_match': (x['mean_abs_diff'] == 0).sum(),
                    'cv_0': ((x['cv'] == 0) & (x['mean_abs_diff'] != 0)).sum(),
                    'cv_gt_0': (x['cv'] > 0).sum()
                    })
                    ).reset_index()
    cv_don_summary.fillna(0, inplace=True)
    
    cv_acc_summary = cv_acc_DF.groupby(['sampleID',exp_factor]).apply(
            lambda x: pd.Series({
                    'ref_match': (x['mean_abs_diff'] == 0).sum(),
                    'cv_0': ((x['cv'] == 0) & (x['mean_abs_diff'] != 0)).sum(),
                    'cv_gt_0': (x['cv'] > 0).sum()
                    })
                    ).reset_index()
    cv_acc_summary.fillna(0, inplace=True)

    #Fill zeroes for samples that have no acceptors or donors meeting threshold 
    for sampleID in exp_factor_DF['sampleID']:
        if sampleID not in cv_don_summary['sampleID'].values:
            row = {column: 0 for column in cv_don_summary.columns}
            row['sampleID'] = sampleID  # Set the sampleID for the new row
            row[exp_factor] = exp_factor_DF.loc[exp_factor_DF['sampleID'] == sampleID, exp_factor].iloc[0]
            row = pd.DataFrame([row])
            cv_don_summary = pd.concat([cv_don_summary, row], ignore_index=True)
            print(f"Note: {sampleID} has no donors in annotated genes with > {jxn_exp} reads")
            
        if sampleID not in cv_acc_summary['sampleID'].values:
            row = {column: 0 for column in cv_acc_summary.columns}
            row['sampleID'] = sampleID  # Set the sampleID for the new row
            row[exp_factor] = exp_factor_DF.loc[exp_factor_DF['sampleID'] == sampleID, exp_factor].iloc[0]
            row = pd.DataFrame([row])
            cv_acc_summary = pd.concat([cv_acc_summary, row], ignore_index=True)
            print(f"Note: {sampleID} has no acceptors in annotated genes with > {jxn_exp} reads")

    
    
    #Make nov_can_perc_DF - plot 35
    
    nov_can_perc_DF = nov_can_DF.copy()
    categories = [col for col in nov_can_perc_DF.columns if col not in ['sampleID', exp_factor]]
    nov_can_perc_DF['total'] = nov_can_perc_DF[categories].sum(axis=1)
    for category in categories:
        nov_can_perc_DF[category] = (nov_can_perc_DF[category] / nov_can_perc_DF['total']) * 100
    nov_can_perc_DF.drop('total', axis=1, inplace=True)
    nov_can_perc_DF['sampleID'] = nov_can_DF['sampleID']
    nov_can_perc_DF[exp_factor] = nov_can_DF[exp_factor]
    
    return (all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF, annot_gene_percs_pivot_DF, gene_agg_DF, 
             gene_percs_unstacked, melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF, 
             length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio, 
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_NNC_perc_DF, nov_can_DF, nov_can_perc_DF, length_DF2)
    
    
    
def plot_pdf_by_factor(out_path, all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF, annot_gene_percs_pivot_DF, gene_agg_DF, 
             gene_percs_unstacked, melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF, 
             length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio, 
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_NNC_perc_DF,nov_can_DF, nov_can_perc_DF,length_DF2 ):
    
    plt.rcParams.update({'font.size': 12})
    
    exp_factor = args.inFACTOR

    def plot_stacked_bars(*args, **kwargs):
        data = kwargs.pop('data')
        ax = plt.gca()
        bottom = np.zeros(len(data))
        for idx, category in enumerate(categories):
            values = data[category].values
            # Only add bars for non-zero values
            non_zero_indices = values != 0
            ax.bar(data['sampleID'][non_zero_indices], values[non_zero_indices], bottom=bottom[non_zero_indices], color=palette[idx], label=category)
            # Update bottom only with non-zero values
            bottom += values

    def plot_stacked_bars_custom_palette(data, color_palette, *args, **kwargs):
        ax = plt.gca()
        bottom = np.zeros(len(data))
        
        kwargs.pop('color', None)
        categories = [col for col in data.columns if col not in ['sampleID', exp_factor]]
        for category in categories:
            values = data[category].values
            if category in color_palette:
                # Only add bars for non-zero values
                non_zero_indices = values != 0
                ax.bar(data['sampleID'][non_zero_indices], values[non_zero_indices], bottom=bottom[non_zero_indices], color=color_palette[category], label=category, **kwargs)
                # Update bottom only with non-zero values
                bottom += values
            else:
                print(f"Color for {category} not found in palette.")
 
    def plot_side_by_side_bars(*args, **kwargs):
        data = kwargs.pop('data')
        ax = plt.gca()
        bar_width = 0.35
        num_samples = len(data['sampleID'].unique())
        gap_width = 0.1  # Width of the gap between groups
        total_bar_width = (len(categories) * bar_width) + gap_width
    
        # Calculate initial positions for each group, ensuring gaps between groups
        positions = np.arange(num_samples) * total_bar_width
    
        for idx, category in enumerate(categories):
            # Offset positions within each group for each category
            category_positions = positions + idx * bar_width
            values = data[category].values
            # Filter out zero values to maintain visualization integrity
            non_zero_indices = values != 0
            ax.bar(category_positions[non_zero_indices], values[non_zero_indices], width=bar_width, color=palette[idx], label=category)
    
        # Centralize the x-ticks for each group
        central_offset = bar_width * (len(categories) - 1) / 2
        ax.set_xticks(positions + central_offset)
        ax.set_xticklabels(data['sampleID'].unique(), rotation=90)
        
    num_factors = all_gene_percs_long_DF[exp_factor].nunique()
    num_samples = all_gene_percs_long_DF['sampleID'].nunique()
    
    #Define category color palette
    
    category_color_palette = {
    "FSM": "#6BAED6",
    "ISM": "#FC8D59",
    "NIC": "#78C679",
    "NNC": "#EE6A50",
    "GENIC": "#969696",
    "AS": "#66C2A4",
    "FUS": "#ffc125",  
    "INTER": "#e9967a",  
    "GI": "#41B6C4"
    }   
    
    subcat_color_palette = {
    "alternative_5end": '#02314d',
    "alternative_3end": '#0e5a87',
    "alternative_3end5end": '#7ccdfc',
    'reference_match': '#c4e1f2',
    "3prime_fragment": '#c4531d',
    "internal_fragment": '#e37744',
    "5prime_fragment": '#e0936e',
    "combination_of_known_junctions": '#014d02',
    "combination_of_known_splicesites": '#379637',
    "intron_retention": '#81eb82',
    "at_least_one_novel_splicesite": '#6ec091', ## Changed this from Not comb. of annot. junctions
    "mono-exon_by_intron_retention": '#4aaa72',
    "At least 1 annot. don./accept.": '#32734d',
    "mono-exon": '#cec2d2',
    "multi-exon": '#876a91',
    "mono_in_multi": '#aec6cf'
    }
    
    #Define sample color palette
    
    unique_sampleIDs = all_gene_percs_long_DF['sampleID'].unique()
    sample_colors = plt.cm.rainbow(np.linspace(0, 1, num_samples))
    sample_color_palette = {sampleID: color for sampleID, color in zip(unique_sampleIDs, sample_colors)}
    
    with PdfPages(out_path) as pdf:
        ### Plot 0 - Vertical Scatterplot -Structural Category - All genes
        
        palette = sample_color_palette
        g = sns.catplot(
            data=all_gene_percs_long_DF,
            x='Category',
            y='Percentage',
            hue='sampleID',
            col=exp_factor,  
            kind='strip',
            col_wrap=num_factors, 
            jitter=False,  
            palette=palette,  
            linewidth=0,
            size=20,  
            alpha=0.6, 
            legend=True,
            height=8, 
            aspect=1.6 
        )
        g.set_xticklabels(rotation=0)
        g.set_axis_labels("Structural Category", "Percentages")
        g.fig.suptitle("% Reads in Each Structural Category - All Genes", y=1.02, fontsize=16)
        title = g.fig.suptitle("% Reads in Each Structural Category - All Genes", y=1.02, fontsize=16)
        legend = g._legend
        legend.set_frame_on(True)
        legend.get_frame().set_facecolor((1, 1, 1, 0.5))
        #handles, labels = g.fig.axes[0].get_legend_handles_labels()
        #lgd = g.fig.legend(handles, labels, title='Sample ID', loc='upper left', bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure)
        plt.subplots_adjust(right=0.8)
        plt.tight_layout(pad=3)
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        ### Plot 1 - Vertical Scatterplot -Structural Category - Annotated genes
        g = sns.catplot(
            data=annot_gene_percs_long_DF,
            x='Category',
            y='Percentage',
            hue='sampleID',
            col=exp_factor,  
            kind='strip',
            col_wrap=num_factors,  
            jitter=False,
            palette=palette,  
            linewidth=0,
            size=20, 
            alpha=0.6, 
            legend=True,
            height=8, 
            aspect=1.6 
        )
        g.set_xticklabels(rotation=0)
        g.set_axis_labels("Structural Category", "Percentages")
        title = g.fig.suptitle("% Reads in Each Structural Category - Annotated Genes", y=1.02, fontsize=16)
        #handles, labels = g.fig.axes[0].get_legend_handles_labels()
        #lgd = g.fig.legend(handles, labels, title='Sample ID', loc='upper left', bbox_to_anchor=(1.05, 1), bbox_transform=plt.gcf().transFigure)
        legend = g._legend
        legend.set_frame_on(True)
        legend.get_frame().set_facecolor((1, 1, 1, 0.5))
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        ## Plot 2 - Barplot - structural category % - ALL genes
        
        all_gene_percs_pivot_DF=all_gene_percs_pivot_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in all_gene_percs_pivot_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(all_gene_percs_pivot_DF, col=exp_factor, col_wrap=num_factors, height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = category_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), all_gene_percs_pivot_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Percentages")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Percentage Reads in Each Structural Category - All Genes", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        
       
        
         ## Plot 3 - Barplot - structural category % - Annotated genes
         
        annot_gene_percs_pivot_DF =annot_gene_percs_pivot_DF.sort_values(by=[exp_factor, 'sampleID'])
        categories = [cat for cat in ['FSM', 'ISM', 'NIC', 'NNC','GI','GENIC'] if cat in annot_gene_percs_pivot_DF.columns]
        g = sns.FacetGrid(annot_gene_percs_pivot_DF, col=exp_factor, col_wrap=num_factors, height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = category_color_palette)
        g.set_axis_labels("Sample ID", "Percentages")
        g.set_titles(exp_factor + " = {col_name}")
        for ax, (name, group) in zip(g.axes.flatten(), annot_gene_percs_pivot_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        title = g.fig.suptitle("Percentage Reads in Each Structural Category - Annotated Genes", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
             ## Barplot subcategories
             #Plot 4
        FSM_DF = FSM_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in FSM_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(FSM_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = subcat_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), FSM_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of reads")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Number of Reads in Each Sub-Category - FSM ", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 5
        FSM_perc_DF =FSM_perc_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in FSM_perc_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(FSM_perc_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = subcat_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), FSM_perc_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Percentage")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Percentage of FSM Reads in Each Sub-Category ", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 6
        ISM_DF =ISM_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in ISM_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(ISM_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = subcat_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), ISM_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of reads")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Number of Reads in Each Sub-Category - ISM ", y=1.02, fontsize=16)
        lgd =  plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 7
        ISM_perc_DF =ISM_perc_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in ISM_perc_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(ISM_perc_DF, col=exp_factor, col_wrap=num_factors, height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = subcat_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), ISM_perc_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Percentage")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Percentage of ISM Reads in Each Sub-Category ", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 8
        NIC_NNC_DF =NIC_NNC_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in NIC_NNC_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(NIC_NNC_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = subcat_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), NIC_NNC_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of reads")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Number of Reads in Each Sub-Category - NIC and NNC ", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 9
        NIC_NNC_perc_DF =NIC_NNC_perc_DF.sort_values(by=[exp_factor, 'sampleID']) 
        categories =  [col for col in NIC_NNC_perc_DF.columns if col not in ['sampleID', exp_factor]]   
        g = sns.FacetGrid(NIC_NNC_perc_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6,sharex=False, sharey=True)
        g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = subcat_color_palette)
        for ax, (name, group) in zip(g.axes.flatten(), NIC_NNC_perc_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Percentage")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Percentage of NIC/NNC Reads in Each Sub-Category ", y=1.02, fontsize=16)
        lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        
        ## Plot 10 Barplot - Genes detected counts
        gene_agg_DF = gene_agg_DF.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['1 read','2-10 reads', '11-50 reads', '50-100 reads','100+ reads']
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(gene_agg_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_stacked_bars)
        g.set_axis_labels("Sample ID", "Number of Genes")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Number of Genes Detected", y=1.02, fontsize=16)
        lgd = plt.legend(title='Read Count', bbox_to_anchor=(1.05, 1), loc='upper left')
        for ax, (name, group) in zip(g.axes.flatten(), gene_agg_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        
        #Plot 11 Barplot - Genes detected percentages
        gene_percs_unstacked = gene_percs_unstacked.sort_values(by=[exp_factor, 'sampleID'])
        g = sns.FacetGrid(gene_percs_unstacked, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)    
        
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(), gene_percs_unstacked.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Percentage")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Number of Genes Detected", y=1.02, fontsize=16)
        lgd = plt.legend(title='Read Count', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
    
            
        
            
        # Plots 12 - 15 Boxplots Distribution of % structural category (FSM ISM NIC NNC) 
        melted_annotated_gene_DF = melted_annotated_gene_DF.sort_values(by=[exp_factor, 'sampleID'])
        palette = sample_color_palette
        for category in melted_annotated_gene_DF['category'].unique():
            plt.figure(figsize=(10, 7.5))  
            for i, exp_factor_val in enumerate(melted_annotated_gene_DF[exp_factor].unique(), start=1):
                # Filter the DataFrame for both category and exp_factor
                df_filtered = melted_annotated_gene_DF[(melted_annotated_gene_DF['category'] == category) &
                                                       (melted_annotated_gene_DF[exp_factor] == exp_factor_val)]
                
                # Create a subplot for this exp_factor
                plt.subplot(1, num_factors, i)
                sns.boxplot(data=df_filtered, x='sampleID', y='percentage', palette=palette)
                plt.title(exp_factor + f" = {exp_factor_val}")
                plt.xticks(rotation=90)
            
            title = plt.suptitle(f'Gene distribution - {category}', y= 1.02)
            plt.tight_layout()
            pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
            plt.show()
            plt.close()
        
    
    
        # Plots 16 -19  - UJC barplots - percent and counts, structural 
        for stack_by in ['read_category','structural_category']:
            
            ujc_cnts_dct[stack_by] =  ujc_cnts_dct[stack_by].sort_values(by=[exp_factor, 'sampleID'])
            ujc_percs_dct[stack_by] =  ujc_percs_dct[stack_by].sort_values(by=[exp_factor, 'sampleID'])
            
            if stack_by == 'read_category':
                categories = ['1 read','2-10 reads', '11-50 reads', '50-100 reads','100+ reads']  
                lab = 'Read count'
                
                palette = sns.color_palette("rainbow", len(categories))
                g = sns.FacetGrid(ujc_cnts_dct[stack_by], col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
                
                #Barplot - number of UJCs per sample, colured by stackby
                
                g.map_dataframe(plot_stacked_bars)
                for ax, (name, group) in zip(g.axes.flatten(), ujc_cnts_dct[stack_by].groupby(exp_factor)):
                    ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
                    ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
                g.set_axis_labels("Sample ID","Number of UJCs" )
                g.set_titles(exp_factor + " = {col_name}")
                title = g.fig.suptitle("Number of UJCs Detected", y=1.02, fontsize=16)
                lgd = plt.legend(title=lab, bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
                plt.show()
                plt.close()
                
                g = sns.FacetGrid(ujc_percs_dct[stack_by], col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
                #Barplot - number of UJCs per sample, colured by stackby
                g.map_dataframe(plot_stacked_bars)
                for ax, (name, group) in zip(g.axes.flatten(),ujc_percs_dct[stack_by].groupby(exp_factor)):
                    ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
                    ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
                g.set_axis_labels("Sample ID", "Percentage of UJCs")
                g.set_titles(exp_factor + " = {col_name}")
                title = g.fig.suptitle("UJCs Detected", y=1.02, fontsize=16)
                lgd = plt.legend(title=lab, bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
                plt.show()
                plt.close()
                
            elif stack_by == 'structural_category':
                categories = [col for col in ujc_cnts_dct[stack_by].columns if col not in ['sampleID', exp_factor]]
                lab = 'Structural category'
                
                g = sns.FacetGrid(ujc_cnts_dct[stack_by], col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
                g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = category_color_palette)
                for ax, (name, group) in zip(g.axes.flatten(), ujc_cnts_dct[stack_by].groupby(exp_factor)):
                    ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
                    ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
                g.set_axis_labels("Sample ID", "Number of UJCs")
                g.set_titles(exp_factor + " = {col_name}")
                title = g.fig.suptitle("Number of UJCs detected", y=1.02, fontsize=16)
                lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
                plt.show()
                plt.close()
                
            
                g = sns.FacetGrid(ujc_percs_dct[stack_by], col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)
                g.map_dataframe(plot_stacked_bars_custom_palette, color_palette = category_color_palette)
                for ax, (name, group) in zip(g.axes.flatten(), ujc_percs_dct[stack_by].groupby(exp_factor)):
                    ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
                    ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
                g.set_axis_labels("Sample ID", "Number of UJCs")
                g.set_titles(exp_factor + " = {col_name}")
                title = g.fig.suptitle("Percentage of UJCs detected", y=1.02, fontsize=16)
                lgd = plt.legend(title='Structural Category', bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
                plt.show()
                plt.close()
                
                
                
                
        ## Plot 20  Total Mapped reads vs % reads gt 1kb
    
        palette = sample_color_palette
        plt.figure(figsize=(10, 7.5))
        sns.scatterplot(data=length_DF, y='total_reads', x='perc_reads_gt_1kb', hue='sampleID', style=exp_factor,
                        palette=palette, legend='full', s=100)
        
        plt.title('Total Reads vs Percentage of Reads > 1kb')
        plt.xlabel('Percentage of Reads > 1kb')
        plt.ylabel('Total Reads')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Legend')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        ## Plot 21 - Bar graph read length categories - counts
        length_cnts_agg = length_cnts_agg.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['reads_lt_1kb','reads_1kb_to_2kb', 'reads_2kb_to_3kb','reads_gt_3kb' ]
    
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(length_cnts_agg, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)    
        
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(), length_cnts_agg.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of Reads")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Lengths of All Mapped Reads", y=1.02, fontsize=16)
        lgd = plt.legend(title='Read Count', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
         
        ## Plot 22 Bar graph read length categories - %
        length_percs_agg = length_percs_agg.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['reads_lt_1kb_perc','reads_1kb_to_2kb_perc', 'reads_2kb_to_3kb_perc','reads_gt_3kb_perc']
        g = sns.FacetGrid(length_percs_agg, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6,sharex=False, sharey=True)    
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(), length_percs_agg.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Percentage")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Lengths of All Mapped Reads", y=1.02, fontsize=16)
        lgd = plt.legend(title='Read Count', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.legend(title='Read Count', bbox_to_anchor=(1.05, 1), loc='upper left')
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 23 - Violin plots of lengths
        g = sns.FacetGrid(length_DF2, col=exp_factor, col_wrap=num_factors, height=8, aspect=1.6, sharex=False, sharey=True)
        g.map_dataframe(sns.violinplot, x='sampleID', y='length', palette=sample_color_palette)
        g.set_axis_labels("Sample ID", "Length")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Read Length Distribution", y=1.02, fontsize=16)
        #handles = [plt.Rectangle((0,0),1,1, color=sample_color_palette[sampleID]) for sampleID in unique_sampleIDs]
        #labels = list(unique_sampleIDs)
        #lgd = g.fig.legend(handles, labels, title='Sample ID', bbox_to_anchor=(1.05, 1), loc='upper left')
        for ax in g.axes.flatten():
            # Set x-tick labels with sample IDs, rotating for better visibility
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.show()
        plt.close()

        ## Plots 24 - 29 Scatterplot % reads > 1kb vs % structural category
        categories =[cat for cat in ['FSM', 'ISM', 'NIC', 'NNC','GI','GENIC'] if cat in length_DF.columns]
        
        palette = sample_color_palette
        for category in categories:
            plt.figure(figsize=(10, 7.5))
            sns.scatterplot(data=length_DF, y='perc_reads_gt_1kb', x=category, 
                        hue='sampleID', style=exp_factor,
                        palette=palette, legend='full', s=100)
        
            plt.title('Percentage of Reads > 1kb vs %' + category)
            plt.xlabel('% ' + category)
            plt.ylabel('Percentage of Reads > 1kb')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Legend')
            plt.tight_layout()
            pdf.savefig()
            plt.show()
            plt.close()
            
        # Plot 30 - Bar graph RTS/intra priming
        err_DF = err_DF.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['num_reads_RTS','num_reads_intrapriming', 'num_reads_non-canonical']
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(err_DF, col=exp_factor, col_wrap=num_factors,  height=8, aspect = 1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_side_by_side_bars)
        for ax, (name, group) in zip(g.axes.flatten(), err_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of Reads")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Number of Artefact Reads", y=1.02, fontsize=16)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
    
        #Plot 31 - Bar graph RTS/intra priming (percentages)
        categories = ['perc_reads_RTS','perc_reads_intrapriming', 'perc_reads_non-canonical']
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(err_DF, col=exp_factor, col_wrap=num_factors, height=8, aspect=1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_side_by_side_bars)

        g.set_axis_labels("Sample ID", "Percentage of Reads")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Percentage of Artefact Reads", y=1.02, fontsize=16)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 32 - PCA plot
        palette = sample_color_palette
        plt.figure(figsize=(10, 7.5))
        sns.scatterplot(x=0, y=1, hue='sampleID', data=pca_DF, 
                        style=exp_factor,palette=palette, s = 100)
        title = plt.title('PCA Plot Based on QC metrics')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Legend')
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
                
        #Plot 33 - Scree plot and loadings heatmap
        # Calculate the cumulative variance and determine the number of components to use
        cumulative_variance = np.cumsum(variance_ratio)
        n_components = np.argmax(cumulative_variance >= 0.85) + 1
        
         # Create the plots
        fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex='col', gridspec_kw={'width_ratios': [10, 3], 'height_ratios': [3, 10]})
        loadings_DF = loadings_DF.iloc[:, :n_components]
        link = linkage(loadings_DF, method='average')
        sorted_idx = leaves_list(link)
        loadings_DF = loadings_DF.iloc[sorted_idx]
        
        # Bar plot for explained variance (Scree Plot)
        x_tick_pos = [i + 0.5 for i in range(n_components)]
        ax[0, 0].bar(x_tick_pos, variance_ratio[:n_components], align='center', label='Individual explained variance')
        ax[0, 0].step(x_tick_pos, cumulative_variance[:n_components], where='mid', label='Cumulative explained variance')
        ax[0, 0].set_xticks(x_tick_pos)
        ax[0, 0].set_xticklabels([])  # Clear x tick labels here
        ax[0, 0].set_ylabel('Variance Explained')
        
        # Set x tick labels for the heatmap
        x_ticks = [f'PC{i+1}' for i in range(n_components)]
        sns.heatmap(loadings_DF, cmap="coolwarm", ax=ax[1, 0], cbar_ax=ax[1, 1], xticklabels=x_ticks)
        ax[1, 0].set_yticks(np.arange(loadings_DF.shape[0]) + 0.5)
        ax[1, 0].set_yticklabels(loadings_DF.index, rotation=0)
        ax[1, 0].set_xlabel('Principal Components')
        
        # Use the ax[0,1] for legend
        ax[0, 1].axis('off')  # Turn off the axis lines and labels
        handles, labels = ax[0, 0].get_legend_handles_labels()
        ax[0, 1].legend(handles, labels, loc='center')  # Place legend at the center of ax[0, 1]
        
        title = fig.suptitle('Screeplot and Heatmap of PC loadings',y=1.02, fontsize=16)
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        
        #Plot 34 - Barplot junction category counts
        nov_can_DF = nov_can_DF.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(nov_can_DF, col=exp_factor, col_wrap=num_factors, height=9, aspect = 1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(),nov_can_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
        g.set_axis_labels("Sample ID", "Number of Junctions")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Junctions by Category", y=1.02, fontsize=16)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 35 - Barplot junction category percentages
        nov_can_perc_DF = nov_can_perc_DF.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(nov_can_perc_DF, col=exp_factor, col_wrap=num_factors, height=9, aspect = 1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(),nov_can_perc_DF.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
        g.set_axis_labels("Sample ID", "Percentage")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Junctions by Category", y=1.02, fontsize=16)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 36 - Barplot cv junction categories acceptors
        cv_acc_summary = cv_acc_summary.sort_values(by=[exp_factor, 'sampleID'])
        categories = ['ref_match','cv_0','cv_gt_0']
        palette = sns.color_palette("rainbow", len(categories))
        g = sns.FacetGrid(cv_acc_summary, col=exp_factor, col_wrap=num_factors, height=8, aspect = 1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(),cv_acc_summary.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of Acceptors")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Acceptors with > 3 reads", y=1.02, fontsize=16)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 37 - Barplot cv junction categories donors
        cv_don_summary = cv_don_summary.sort_values(by=[exp_factor, 'sampleID'])
        g = sns.FacetGrid(cv_don_summary, col=exp_factor, col_wrap=num_factors, height=8, aspect = 1.6, sharex=False, sharey=True)    
        g.map_dataframe(plot_stacked_bars)
        for ax, (name, group) in zip(g.axes.flatten(), cv_don_summary.groupby(exp_factor)):
            ax.set_xticklabels(group['sampleID'].unique(), rotation=90)
            ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        g.set_axis_labels("Sample ID", "Number of Donors")
        g.set_titles(exp_factor + " = {col_name}")
        title = g.fig.suptitle("Donors with > 3 reads", y=1.02, fontsize=16)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(lgd,title,), bbox_inches='tight')
        plt.show()
        plt.close()
    
        
        
    

            
def plot_pdf(out_path, all_gene_percs_long_DF, annot_gene_percs_long_DF, all_gene_percs_pivot_DF, annot_gene_percs_pivot_DF, gene_agg_DF, 
             gene_percs_unstacked, melted_annotated_gene_DF, ujc_cnts_dct, ujc_percs_dct, length_DF, 
             length_cnts_agg, length_percs_agg, err_DF, pca_DF, loadings_DF, variance_ratio, 
             cv_acc_summary, cv_don_summary, FSM_DF, ISM_DF, NIC_NNC_DF, FSM_perc_DF, ISM_perc_DF, NIC_NNC_perc_DF,nov_can_DF, nov_can_perc_DF,length_DF2 ):
    
    
    plt.rcParams.update({'font.size': 11})
    
    
    num_samples = all_gene_percs_long_DF['sampleID'].nunique()
    exp_factor = 'temp_factor'
    
 #Define category color palette
    
    category_color_palette = {
    "FSM": "#6BAED6",
    "ISM": "#FC8D59",
    "NIC": "#78C679",
    "NNC": "#EE6A50",
    "GENIC": "#969696",
    "AS": "#66C2A4",
    "FUS": "#ffc125",  
    "INTER": "#e9967a",  
    "GI": "#41B6C4"
    }   
    
    subcat_color_palette = {
    "alternative_5end": '#02314d',
    "alternative_3end": '#0e5a87',
    "alternative_3end5end": '#7ccdfc',
    'reference_match': '#c4e1f2',
    "3prime_fragment": '#c4531d',
    "internal_fragment": '#e37744',
    "5prime_fragment": '#e0936e',
    "combination_of_known_junctions": '#014d02',
    "combination_of_known_splicesites": '#379637',
    "intron_retention": '#81eb82',
    "at_least_one_novel_splicesite": '#6ec091', ## Changed this from Not comb. of annot. junctions
    "mono-exon_by_intron_retention": '#4aaa72',
    "At least 1 annot. don./accept.": '#32734d',
    "mono-exon": '#cec2d2',
    "multi-exon": '#876a91',
    "mono_in_multi": '#aec6cf'
    }
    
    #Define sample color palette
    
    unique_sampleIDs = all_gene_percs_long_DF['sampleID'].unique()
    sample_colors = plt.cm.rainbow(np.linspace(0, 1, num_samples))
    sample_color_palette = {sampleID: color for sampleID, color in zip(unique_sampleIDs, sample_colors)}
    
    
    with PdfPages(out_path) as pdf:
        
        ##Plot 0 - Vertical plot - reads in all genes by structural category
        ### each dot is a sample
        palette = sample_color_palette
        plt.figure(figsize=(10, 7.5))
        sns.stripplot(x='Category', y='Percentage', data=all_gene_percs_long_DF, jitter=0, 
                      alpha=0.6, hue='sampleID', dodge=False, **{'linewidth': 0, 's': 20}, palette = palette)
        plt.xticks(rotation=45)
        plt.ylabel('Percentage')
        plt.title(' % reads in each structural category')
        plt.legend(title='SampleID', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        ##Plot 1 - Vertical plot - reads in annotated genes by structural category
        plt.figure(figsize=(10, 7.5))
        sns.stripplot(x='Category', y='Percentage', data=annot_gene_percs_long_DF, jitter=0, 
                      alpha=0.6, hue='sampleID', dodge=False, **{'linewidth': 0, 's': 20}, palette = palette)
        plt.xticks(rotation=45)
        plt.ylabel('Percentage')
        plt.title(' % reads in each structural category')
        plt.legend(title='SampleID', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        ##Plot 2 - Barplot - reads in all genes by structural category
        
        categories = [col for col in all_gene_percs_pivot_DF.columns if col not in ['sampleID', exp_factor]]
        
        cols = ['sampleID'] + categories
        all_gene_percs_pivot_DF = all_gene_percs_pivot_DF[cols]
        all_gene_percs_pivot_DF = all_gene_percs_pivot_DF.sort_values(by= 'sampleID')
        plt.figure(figsize=(10, 7.5))
        bottom = np.zeros(len(all_gene_percs_pivot_DF))

        for col in categories:
            if col in category_color_palette:
                plt.bar(all_gene_percs_pivot_DF.index, all_gene_percs_pivot_DF[col], bottom=bottom, label=col, color=category_color_palette[col])
                bottom += all_gene_percs_pivot_DF[col].values

        plt.title('Percentage Reads in Each Structural Category - All Genes')
        plt.xlabel('SampleID')
        plt.ylabel('Percentages')
        plt.xticks(ticks=np.arange(len(all_gene_percs_pivot_DF['sampleID'])), labels=all_gene_percs_pivot_DF['sampleID'], rotation=90, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig() 
        plt.show()
        plt.close()  
        
        ##Plots 3 - Barplot - reads in annotated genes by structural category
        categories = [cat for cat in ['FSM', 'ISM', 'NIC', 'NNC', 'GI', 'GENIC'] if cat in annot_gene_percs_pivot_DF.columns]
        palette = [category_color_palette[cat] for cat in categories]  
        annot_gene_percs_pivot_DF = annot_gene_percs_pivot_DF.sort_values(by= 'sampleID')
        plt.figure(figsize=(10, 7.5))
        annot_gene_percs_pivot_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Percentage Reads in Each Structural Category - Annotated Genes')
        plt.xlabel('SampleID')
        plt.ylabel('Percentages')
        plt.xticks(ticks=np.arange(len(annot_gene_percs_pivot_DF['sampleID'])), labels=annot_gene_percs_pivot_DF['sampleID'], rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        ##Subcategory plots
        ##Plot 4
        FSM_DF.set_index('sampleID', inplace=True)
        categories =  [col for col in FSM_DF.columns if col not in ['sampleID', exp_factor]] 
        palette = [subcat_color_palette[cat] for cat in categories]      
        plt.figure(figsize=(10, 7.5))
        FSM_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Number of Reads in Each subcategory - FSM')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Reads')
        plt.xticks(rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        ##Plot 5
        FSM_perc_DF.set_index('sampleID', inplace=True)
        plt.figure(figsize=(10, 7.5))
        categories =  [col for col in FSM_perc_DF.columns if col not in ['sampleID', exp_factor]] 
        palette = [subcat_color_palette[cat] for cat in categories]      
        FSM_perc_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Number of Reads in Each subcategory - FSM')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage')
        plt.xticks(rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 6
        ISM_DF.set_index('sampleID', inplace=True)
        plt.figure(figsize=(10, 7.5))
        categories =  [col for col in ISM_DF.columns if col not in ['sampleID', exp_factor]] 
        palette = [subcat_color_palette[cat] for cat in categories]      
        ISM_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Number of Reads in Each subcategory - ISM')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Reads')
        plt.xticks(rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 7
        ISM_perc_DF.set_index('sampleID', inplace=True)
        plt.figure(figsize=(10, 7.5))
        categories =  [col for col in ISM_perc_DF.columns if col not in ['sampleID', exp_factor]] 
        palette = [subcat_color_palette[cat] for cat in categories]      
        ISM_perc_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Number of Reads in Each subcategory - ISM')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage')
        plt.xticks(rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
       
        #Plot 8
        NIC_NNC_DF.set_index('sampleID', inplace=True)
        plt.figure(figsize=(10, 7.5))
        categories =  [col for col in NIC_NNC_DF.columns if col not in ['sampleID', exp_factor]] 
        palette = [subcat_color_palette[cat] for cat in categories]      
        NIC_NNC_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Number of Reads in Each subcategory - NIC/NNC')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Reads')
        plt.xticks(rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 9
        NIC_NNC_perc_DF.set_index('sampleID', inplace=True)
        plt.figure(figsize=(10, 7.5))
        categories =  [col for col in NIC_NNC_perc_DF.columns if col not in ['sampleID', exp_factor]] 
        palette = [subcat_color_palette[cat] for cat in categories]      
        NIC_NNC_perc_DF.plot(kind='bar', stacked=True, color=palette)
        plt.title('Number of Reads in Each subcategory - NIC/NNC')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage')
        plt.xticks(rotation=90, ha='right')
        plt.legend(labels=categories, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 10 - Barplot genes detyected colored by number of reads in gene
        categories = ['1 read','2-10 reads', '11-50 reads', '50-100 reads','100+ reads']
        cols = ['sampleID'] + categories
        gene_agg_DF = gene_agg_DF[cols]
        gene_agg_DF.set_index('sampleID', inplace=True)
        palette = sns.color_palette("rainbow", len(categories))
     
        plt.figure(figsize=(10, 7.5))
        gene_agg_DF.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Genes detected')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Genes detected')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 11 - Barplot - % genes per sample coloured by number of reads in gene
        
        gene_percs_unstacked = gene_percs_unstacked[cols]
        gene_percs_unstacked.set_index('sampleID', inplace=True)
        plt.figure(figsize=(10, 7.5))
        gene_percs_unstacked.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Genes detected')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage')
        plt.xticks(rotation=90, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 12 -15 - Boxplots Distribution of % structural category (FSM ISM NIC NNC) 
        

        for category in melted_annotated_gene_DF['category'].unique():
             plt.figure(figsize=(10, 6))
             sns.boxplot(x='sampleID', y='percentage', hue='sampleID' ,data=melted_annotated_gene_DF[melted_annotated_gene_DF['category'] == category],
                    palette=sample_color_palette)
             plt.title(f'Gene distribution - {category}')
             plt.xticks(rotation=45)
             plt.tight_layout()
             pdf.savefig()
             plt.show()
             plt.close()


        
        #Plots 16 - 19 - UJC barplots
        for stack_by in ['read_category', 'structural_category']:
            # Create a figure and a set of subplots
            ujc_cnts_dct[stack_by] =  ujc_cnts_dct[stack_by].sort_values(by='sampleID')
            ujc_cnts_dct[stack_by] =  ujc_cnts_dct[stack_by].drop(columns=[exp_factor])
            ujc_percs_dct[stack_by] =  ujc_percs_dct[stack_by].sort_values(by = 'sampleID')
            ujc_percs_dct[stack_by] =  ujc_percs_dct[stack_by].drop(columns=[exp_factor])
            
            if stack_by == 'read_category':
                fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))  # Adjust figsize as needed
            
                # Plot number of UJCs per sample on the first subplot
                ujc_cnts_dct[stack_by].plot(kind='bar', stacked=True, colormap='rainbow', ax=axes[0])
                axes[0].set_title('Number of UJCs detected')
                axes[0].set_xticks(np.arange(len(ujc_cnts_dct[stack_by]['sampleID'])))
                axes[0].set_xticklabels(ujc_cnts_dct[stack_by]['sampleID'], rotation=90, ha='right')
                axes[0].set_xlabel('SampleID')
                axes[0].set_ylabel('Number of UJCs')
                axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
                # Plot % UJCs per sample on the second subplot
                ujc_percs_dct[stack_by].plot(kind='bar', stacked=True, colormap='rainbow', ax=axes[1])
                axes[1].set_title('Percentage of UJCs detected')
                axes[1].set_xticks(np.arange(len(ujc_percs_dct[stack_by]['sampleID'])))
                axes[1].set_xticklabels(ujc_percs_dct[stack_by]['sampleID'], rotation=90, ha='right')
                axes[1].set_xlabel('SampleID')
                axes[1].set_ylabel('Percentage of UJCs')
                axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
                plt.tight_layout()
                pdf.savefig()  # If saving to PDF, uncomment this line
                plt.show()
                plt.close(fig)  # Close the figure to free up memory
            elif stack_by == 'structural_category':
                
                categories = [col for col in ujc_cnts_dct[stack_by].columns if col not in ['sampleID', exp_factor]]
                colors = [category_color_palette[cat] for cat in categories]
                
                fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))  # Adjust figsize as needed
            
                # Plot number of UJCs per sample on the first subplot
                ujc_cnts_dct[stack_by].plot(kind='bar', stacked=True, color=colors, ax=axes[0])
                axes[0].set_title('Number of UJCs detected')
                axes[0].set_xticks(np.arange(len(ujc_cnts_dct[stack_by]['sampleID'])))
                axes[0].set_xticklabels(ujc_cnts_dct[stack_by]['sampleID'], rotation=90, ha='right')
                axes[0].set_xlabel('SampleID')
                axes[0].set_ylabel('Number of UJCs')
                axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
                # Plot % UJCs per sample on the second subplot
                ujc_percs_dct[stack_by].plot(kind='bar', stacked=True, color = colors, ax=axes[1])
                axes[1].set_title('Percentage of UJCs detected')
                axes[1].set_xticks(np.arange(len(ujc_percs_dct[stack_by]['sampleID'])))
                axes[1].set_xticklabels(ujc_percs_dct[stack_by]['sampleID'], rotation=90, ha='right')
                axes[1].set_xlabel('SampleID')
                axes[1].set_ylabel('Percentage of UJCs')
                axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
                plt.tight_layout()
                pdf.savefig() 
                plt.show()
                plt.close(fig)  
                
            
        ##Plot 20 - num reads vs % reads gt 1kb
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=length_DF, x='perc_reads_gt_1kb', y='total_reads', 
                        hue='sampleID', palette=sample_color_palette, legend='full', s=100)
        plt.title('Total Reads vs Percentage of Reads > 1kb')
        plt.xlabel('Percentage of Reads > 1kb')
        plt.ylabel('Total Reads')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        

        # Plot 21 : Barplot lengths of reads counts
        length_cnts_agg = length_cnts_agg.sort_values(by= 'sampleID')
        length_cnts_agg =  length_cnts_agg.drop(columns=[exp_factor])
        plt.figure(figsize=(10, 6))
        length_cnts_agg.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Number of reads')
        plt.xticks(ticks=np.arange(len(length_cnts_agg['sampleID'])), labels=length_cnts_agg['sampleID'], rotation=90, ha='right')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Reads')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Read Length Category')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()

        # Plot 22 - Barplot lengths of reads %
        length_percs_agg = length_percs_agg.sort_values(by= 'sampleID')
        length_percs_agg =  length_percs_agg.drop(columns=[exp_factor])
        plt.figure(figsize=(10, 6))
        length_percs_agg.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Number of reads')
        plt.xticks(ticks=np.arange(len(length_percs_agg['sampleID'])), labels=length_percs_agg['sampleID'], rotation=90, ha='right')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage of Reads')
        plt.legend(title='Read Length Category', loc='upper left', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 23 Violin plots
        sns.violinplot(x='sampleID', y='length', data=length_DF2, palette = sample_color_palette)
        plt.xlabel('Sample ID')
        plt.ylabel('Length')
        plt.title('Read Length Distribution')
        pdf.savefig()
        plt.show()
        plt.close()
        
        ##Plots 24 - 29 -  % structural category vs %reads greater than 1kb
        categories = [cat for cat in ['FSM', 'ISM', 'NIC', 'NNC','GI','GENIC'] if cat in length_DF.columns]
        palette = sample_color_palette
        for category in categories:
            plt.figure(figsize=(10, 6))
            sns.scatterplot(data=length_DF,y='perc_reads_gt_1kb', x=category, 
                            hue='sampleID', palette=palette, legend='full', s=100)
            plt.title('Percentage of Reads > 1kb vs %' + category )
            plt.ylabel('Percentage of Reads > 1kb')
            plt.xlabel('%' + category)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            pdf.savefig() 
            plt.show() 
            plt.close()
        
        ##Plot 30 - Artefact reads counts
        cnt_categories = ['num_reads_RTS','num_reads_intrapriming','num_reads_non-canonical']
        perc_categories = ['perc_reads_RTS','perc_reads_intrapriming','perc_reads_non-canonical' ]
        cnt_cols = ['sampleID'] + cnt_categories
        perc_cols =  ['sampleID'] + perc_categories
        err_cnt_DF = err_DF[cnt_cols]
        err_perc_DF = err_DF[perc_cols]
        err_cnt_DF.set_index('sampleID', inplace=True)
        err_perc_DF.set_index('sampleID', inplace=True)
        
        palette = sns.color_palette("rainbow", 3)
        
        
        plt.figure(figsize=(10, 6))
        err_cnt_DF.plot(kind='bar', stacked=False, colormap='rainbow')
        plt.title('Number of Artefact reads')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Reads')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        ##Plot 31 - Artefact reads %
        plt.figure(figsize=(10, 6))
        err_perc_DF.plot(kind='bar', stacked=False, colormap='rainbow')
        plt.title('Percent Artefact reads')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 32 - PCA plot
        palette = sample_color_palette
        plt.figure(figsize=(10, 7))
        sns.scatterplot(x=0, y=1, hue='sampleID', data=pca_DF, palette=sample_color_palette, s=50)
        plt.title('PCA Plot Based on sampleID')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()
        pdf.savefig()
        plt.show()
        plt.close()
        

        ##Plot 33 - Screeplot and Loadings heatmap
        cumulative_variance = np.cumsum(variance_ratio)
        n_components = np.argmax(cumulative_variance >= 0.85) + 1
        
        # Create the plots
        fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex='col', gridspec_kw={'width_ratios': [10, 3], 'height_ratios': [3, 10]})
        loadings_DF = loadings_DF.iloc[:, :n_components]
        link = linkage(loadings_DF, method='average')
        sorted_idx = leaves_list(link)
        loadings_DF = loadings_DF.iloc[sorted_idx]
        
        # Bar plot for explained variance (Scree Plot)
        x_tick_pos = [i + 0.5 for i in range(n_components)]
        ax[0, 0].bar(x_tick_pos, variance_ratio[:n_components], align='center', label='Individual explained variance')
        ax[0, 0].step(x_tick_pos, cumulative_variance[:n_components], where='mid', label='Cumulative explained variance')
        ax[0, 0].set_xticks(x_tick_pos)
        ax[0, 0].set_xticklabels([])  # Clear x tick labels here
        ax[0, 0].set_ylabel('Variance Explained')
        
        # Set x tick labels for the heatmap
        x_ticks = [f'PC{i+1}' for i in range(n_components)]
        sns.heatmap(loadings_DF, cmap="coolwarm", ax=ax[1, 0], cbar_ax=ax[1, 1], xticklabels=x_ticks)
        ax[1, 0].set_yticks(np.arange(loadings_DF.shape[0]) + 0.5)
        ax[1, 0].set_yticklabels(loadings_DF.index, rotation=0)
        ax[1, 0].set_xlabel('Principal Components')
        
        # Use the ax[0,1] for legend
        ax[0, 1].axis('off')  # Turn off the axis lines and labels
        handles, labels = ax[0, 0].get_legend_handles_labels()
        ax[0, 1].legend(handles, labels, loc='center')  # Place legend at the center of ax[0, 1]
        
        title = fig.suptitle('Screeplot and Heatmap of PC loadings',y=1.02, fontsize=16)
        plt.tight_layout()
        pdf.savefig(bbox_extra_artists=(title,), bbox_inches='tight')
        plt.show()
        plt.close()
        
        #Plot 34 - Barplot junction category counts
        nov_can_DF = nov_can_DF.sort_values(by= 'sampleID')
        nov_can_DF  =  nov_can_DF.drop(columns=[exp_factor])
        plt.figure(figsize=(10, 6))
        nov_can_DF.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Junctions by Category')
        plt.xticks(ticks=np.arange(len(nov_can_DF['sampleID'])), labels=nov_can_DF['sampleID'], rotation=90, ha='right')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Junctions')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 35 - Barplot junction category percentages
        nov_can_perc_DF = nov_can_perc_DF.sort_values(by= 'sampleID')
        nov_can_perc_DF  =  nov_can_perc_DF.drop(columns=[exp_factor])
        plt.figure(figsize=(10, 6))
        nov_can_perc_DF.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Junctions by Category')
        plt.xticks(ticks=np.arange(len(nov_can_perc_DF['sampleID'])), labels=nov_can_perc_DF['sampleID'], rotation=90, ha='right')
        plt.xlabel('SampleID')
        plt.ylabel('Percentage')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 36 - Barplot cv junction categories acceptors
        cv_acc_summary = cv_acc_summary.sort_values(by= 'sampleID')
        cv_acc_summary  =  cv_acc_summary.drop(columns=[exp_factor])
        plt.figure(figsize=(10, 6))
        cv_acc_summary.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Number of Acceptors > 3 reads')
        plt.xticks(ticks=np.arange(len(cv_acc_summary['sampleID'])), labels=cv_acc_summary['sampleID'], rotation=90, ha='right')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Acceptors')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
        
        #Plot 37 - Barplot cv junction categories donors
        cv_don_summary = cv_don_summary.sort_values(by= 'sampleID')
        cv_don_summary  =  cv_don_summary.drop(columns=[exp_factor])
        plt.figure(figsize=(10, 6))
        cv_don_summary.plot(kind='bar', stacked=True, colormap='rainbow')
        plt.title('Number of Donors > 3 reads')
        plt.xticks(ticks=np.arange(len(cv_don_summary['sampleID'])), labels=cv_don_summary['sampleID'], rotation=90, ha='right')
        plt.xlabel('SampleID')
        plt.ylabel('Number of Donors')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()
        plt.close()
    
    
def main():
    ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs, fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct = proc_samples(args.inDESIGN, args.inREF)
    
    gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF,FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct = prep_tables(ref_DF, gene_count_dfs,ujc_count_dfs,length_dfs,cv_dfs, err_dfs,
                                                                                                                            fsm_dfs, ism_dfs, nic_nnc_dfs, nov_can_dfs,length_Dct)
    identify_cand_underannot(gene_count_DF)
    
    dfs_for_plotting = prep_data_4_plots( gene_count_DF, ujc_count_DF, length_DF, cv_DF, err_DF, FSM_DF, ISM_DF, NIC_NNC_DF, nov_can_DF, length_Dct )
    
    if args.GENERATE_PLOTS:
        if args.inFACTOR == None:
            
            plot_pdf(os.path.join(args.OUT, args.PREFIX + '_plots.pdf'), *dfs_for_plotting)
        else: 
            plot_pdf_by_factor(os.path.join(args.OUT, args.PREFIX + '_plots.pdf'), *dfs_for_plotting)
        
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

