

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

    
    
    #reassign counts to surviving isoforms
    rescued = rescued[rescued['rescue_candidate'].isin(old_counts.keys())]
    inclusion_df = pd.concat([inclusion_list, not_rescued['transcript_id']])
    inclusion_df.apply(lambda x: select_hit(x.iloc[0], rescued, new_counts, old_counts), axis = 1)



    # Save transcript divergency relations
    not_rescued.to_csv(f"{prefix}_transcript_divergency_relations.tsv", header = True, index = False, sep = '\t')
