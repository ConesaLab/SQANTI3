import os, re

__author__  = "etseng@pacb.com"
__version__ = '5.5.1'  # Python 3.7
utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
default_json = os.path.abspath(utilitiesPath + "/filter/filter_default.json")


FIELDS_JUNC = ['isoform', 'chrom', 'strand', 'junction_number', 'genomic_start_coord',
                   'genomic_end_coord', 'transcript_coord', 'junction_category',
                   'start_site_category', 'end_site_category', 'diff_to_Ref_start_site',
                   'diff_to_Ref_end_site', 'bite_junction', 'splice_site', 'canonical',
                   'RTS_junction', 'indel_near_junct',
                   'phyloP_start', 'phyloP_end', 'sample_with_cov', "total_coverage_unique", "total_coverage_multi"] #+coverage_header

FIELDS_CLASS = ['isoform', 'chrom', 'strand', 'length',  'exons',  'structural_category',
                'associated_gene', 'associated_transcript',  'ref_length', 'ref_exons',
                'diff_to_TSS', 'diff_to_TTS', 'diff_to_gene_TSS', 'diff_to_gene_TTS',
                'subcategory', 'RTS_stage', 'all_canonical',
                'min_sample_cov', 'min_cov', 'min_cov_pos',  'sd_cov', 'FL', 'n_indels',
                'n_indels_junc',  'bite',  'iso_exp', 'gene_exp',  'ratio_exp',
                'FSM_class',   'coding', 'ORF_length', 'CDS_length', 'CDS_start',
                'CDS_end', 'CDS_genomic_start', 'CDS_genomic_end', 'predicted_NMD',
                'perc_A_downstream_TTS', 'seq_A_downstream_TTS',
                'dist_to_CAGE_peak', 'within_CAGE_peak',
                'dist_to_polyA_site', 'within_polyA_site',
                'polyA_motif', 'polyA_dist', 'polyA_motif_found', 'ORF_seq', 'ratio_TSS']

# Sequence names
seqid_rex1 = re.compile(r'PB\.(\d+)\.(\d+)$')
seqid_rex2 = re.compile(r'PB\.(\d+)\.(\d+)\|\S+')
seqid_fusion = re.compile(r'PBfusion\.(\d+)\.(\d+)\S*')

EXP_KALLISTO_HEADERS = ['target_id', 'length', 'eff_length', 'est_counts', 'tpm']
EXP_RSEM_HEADERS = ['transcript_id', 'length', 'effective_length', 'expected_count', 'TPM']
