from csv import DictWriter
import pandas as pd
import pytest,sys,os
from Bio import SeqIO

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.classification_steps import classify_isoform
from src.parsers import isoforms_parser, reference_parser
from src.classification_steps import write_junction_info
from src.config import  FIELDS_JUNC

# Test rts function

@pytest.fixture
def junction_file():
    return os.path.join(main_path, "test/test_data/test_junctions.txt")

@pytest.fixture
def genome():
    return os.path.join(main_path, "test/test_data/genome_test.fasta")

@pytest.fixture
def genome_dict(genome):
    return dict((r.name, r) for r in SeqIO.parse(open(genome), 'fasta'))

@pytest.fixture
def junctions_dicts(genome_dict):
    annotation = os.path.join(main_path, "test/test_data/test_reference.gtf")
    refs_1exon_by_chr, refs_exons_by_chr, \
        junctions_by_chr, junctions_by_gene, \
            start_ends_by_gene = reference_parser(annotation,os.path.join(main_path,"test/test_data"),
                                                                               "test",list(genome_dict.keys()))
    
    return refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene

## Test isoforms_info
@pytest.fixture
def isoforms_info(genome_dict,junctions_dicts):
    isoforms_res = {}
    refs_1exon_by_chr, refs_exons_by_chr, \
        junctions_by_chr, junctions_by_gene, \
            start_ends_by_gene = junctions_dicts

    for _, records in isoforms_dict.items():
        for record in records:
            print(record.id,record.junctions)
            isoforms_res[record.id]  = classify_isoform(record, refs_1exon_by_chr, refs_exons_by_chr, \
                                                        junctions_by_chr, junctions_by_gene, \
                                                            start_ends_by_gene , genome_dict)
            print(isoforms_res[record.id].junctions)
    return isoforms_res

@pytest.fixture
def result_dataframe():
    return pd.read_csv(os.path.join(main_path,"test/test_data/test_junctions.txt"),sep="\t")

def test_initial_junctions(junctions_dicts,genome_dict,result_dataframe):
    outputJuncPath = os.path.join(main_path,"test","test_data","intermediates","junctions_temporal.txt")
    try:
        os.remove(outputJuncPath)
    except FileNotFoundError:
        pass
    # Defining function imputs
    _, _, junctions_by_chr, _, _ = junctions_dicts
    sites = "ATAC,GCAG,GTAG"
    accepted_cannonical_sites = list(sites.split(","))
    indelsJunc = None
    # Creating file header
    fields_junc_cur = FIELDS_JUNC
    
    handle_junc = open(outputJuncPath, "w")
    fout_junc = DictWriter(handle_junc, fieldnames=fields_junc_cur, delimiter='\t')
    fout_junc.writeheader()
    isoforms_dict = isoforms_parser(os.path.join(main_path, "test/test_data/test_isoforms.genePred"))

    for rec in isoforms_dict["chr22"]:
        write_junction_info(rec,junctions_by_chr,accepted_cannonical_sites,indelsJunc,genome_dict,fout_junc)
    handle_junc.close()
    junctions_res = pd.read_csv(outputJuncPath,sep="\t")
    
    #Order junctions
    # Add a temporary order column to result_dataframe
    result_dataframe['temp_order'] = range(len(result_dataframe))

    # Merge junctions_res with result_dataframe
    merged_df = pd.merge(result_dataframe[['isoform','junction_number', 'temp_order']], 
                        junctions_res, 
                        on=['isoform','junction_number'], 
                        how='left')

    # Sort by the temporary order
    merged_df = merged_df.sort_values('temp_order')

    # Drop the temporary order column
    junctions_res = merged_df.drop('temp_order', axis=1)

    # Reset the index if needed
    junctions_res = junctions_res.reset_index(drop=True)

    columns_to_compare = [
        'junction_category', 'start_site_category', 'end_site_category',
        'diff_to_Ref_start_site', 'diff_to_Ref_end_site', 'bite_junction',
        'splice_site', 'canonical'
    ]
    
    for index, row in junctions_res.iterrows():
        for column in columns_to_compare:
            try:
                assert row[column] == result_dataframe.loc[index, column], \
                    f"Mismatch in row {index}, column {column}: {row[column]} != {result_dataframe.loc[index, column]}"
            except AssertionError as e:
                isoform = row['isoform']
                junction_number = row['junction_number']
                error_message = f"Error in isoform: {isoform}, junction_number: {junction_number}\n{str(e)}"
                pytest.fail(error_message)
    