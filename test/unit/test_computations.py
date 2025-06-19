import pytest
import os, sys
import pandas as pd

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append(main_path)
data_path = os.path.join(main_path, 'test','test_data')
from src.qc_computations import (
    process_rts
)

@pytest.fixture
def results_dataframe(data_path):
    pd.read_csv(os.path.join(data_path,"test_isoforms_classification.csv"))

@pytest.fixture
def genome_dict(data_path):
    return dict((r.name, r) for r in SeqIO.parse(open(os.path.join(data_path,"genome_test.fasta")), 'fasta'))

@pytest.fixture
def reference_data(data_path,genome_dict):
    annotation = os.path.join(data_path,"test_isoforms.gtf")
    ref = reference_parser(annotation,data_path,"test",list(genome_dict.keys()))
    
    isoform = os.path.join(data_path,"test_isoforms.genePred")
    isoforms_dict = isoforms_parser(isoform)

    for _, records in isoform_data.items():
    
        for record in records:
            result, novel_gene_index = classify_isoform(record, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr,
                                      junctions_by_gene, start_ends_by_gene, genome_dict,novel_gene_index)
            
            isoforms_info = isoformClassification(isoforms_dict, ref)   
