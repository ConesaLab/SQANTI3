import os, sys, pytest
from Bio import SeqIO
import pandas as pd
main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.classification_main import (
    classify_isoform
)

from src.parsers import (
    reference_parser, isoforms_parser
)
novel_gene_index = 1

@pytest.fixture
def data_path():
    return os.path.join(main_path,"test","test_data")

@pytest.fixture
def result_dataframe(data_path):
    return pd.read_csv(os.path.join(data_path,"test_isoforms_classification.csv"))

@pytest.fixture
def genome_dict(data_path):
    return dict((r.name, r) for r in SeqIO.parse(open(os.path.join(data_path,"genome_test.fasta")), 'fasta'))


@pytest.fixture
def reference_data(data_path,genome_dict):
    annotation = os.path.join(data_path,"test_isoforms.gtf")
    
    return reference_parser(annotation,data_path,"test",list(genome_dict.keys()))
     
@pytest.fixture
def isoform_data(data_path):
    isoform = os.path.join(data_path,"test_isoforms.genePred")
    isoforms_dict = isoforms_parser(isoform)
    return isoforms_dict

def test_isoformClassification(reference_data, result_dataframe, isoform_data, genome_dict):
    novel_gene_index = 1
    refs_1exon_by_chr, refs_exons_by_chr, \
    junctions_by_chr, junctions_by_gene, start_ends_by_gene = reference_data

    expected_results = result_dataframe

    for _, records in isoform_data.items():
    
        for record in records:
            result, novel_gene_index = classify_isoform(record, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr,
                                      junctions_by_gene, start_ends_by_gene, genome_dict,novel_gene_index)
            
            expected_row = expected_results.loc[expected_results["isoform"] == result.id]
            expected_str_class = expected_row["structural_category"].values[0]
            expected_subtype = expected_row["subcategory"].values[0]
            expected_genes = expected_row["associated_gene"].values[0]
            assert result.str_class == expected_str_class, (
                            f"Isoform {result.id}. Expected {expected_str_class};{expected_subtype}, got {result.str_class};{result.subtype}"
            )
            assert result.subtype == expected_subtype, (
                f"Isoform {result.id}. Expected {expected_str_class};{expected_subtype}, got {result.str_class};{result.subtype}"
            )
            detected_genes = "_".join(sorted(set(result.genes)))
            assert detected_genes == expected_genes, (
                f"Isoform {result.id}. Expected {expected_genes}, got {result.genes}"
            )

