import pandas as pd
import pytest,sys,os
from Bio import SeqIO

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.classification_steps import classify_isoform
from src.parsers import isoforms_parser, reference_parser
from src.qc_computations import (
    process_rts
)
from src.utilities.rt_switching import rts

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

def test_result(junction_file,genome,genome_dict):
    RTS_info = rts([junction_file, genome, "-a"], genome_dict)
    assert RTS_info['PB.137296.2'] == ['junction_1','junction_3']
    assert RTS_info['PB.137356.1'] == ['junction_1']
    # test that should raise and error
    with pytest.raises(KeyError):
        assert RTS_info['PB.103688.1']

def test_result_shared(junction_file,genome,genome_dict):
    RTS_info = rts([junction_file, genome, "-a"], genome_dict)
    assert RTS_info['PB.124736.2'] == ['junction_3']
    assert RTS_info['PB.124736.3'] == ['junction_5']

## Test isoforms_info
@pytest.fixture
def isoforms_info(genome_dict):
    isoforms_dict = isoforms_parser(os.path.join(main_path, "test/test_data/test_isoforms.genePred"))

    annotation = os.path.join(main_path, "test/test_data/test_reference.gtf")
    refs_1exon_by_chr, refs_exons_by_chr, \
    junctions_by_chr, junctions_by_gene, start_ends_by_gene = reference_parser(annotation,os.path.join(main_path,"test/test_data"),
                                                                               "test",list(genome_dict.keys()))
    
    isoforms_res = {}
    for _, records in isoforms_dict.items():
        for record in records:
            isoforms_res[record.id]  = classify_isoform(record, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr,
                                      junctions_by_gene, start_ends_by_gene, genome_dict)

    return isoforms_res

@pytest.fixture
def result_dataframe():
    return pd.read_csv(os.path.join(main_path,"test/test_data/test_isoforms_classification.tsv"),sep="\t")
    

def test_isoforms_output(isoforms_info,junction_file,genome,genome_dict,result_dataframe):
    isoforms_res, _ = process_rts(isoforms_info, junction_file,
                                                   genome,genome_dict,"")
    for _, result in isoforms_res.items():
        expected_row = result_dataframe.loc[result_dataframe["isoform"] == result.id]
        print(result.id,result.RT_switching)
        assert result.RT_switching.upper() == str(expected_row["RTS_stage"].values[0]).upper(), f"Failed for {result.id}. Expected {expected_row['RTS_stage'].values[0]}, got {result.RT_switching}"

    
# Junction output file
