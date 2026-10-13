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
from src.qc_output import write_junction_output

# Test rts function

@pytest.fixture
def junctions_file():
    return os.path.join(main_path, "test/test_data/junctions/test_junctions.txt")

@pytest.fixture
def genome():
    return os.path.join(main_path, "test/test_data/genome/genome_test.fasta")

@pytest.fixture
def genome_dict(genome):
    return dict((r.name, r) for r in SeqIO.parse(open(genome), 'fasta'))

def test_result(junctions_file,genome,genome_dict):
    RTS_info = rts([junctions_file, genome, "-a"], genome_dict)
    assert RTS_info['PB.137296.2'] == ['junction_1','junction_3']
    assert RTS_info['PB.137356.1'] == ['junction_1']
    # test that should raise and error
    with pytest.raises(KeyError):
        assert RTS_info['PB.103688.1']

def test_result_shared(junctions_file,genome,genome_dict):
    RTS_info = rts([junctions_file, genome, "-a"], genome_dict)
    assert RTS_info['PB.124736.2'] == ['junction_3']
    assert RTS_info['PB.124736.3'] == ['junction_5']

## Test isoforms_info
@pytest.fixture
def isoforms_info(genome_dict):
    isoforms_dict = isoforms_parser(os.path.join(main_path, "test/test_data/isoforms/test_isoforms.genePred"))

    annotation = os.path.join(main_path, "test/test_data/reference/test_reference.gtf")
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
    return pd.read_csv(os.path.join(main_path,"test/test_data/isoforms/test_isoforms_classification.tsv"),sep="\t")
    

def test_isoforms_output(isoforms_info,junctions_file,genome,genome_dict,result_dataframe):
    isoforms_res, _ = process_rts(isoforms_info, junctions_file,
                                                   genome,genome_dict,"")
    for _, result in isoforms_res.items():
        expected_row = result_dataframe.loc[result_dataframe["isoform"] == result.isoform]
        assert result.RTS_stage.upper() == str(expected_row["RTS_stage"].values[0]).upper(), f"Failed for {result.isoform}. Expected {expected_row['RTS_stage'].values[0]}, got {result.RTS_stage}"

    
# Junction output file

def test_junction_output_rts_column(junctions_file, genome, genome_dict, tmp_path):
    """Every junction gets TRUE/FALSE in RTS_junction, including isoforms with no RT-switching junction."""
    junc = pd.read_csv(junctions_file, sep="\t")
    junc["RTS_junction"] = "????"  # placeholder written to the _tmp file by the QC
    out_path = str(tmp_path / "test_junctions.txt")
    junc.to_csv(out_path + "_tmp", sep="\t", index=False)

    _, RTS_info = process_rts({}, out_path, genome, genome_dict)
    write_junction_output(out_path, RTS_info, list(junc.columns))

    result = pd.read_csv(out_path, sep="\t", dtype={"RTS_junction": str})
    assert len(result) == len(junc)
    assert set(result["RTS_junction"]) == {"TRUE", "FALSE"}
    # isoforms without any RT-switching junction must be FALSE, not left as ????
    assert "PB.103724.1" not in RTS_info
    assert (result.loc[result["isoform"] == "PB.103724.1", "RTS_junction"] == "FALSE").all()
    # isoforms with RT-switching junctions: only those junctions are TRUE
    iso = result[result["isoform"] == "PB.137296.2"].set_index("junction_number")["RTS_junction"]
    assert iso["junction_1"] == "TRUE" and iso["junction_3"] == "TRUE"
    assert (iso.drop(["junction_1", "junction_3"]) == "FALSE").all()
    rts_true = {(i, j) for i, js in RTS_info.items() for j in js}
    assert set(map(tuple, result.loc[result["RTS_junction"] == "TRUE", ["isoform", "junction_number"]].values)) == rts_true
