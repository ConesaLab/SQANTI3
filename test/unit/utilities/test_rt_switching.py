import pytest,sys,os
import csv
from Bio import SeqIO

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)
from src.utilities.rt_switching import (
    SpliceJunctions, loadSpliceJunctions,
    checkForRepeatPat, checkSJforRTS, FIELDS_RTS
)

@pytest.fixture
def junctions_file():
    return os.path.join(main_path, "test/test_data/junctions/junctions_test.txt")


@pytest.fixture
def loaded_data(junctions_file):
    return loadSpliceJunctions(junctions_file)

def test_loadSpliceJunctions_return_types(loaded_data):
    sj_dict, sj_seen_counts = loaded_data
    assert isinstance(sj_dict, dict)
    assert isinstance(sj_seen_counts, dict)

def test_sj_dict_structure(loaded_data):
    sj_dict, _ = loaded_data
    assert len(sj_dict) > 0  # Replace with the actual number of transcripts
    sample_transcript = list(sj_dict.keys())[0]
    assert isinstance(sj_dict[sample_transcript], list)
    assert isinstance(sj_dict[sample_transcript][0], SpliceJunctions)

def test_splice_junctions_object(loaded_data):
    sj_dict, _ = loaded_data
    sample_transcript = list(sj_dict.keys())[0]
    sample_sj = sj_dict[sample_transcript][0]
    assert sample_sj.trans == sample_transcript
    assert sample_sj.chromo == "chr22"  # Replace with actual chromosome
    assert sample_sj.strand in ("+", "-")
    assert isinstance(sample_sj.sjn, str)
    assert isinstance(sample_sj.strpos, int)
    assert isinstance(sample_sj.endpos, int)
    assert sample_sj.transpos is None
    assert sample_sj.category in ("novel", "known")
    assert sample_sj.startCat in ("known", "novel")  # Replace with actual categories
    assert sample_sj.endCat in ("known", "novel")  # Replace with actual categories
    assert sample_sj.type in ("canonical", "non_canonical")

def test_sj_seen_counts_structure(loaded_data):
    _, sj_seen_counts = loaded_data
    assert len(sj_seen_counts) > 0  # Replace with the actual number of unique junctions
    sample_junction = list(sj_seen_counts.keys())[0]
    assert isinstance(sample_junction, tuple)
    assert len(sample_junction) == 4  # (chrom, strand, start, end)
    assert isinstance(sj_seen_counts[sample_junction][0], str)

def test_no_duplicate_junctions(loaded_data):
    sj_dict, sj_seen_counts = loaded_data
    unique_junctions = sum(len(junctions) for junctions in sj_dict.values())
    total_junctions = 0
    for val in sj_seen_counts.values():
        total_junctions += len(val)
    assert total_junctions >= unique_junctions

def test_all_junctions_in_seen_counts(loaded_data):
    sj_dict, sj_seen_counts = loaded_data
    for transcript_junctions in sj_dict.values():
        for sj in transcript_junctions:
            junction_key = (sj.chromo, sj.strand, str(sj.strpos), str(sj.endpos))
            assert junction_key in sj_seen_counts

def test_valid_junction_categories(loaded_data):
    sj_dict, _ = loaded_data
    for transcript_junctions in sj_dict.values():
        for sj in transcript_junctions:
            assert sj.type in ("canonical", "non_canonical")
            assert sj.category in ("novel", "known")

def test_valid_junction(loaded_data):
    sj_dict, sj_dict_count = loaded_data
    assert sj_dict["PB.103684.1"] == []
    assert len(sj_dict["PB.103684.2"]) == 1
    assert len(sj_dict_count[('chr22', '+', '15815567', '15817907')]) == 2
    assert sj_dict_count[('chr22', '+', '15815567', '15817907')] == ['PB.8631.2-junction_5', 'PB.3796.2-junction_7']

### checkForRepeatPat ###

def test_exact_match():
    seq_exon = "ATCGATCG"
    seq_intron = "GGGATCGATCGGG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == True
    assert matchLen == 8
    assert matchPattern == "ATCGATCG"
    assert mismatch == 0

def test_no_match():
    seq_exon = "ATCGATCG"
    seq_intron = "GGGTTTAAAGG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == False
    assert matchLen == None
    assert matchPattern == None
    assert mismatch == None

def test_match_with_mismatch():
    seq_exon = "ATCGAACTG"
    seq_intron = "GGGATCGTACTGGG"
    min_match = 9
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == True
    assert matchLen == 9
    assert matchPattern == "ATCGAACTG" # TODO: If the match length is set to be bigger, there will be some issues with getting the right match sequence. It is because it fails to move a bit backwards
    assert mismatch == 1

def test_match_without_mismatch_allowed():
    seq_exon = "ATCGATCG"
    seq_intron = "GGGATCGXTCGGG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match, allow_mismatch=False)
    
    assert is_RTS == False
    assert matchLen == None
    assert matchPattern == None
    assert mismatch == None

def test_partial_match():
    seq_exon = "ATCGATCG"
    seq_intron = "GGGATCGGG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == False
    assert matchLen == None
    assert matchPattern == None
    assert mismatch == None

def test_match_at_start():
    seq_exon = "ATCGATCG"
    seq_intron = "ATCGATCGGG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == True
    assert matchLen == 8
    assert matchPattern == "ATCGATCG"
    assert mismatch == 0

def test_match_at_end():
    seq_exon = "ATCGATCG"
    seq_intron = "GGGATCGATCG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == True
    assert matchLen == 8
    assert matchPattern == "ATCGATCG"
    assert mismatch == 0

def test_long_sequences():
    seq_exon = "A" * 1000 + "ATCGATCG" + "A" * 1000
    seq_intron = "G" * 1000 + "ATCGATCG" + "G" * 1000
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == True
    assert matchLen == 8
    assert matchPattern == "ATCGATCG"
    assert mismatch == 0

def test_min_match_larger_than_sequence():
    seq_exon = "ATCG"
    seq_intron = "ATCG"
    min_match = 8
    
    is_RTS, matchLen, matchPattern, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match)
    
    assert is_RTS == False
    assert matchLen == None
    assert matchPattern == None
    assert mismatch == None


### checkSJfor RTS ###

@pytest.fixture
def mock_output_file(tmp_path):
    return tmp_path / "output.tsv"

@pytest.fixture
def mock_sj_dict(loaded_data):
    return loaded_data[0]

@pytest.fixture
def mock_sj_relations_dict(loaded_data):
    return loaded_data[1]

@pytest.fixture
def mock_genome_dict():
    genome_path = os.path.join(main_path, "test/test_data/genome/genome_test.fasta")
    return dict((r.id, r) for r in SeqIO.parse(open(genome_path), 'fasta'))

# Tests
def test_checkSJforRTS_basic(mock_output_file,mock_sj_relations_dict,
                             mock_sj_dict,mock_genome_dict):
    result = checkSJforRTS(
        mock_sj_dict,
        mock_sj_relations_dict,
        mock_genome_dict,
        wiggle_count=5,
        include_category='a',
        include_type='a',
        min_match=6,
        allow_mismatch=True,
        output_filename=str(mock_output_file)
    )
    
    assert isinstance(result, dict)
    assert 'PB.103684.2' in result
    assert isinstance(result['PB.103684.2'], list)

def test_checkSJforRTS_output_file(mock_output_file,mock_sj_relations_dict,
                                   mock_sj_dict,mock_genome_dict):
    checkSJforRTS(
        mock_sj_dict,
        mock_sj_relations_dict,
        mock_genome_dict,
        wiggle_count=5,
        include_category='a',
        include_type='a',
        min_match=6,
        allow_mismatch=True,
        output_filename=str(mock_output_file)
    )
    
    assert mock_output_file.exists()
    with open(mock_output_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames
        assert headers == FIELDS_RTS  # Assuming FIELDS_RTS is defined in your module

def test_checkSJforRTS_filtering(mock_output_file,mock_sj_relations_dict,
                                 mock_sj_dict,mock_genome_dict):
    result = checkSJforRTS(
        mock_sj_dict,
        mock_sj_relations_dict,
        mock_genome_dict,
        wiggle_count=5,
        include_category='k',  # Only novel
        include_type='n',  # Only canonical
        min_match=6,
        allow_mismatch=True,
        output_filename=str(mock_output_file)
    )
    
    with pytest.raises(KeyError):
        assert result['PB.103684.2']

    assert len(result['PB.15672.2']) > 0

def test_checkSJforRTS_strand_handling(mock_output_file,mock_sj_relations_dict,
                                       mock_sj_dict,mock_genome_dict):
    result = checkSJforRTS(
        mock_sj_dict,
        mock_sj_relations_dict,
        mock_genome_dict,
        wiggle_count=5,
        include_category='a',
        include_type='a',
        min_match=6,
        allow_mismatch=True,
        output_filename=str(mock_output_file)
    )
    
    # Check if both strands are processed
    assert len(result['PB.83093.1']) > 0
    assert len(result['PB.137289.1']) > 0

    # Verify strand-specific processing in output file
    with open(mock_output_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        strands = set(row['strand'] for row in reader)
        assert strands == {'+', '-'}