import pytest,sys,os

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)
from src.utilities.rt_switching import (
    SpliceJunctions, loadSpliceJunctions,
    checkForRepeatPat
)

@pytest.fixture
def junctions_file():
    return os.path.join(main_path, "test/test_data/junctions_test.txt")


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
    assert isinstance(sj_seen_counts[sample_junction], int)

def test_no_duplicate_junctions(loaded_data):
    sj_dict, sj_seen_counts = loaded_data
    unique_junctions = sum(len(junctions) for junctions in sj_dict.values())
    total_junctions = sum(sj_seen_counts.values())
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
    assert sj_dict_count[ ('chr22', '+', '15815567', '15817907')] == 2

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