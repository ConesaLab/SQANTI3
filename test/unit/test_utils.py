import pytest,sys,os, math
from collections.abc import Iterable

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
# Assuming the functions are in a file named 'utils.py'
from src.utils import mergeDict, flatten, pstdev, find_polyA_motif

# Tests for mergeDict function
def test_merge_non_overlapping_dicts():
    dict1 = {'a': 1, 'b': 2}
    dict2 = {'c': 3, 'd': 4}
    result = mergeDict(dict1, dict2)
    assert result == {'a': 1, 'b': 2, 'c': 3, 'd': 4}

def test_merge_overlapping_dicts():
    dict1 = {'a': 1, 'b': 2, 'c': 3}
    dict2 = {'b': 4, 'c': 5, 'd': 6}
    result = mergeDict(dict1, dict2)
    assert result == {'a': 1, 'b': [4, 2], 'c': [5, 3], 'd': 6}

def test_merge_empty_dicts():
    dict1 = {}
    dict2 = {}
    result = mergeDict(dict1, dict2)
    assert result == {}

def test_merge_one_empty_dict():
    dict1 = {'a': 1, 'b': 2}
    dict2 = {}
    result = mergeDict(dict1, dict2)
    assert result == {'a': 1, 'b': 2}

# Tests for flatten function
def test_flatten_list():
    nested_list = [1, [2, 3, [4, 5]], 6, [7, 8]]
    result = list(flatten(nested_list))
    assert result == [1, 2, 3, 4, 5, 6, 7, 8]

def test_flatten_with_strings():
    nested_list = [1, [2, "hello", [4, 5]], 6, ["world"]]
    result = list(flatten(nested_list))
    assert result == [1, 2, "hello", 4, 5, 6, "world"]

def test_flatten_empty_list():
    empty_list = []
    result = list(flatten(empty_list))
    assert result == []

def test_flatten_non_nested_list():
    non_nested_list = [1, 2, 3, 4, 5]
    result = list(flatten(non_nested_list))
    assert result == [1, 2, 3, 4, 5]

# Tests for pstdev function
def test_pstdev_basic():
    data = [1, 2, 3, 4, 5]
    result = pstdev(data)
    expected = math.sqrt(2)  # The population standard deviation of [1,2,3,4,5] is sqrt(2)
    assert math.isclose(result, expected, rel_tol=1e-9)

def test_pstdev_constant():
    data = [7, 7, 7, 7, 7]
    result = pstdev(data)
    assert result == 0

def test_pstdev_single_element():
    data = [42]
    result = pstdev(data)
    assert result == 0

def test_pstdev_empty_list():
    with pytest.raises(ZeroDivisionError):
        pstdev([])

# Tests for find_polyA_motif function
def test_find_polyA_motif_found():
    genome_seq = "ACGTAATAAA"
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("AATAAA", -1, "TRUE")

def test_find_polyA_motif_not_found():
    genome_seq = "ACGTCGTACG"
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("NA", "NA", "FALSE")

def test_find_polyA_motif_multiple_motifs():
    genome_seq = "ACGTATTAAATAAA"
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("ATTAAA", -6, "TRUE")

def test_find_polyA_motif_at_start():
    genome_seq = "AATAAACGT"
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("AATAAA", -3, "TRUE")

def test_find_polyA_motif_at_end():
    genome_seq = "ACGTAATAAA"
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("AATAAA", -1, "TRUE")

def test_find_polyA_motif_empty_sequence():
    genome_seq = ""
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("NA", "NA", "FALSE")

def test_find_polyA_motif_empty_motif_list():
    genome_seq = "ACGTAATAAA"
    polyA_motif_list = []
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("NA", "NA", "FALSE")
