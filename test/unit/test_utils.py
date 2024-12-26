import pytest,sys,os, math
from collections.abc import Iterable

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
# Assuming the functions are in a file named 'utils.py'
from src.utils import get_files_from_dir, mergeDict, flatten, pstdev, find_polyA_motif

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
    polyA_motif_list = ["ACTAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("ATTAAA", -5, "TRUE")

def test_find_polyA_motif_at_start():
    genome_seq = "AATAAACGT"
    polyA_motif_list = ["AATAAA", "ATTAAA"]
    result = find_polyA_motif(genome_seq, polyA_motif_list)
    assert result == ("AATAAA", -4, "TRUE")

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
    
@pytest.fixture
def setup_test_files(tmp_path):
    # Create a temporary directory and files for testing
    test_dir = tmp_path / "test_dir"
    test_dir.mkdir()

    # Create some test files with different extensions
    (test_dir / "file1.txt").write_text("content1")
    (test_dir / "file2.txt").write_text("content2")
    (test_dir / "file3.csv").write_text("content3")

    # Create a file containing a list of files
    file_list = tmp_path / "file_list.txt"
    file_list.write_text(str(test_dir / "file1.txt") + "\n" + str(test_dir / "file2.txt") + "\n")

    return test_dir, file_list

def test_get_files_from_dir_with_directory(setup_test_files):
    test_dir, _ = setup_test_files
    result = get_files_from_dir(test_dir, ".txt")
    expected = [str(test_dir / "file1.txt"), str(test_dir / "file2.txt")]
    assert sorted(result) == sorted(expected)

def test_get_files_from_dir_with_file(setup_test_files):
    _, file_list = setup_test_files
    result = get_files_from_dir(file_list, ".txt")
    expected = [str(file_list.parent / "test_dir" / "file1.txt"), str(file_list.parent / "test_dir" / "file2.txt")]
    assert sorted(result) == sorted(expected)