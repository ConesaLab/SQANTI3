import re, sys, os
# If the path to where the main sqanti3 directory is not in the system path, our modules wont be loaded
main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.config import (
    FIELDS_JUNC, FIELDS_CLASS, seqid_rex1, seqid_rex2, seqid_fusion,
    EXP_KALLISTO_HEADERS, EXP_RSEM_HEADERS
)

def test_FIELDS_JUNC_structure():
    assert isinstance(FIELDS_JUNC, list), "FIELDS_JUNC should be a list"
    assert all(isinstance(field, str) for field in FIELDS_JUNC), "All elements in FIELDS_JUNC should be strings"

def test_FIELDS_CLASS_structure():
    assert isinstance(FIELDS_CLASS, list), "FIELDS_CLASS should be a list"
    assert all(isinstance(field, str) for field in FIELDS_CLASS), "All elements in FIELDS_CLASS should be strings"

def test_regex_patterns():
    assert isinstance(seqid_rex1, re.Pattern), "seqid_rex1 should be a compiled regex pattern"
    assert isinstance(seqid_rex2, re.Pattern), "seqid_rex2 should be a compiled regex pattern"
    assert isinstance(seqid_fusion, re.Pattern), "seqid_fusion should be a compiled regex pattern"

def test_EXP_KALLISTO_HEADERS_structure():
    assert isinstance(EXP_KALLISTO_HEADERS, list), "EXP_KALLISTO_HEADERS should be a list"
    assert EXP_KALLISTO_HEADERS == ['target_id', 'length', 'eff_length', 'est_counts', 'tpm'], \
        "EXP_KALLISTO_HEADERS does not match the expected structure"

def test_EXP_RSEM_HEADERS_structure():
    assert isinstance(EXP_RSEM_HEADERS, list), "EXP_RSEM_HEADERS should be a list"
    assert EXP_RSEM_HEADERS == ['transcript_id', 'length', 'effective_length', 'expected_count', 'TPM'], \
        "EXP_RSEM_HEADERS does not match the expected structure"

def test_regex_matching():
    # Test seqid_rex1
    match = seqid_rex1.match("PB.1.2")
    assert match is not None, "seqid_rex1 failed to match valid input"
    assert match.groups() == ('1', '2'), "seqid_rex1 did not extract groups correctly"

    # Test seqid_rex2
    match = seqid_rex2.match("PB.3.4|sample")
    assert match is not None, "seqid_rex2 failed to match valid input"
    assert match.groups() == ('3', '4'), "seqid_rex2 did not extract groups correctly"

    # Test seqid_fusion
    match = seqid_fusion.match("PBfusion.5.6sample")
    assert match is not None, "seqid_fusion failed to match valid input"
    assert match.groups() == ('5', '6'), "seqid_fusion did not extract groups correctly"
