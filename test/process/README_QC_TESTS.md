# QC Computations Test Suite

## Overview
This document describes the process-level test infrastructure created for testing `qc_computations.py` functions.

## Test Strategy
The tests are classified as **process tests** rather than unit tests because:
1. They require complex `isoforms_info` objects from the classification pipeline
2. They test full workflow steps in the QC pipeline
3. They use realistic test data and fixtures

## Test Infrastructure

### Fixtures (`test/process/conftest.py`)
Shared fixtures provide reusable test data:

- **`sample_isoforms_info`**: Creates a dictionary of 5 realistic `myQueryTranscripts` objects matching test isoforms:
  - PB.124830.1 (novel_in_catalog, 2 exons)
  - PB.103714.1 (full-splice_match, 1 exon)
  - PB.103724.1 (novel_not_in_catalog, 4 exons)
  - PB.103781.1 (full-splice_match, 1 exon)
  - PB.23068.2 (full-splice_match, 2 exons)

- **`fields_class_cur`**: Base classification field list that gets extended by QC functions

- **File path fixtures**: Point to organized test data:
  - `fl_count_file_single`: Single-sample abundance file
  - `fl_count_file_multi`: Multi-sample abundance file
  - `junction_file`: Junction data file
  - `genome_file`: Test genome FASTA

### Test File (`test/process/test_qc_computations.py`)
Comprehensive tests for all major QC computation functions:

## Test Coverage

### 1. `ratio_TSS_dict_reading` (4 tests)
- ✅ Basic ratio TSS assignment
- ✅ Missing isoforms get default value of 1
- ✅ NaN values converted to None
- ✅ Extra isoforms in file handled gracefully

### 2. `process_rts` (2 tests)
- ✅ RTS detection results properly assigned (TRUE/FALSE)
- ✅ All isoforms marked FALSE when no RTS detected
- Uses mocking to avoid complex RTS pipeline dependencies

### 3. `classify_fsm` (3 tests)
- ✅ Single isoform per gene → FSM_class "A"
- ✅ Gene with full-splice_match → FSM_class "C"
- ✅ Gene without full-splice_match → FSM_class "B"

### 4. `full_length_quantification` (4 tests)
- ✅ Single-sample FL count format
- ✅ Multi-sample FL count format with field extension
- ✅ Missing isoforms assigned FL = 0
- ✅ Other isoform attributes preserved

### 5. `isoforms_junctions` (7 tests)
- ✅ Canonical junction detection
- ✅ Bite junction detection
- ✅ Indel near junction counting
- ✅ Minimum coverage calculation
- ✅ Minimum sample coverage calculation
- ✅ Standard deviation of junction coverage
- ✅ Isoforms without junctions remain unmodified

## Running the Tests

```bash
# Run all QC computation tests
pytest test/process/test_qc_computations.py -v

# Run specific test class
pytest test/process/test_qc_computations.py::TestFullLengthQuantification -v

# Run specific test
pytest test/process/test_qc_computations.py::TestClassifyFSM::test_classify_fsm_single_isoform_per_gene -v
```

## Test Data Organization
Test data is now organized by type in `test/test_data/`:
- `abundance/`: FL count files (8 test files with various formats)
- `reference/`: Reference annotations
- `isoforms/`: Isoform files and classifications
- `orfs/`: ORF prediction files
- `junctions/`: Junction data
- `genome/`: Genome FASTA files
- `bedfiles/`: BED format files
- `bam/`: BAM alignment files
- `other/`: Miscellaneous test files

## Key Testing Patterns

### 1. Fixture-Based Testing
```python
def test_fl_count_single_sample(self, sample_isoforms_info, fl_count_file_single, fields_class_cur):
    result_info, result_fields = full_length_quantification(
        fl_count_file_single, sample_isoforms_info, fields_class_cur.copy()
    )
    assert result_info["PB.124830.1"].FL == 150
```

### 2. Mock Junction Readers
```python
@pytest.fixture
def mock_junction_reader(self):
    return [
        {'isoform': 'PB.124830.1', 'canonical': 'canonical', ...},
        {'isoform': 'PB.103724.1', 'canonical': 'non_canonical', ...}
    ]
```

### 3. Monkeypatching for Complex Dependencies
```python
def test_rts_detection(self, monkeypatch):
    def mock_rts(args, genome_dict):
        return {"PB.124830.1": [(100, 200)]}
    monkeypatch.setattr(src.qc_computations, "rts", mock_rts)
```

## Benefits of This Approach

1. **Maintainability**: Shared fixtures reduce code duplication
2. **Realistic Testing**: Uses actual data structures from the pipeline
3. **Fast Execution**: All 20 tests run in ~0.5 seconds
4. **Clear Intent**: Each test focuses on one specific behavior
5. **Easy Extension**: New tests can reuse existing fixtures

## Future Extensions

To test additional functions, simply:
1. Add test class to `test_qc_computations.py`
2. Reuse existing fixtures or add new ones to `conftest.py`
3. Follow the established testing patterns

## Summary

**Total Tests**: 20
**Test Coverage**: All major QC computation functions
**Test Type**: Process/Integration tests
**Execution Time**: ~0.5 seconds
**Status**: ✅ All passing
