# Rescue By Mapping Test Suite Documentation

## Overview

This document describes the comprehensive test suite for the `rescue_by_mapping.py` module, which implements the core logic for rescuing artifacts through mapping-based strategies.

## Test Files Created

### Main Test File
- **Location**: `test/unit/utilities/rescue/test_rescue_by_mapping.py`
- **Total Tests**: 24 tests across 7 test classes
- **Coverage**: Unit tests for individual functions and integration tests for the full pipeline

### Test Data Files
All test data files are located in `test/test_data/rescue/`:

1. **`mapping_hits.tsv`** - Simulated mapping results between rescue candidates and targets
   - Contains 12 mapping records
   - Includes various scenarios: single hits, multiple hits, different alignment scores
   - Columns: `rescue_candidate`, `mapping_hit`, `alignment_type`, `alignment_score`

2. **`reference_rules.tsv`** - Reference classification results using rules-based filtering
   - Contains 12 reference transcripts with filter results
   - Columns: `isoform`, `filter_result` (Isoform/Artifact)

3. **`reference_ml.tsv`** - Reference classification results using ML-based filtering
   - Contains 12 reference transcripts with ML probabilities
   - Columns: `isoform`, `POS_MLprob` (0.0-1.0)

4. **`classification.tsv`** - Long-read transcript classification data
   - Contains 11 transcript entries
   - Includes structural categories, filter results, and gene associations
   - Columns: `isoform`, `structural_category`, `filter_result`, `associated_transcript`, `associated_gene`, `POS_MLprob`

5. **`automatic_rescue.tsv`** - Pre-existing automatic rescue results
   - Contains 2 automatically rescued transcripts
   - Columns: `artifact`, `assigned_transcript`, `rescue_mode`, `origin`, `reintroduced`

## Test Classes and Coverage

### 1. TestMergeClassifications (4 tests)
Tests the `merge_classifications()` function which combines reference and long-read classifications.

**Tests:**
- `test_merge_rules_strategy` - Verifies rules-based merging preserves both reference and LR entries
- `test_merge_ml_strategy` - Tests ML strategy with threshold-based filtering (0.7 threshold)
- `test_merge_ml_different_threshold` - Tests ML with different threshold (0.8)
- `test_merge_preserves_lr_classification` - Ensures LR classification data is preserved correctly

**Key Assertions:**
- Correct creation of `origin` column (reference vs lr_defined)
- Filter results based on ML probability thresholds
- Preservation of original classification data

### 2. TestAddFilterResults (3 tests)
Tests the `add_filter_results()` function which adds filter results to mapping hits.

**Tests:**
- `test_add_filter_results_basic` - Verifies basic column addition (hit_filter_result, hit_origin, etc.)
- `test_add_filter_results_correct_mapping` - Checks correct mapping of filter results to hits
- `test_add_filter_results_with_ml` - Tests addition of ML probabilities

**Key Assertions:**
- New columns are added correctly
- Filter results match expected values
- Data integrity is maintained

### 3. TestSelectBestHits (3 tests)
Tests the `select_best_hits()` function which selects optimal hits from multiple alignments.

**Tests:**
- `test_select_highest_score` - Verifies selection of highest scoring hit
- `test_prefer_lr_defined` - Tests preference for lr_defined over reference when scores are equal
- `test_keep_ties` - Ensures ties (equal scores) are kept

**Key Assertions:**
- Correct score-based selection
- Origin preference logic works
- Tie handling is appropriate

### 4. TestFilterMappingHits (3 tests)
Tests the `filter_mapping_hits()` function which filters and scores mapping hits.

**Tests:**
- `test_filter_rules_strategy` - Tests rules-based filtering (only Isoforms kept)
- `test_filter_ml_strategy` - Tests ML-based scoring (alignment_score × ML_prob)
- `test_filter_selects_best_per_candidate` - Verifies best hit selection per candidate

**Key Assertions:**
- Filtering removes Artifacts
- Score calculation is correct
- Best hits are selected per rescue candidate

### 5. TestMergeRescueModes (3 tests)
Tests the `merge_rescue_modes()` function which combines automatic and mapping-based rescue results.

**Tests:**
- `test_merge_basic` - Basic merging of automatic and full rescue DataFrames
- `test_merge_marks_reintroduced` - Tests reintroduction marking (yes/no)
- `test_merge_adds_rescue_mode` - Verifies rescue_mode column is added correctly

**Key Assertions:**
- Both rescue modes are present in result
- Reintroduction status is correct
- Column structure is as expected

### 6. TestFindReintroducedTranscripts (3 tests)
Tests the `find_reintroduced_transcripts()` function which identifies reference transcripts to reintroduce.

**Tests:**
- `test_find_with_automatic_rescue` - Tests with existing automatic rescue
- `test_find_without_automatic_rescue` - Tests with empty automatic rescue
- `test_find_removes_duplicates` - Verifies duplicate removal

**Key Assertions:**
- Only reference transcripts are included
- Automatic rescue is incorporated
- Duplicates are properly removed

### 7. TestRescueByMapping (5 integration tests)
Integration tests for the complete `rescue_by_mapping()` pipeline.

**Tests:**
- `test_rescue_by_mapping_rules` - Full pipeline with rules strategy
- `test_rescue_by_mapping_ml` - Full pipeline with ML strategy (threshold 0.7)
- `test_rescue_by_mapping_different_threshold` - ML with varying thresholds (0.9 vs 0.6)
- `test_rescue_by_mapping_empty_automatic` - Pipeline without automatic rescue
- `test_rescue_by_mapping_preserves_automatic` - Ensures automatic rescue is preserved

**Key Assertions:**
- Inclusion list is generated correctly
- Rescue DataFrame has correct structure
- Rescue modes are properly assigned
- Automatic rescue entries are preserved
- Threshold effects on rescue count

## Running the Tests

### Run All Tests
```bash
pytest test/unit/utilities/rescue/test_rescue_by_mapping.py -v
```

### Run Specific Test Class
```bash
pytest test/unit/utilities/rescue/test_rescue_by_mapping.py::TestMergeClassifications -v
```

### Run Specific Test
```bash
pytest test/unit/utilities/rescue/test_rescue_by_mapping.py::TestRescueByMapping::test_rescue_by_mapping_rules -v
```

### Run with Coverage
```bash
pytest test/unit/utilities/rescue/test_rescue_by_mapping.py --cov=src.utilities.rescue.rescue_by_mapping --cov-report=html
```

## Test Data Design Principles

1. **Reusability**: Test data files are designed to be used across multiple test scenarios
2. **Completeness**: Data covers edge cases (ties, missing values, multiple alignments)
3. **Clarity**: File structure is simple and readable for easy debugging
4. **Realism**: Data reflects actual rescue pipeline scenarios

## Key Testing Patterns

### Fixture Usage
All tests use pytest fixtures for:
- Test data directory paths
- Loading DataFrames from files
- Maintaining consistent test data across tests

### Test Independence
- Each test is independent and can run in any order
- Fixtures ensure fresh data for each test
- No shared mutable state between tests

### Assertion Strategy
- Tests verify both structure (columns, types) and content (values, relationships)
- Edge cases are explicitly tested (empty data, ties, missing values)
- Integration tests verify end-to-end functionality

## Coverage Summary

| Component | Unit Tests | Integration Tests | Total |
|-----------|-----------|-------------------|-------|
| merge_classifications | 4 | - | 4 |
| add_filter_results | 3 | - | 3 |
| select_best_hits | 3 | - | 3 |
| filter_mapping_hits | 3 | - | 3 |
| merge_rescue_modes | 3 | - | 3 |
| find_reintroduced_transcripts | 3 | - | 3 |
| rescue_by_mapping | - | 5 | 5 |
| **Total** | **19** | **5** | **24** |

## Notes

- Tests currently pass for rules-based strategy (18/24 tests passing)
- ML-based tests require the classification DataFrame to include `POS_MLprob` column
- Integration tests verify the complete rescue pipeline from mapping hits to final inclusion list
- Test data is stored in TSV format for easy inspection and modification
