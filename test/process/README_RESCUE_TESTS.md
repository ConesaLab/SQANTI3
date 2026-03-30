# Rescue Module Test Suite

## Overview
This document outlines the test strategy for the rescue module, focusing on `run_automatic_rescue` and its component functions.

## Test Strategy

### Testing Levels

#### 1. **Unit Tests** (`test/unit/utilities/`)
Test individual helper functions in isolation with simple inputs:
- Functions from `automatic_rescue.py`
- Functions from `rescue_helpers.py`
- Small, focused tests with minimal dependencies

#### 2. **Process Tests** (`test/process/`)
Test complete workflows with realistic data structures:
- `run_automatic_rescue` full workflow
- Integration of multiple functions
- Uses fixtures similar to QC module

## Test Plan for `run_automatic_rescue`

### Component Functions to Unit Test

#### From `automatic_rescue.py`:

1. **`get_lost_reference_id(df)`**
   - Input: DataFrame with FSM/ISM classifications
   - Output: Array of reference IDs without "Isoform" filter_result
   - Test cases:
     - All references represented → empty array
     - Some references lost → returns lost IDs
     - No references at all → empty array
     - Mixed FSM/ISM with some lost

2. **`rescue_lost_reference(ref_id, classif)`**
   - Input: Reference ID string, classification DataFrame
   - Output: DataFrame with rescued isoform or None
   - Test cases:
     - Reference with FSM → returns reference ID
     - Reference without FSM → returns None
     - Reference not in classification → returns None

3. **`rescue_fsm_monoexons(df)`**
   - Input: Full classification DataFrame
   - Output: DataFrame of rescued monoexonic FSM transcripts
   - Test cases:
     - Monoexon FSM with lost reference → rescued
     - Monoexon FSM already represented → not rescued
     - No monoexons → empty DataFrame
     - Multiple monoexons for same gene

4. **`save_automatic_rescue(inclusion_df, class_df, prefix)`**
   - Input: Inclusion list, classification, output prefix
   - Output: Rescue table DataFrame
   - Test cases:
     - Normal rescue → proper table format
     - Empty inclusion list → "none" handling
     - Multiple artifacts to same reference → reintroduced column logic
     - Verify column names and types

#### From `rescue_helpers.py`:

5. **`read_classification(filename)`**
   - Simple wrapper, minimal testing needed
   - Test: Can read valid classification file

6. **`get_rescue_gene_targets(df, rescue_candidates)`**
   - Input: Classification DataFrame, list of candidate IDs
   - Output: Array of unique gene IDs
   - Test cases:
     - Single candidate → one gene
     - Multiple candidates, same gene → one gene
     - Multiple candidates, different genes → multiple genes
     - Candidate not in DataFrame → empty

7. **`parse_rescue_gtf(gtf_file)`**
   - Input: GTF file path
   - Output: Dictionary of gene_id → [transcript_ids]
   - Test cases:
     - Simple GTF → correct mapping
     - Gene with multiple transcripts → all transcripts listed
     - GTF with no transcript_id → handled gracefully

8. **`get_rescue_reference_targets(ref_gtf, target_genes)`**
   - Input: GTF path, list of gene IDs
   - Output: Series of transcript IDs
   - Test cases:
     - Target genes in GTF → transcripts returned
     - Target genes not in GTF → empty Series
     - Multiple transcripts per gene → all returned

### Process Test for `run_automatic_rescue`

Test the complete workflow with realistic data:

**Input Fixtures:**
- `sample_rescue_classification`: DataFrame with FSM, ISM, artifacts
- `monoexon_parameter`: "all", "fsm", or None
- `output_prefix`: Temporary directory path

**Test Scenarios:**

1. **Happy Path - FSM Multi-exon Rescue**
   - Setup: Classification with lost FSM references
   - Expected: Rescued transcripts returned, rescue_df created
   - Verify: Correct transcripts in inclusion list, table columns

2. **No Lost References**
   - Setup: All references represented by isoforms
   - Expected: Empty inclusion list, "none" in rescue_df
   - Verify: Log message "No lost references found"

3. **Monoexon Rescue - All**
   - Setup: Classification with lost monoexon FSM
   - Parameter: monoexons="all"
   - Expected: Monoexons included in rescue
   - Verify: Monoexons in final inclusion list

4. **Monoexon Rescue - FSM Only**
   - Setup: Classification with monoexon FSM and ISM
   - Parameter: monoexons="fsm"
   - Expected: Only FSM monoexons rescued
   - Verify: ISM monoexons not in inclusion list

5. **Monoexon Rescue - Disabled**
   - Setup: Classification with monoexon FSM
   - Parameter: monoexons=None
   - Expected: No monoexons rescued
   - Verify: Empty or no monoexons in inclusion list

6. **Multiple Artifacts per Reference**
   - Setup: Multiple ISM artifacts pointing to same reference
   - Expected: Reference rescued once, "reintroduced" column correct
   - Verify: First artifact marked "yes", others "no"

## Test Data Requirements

### Sample Classification DataFrame Schema
Required columns:
```python
[
    'isoform',
    'structural_category',  # FSM, ISM, NIC, NNC
    'exons',
    'filter_result',  # "Isoform" or "Artifact"
    'associated_transcript',
    'associated_gene'
]
```

### Test Data Scenarios

Create small CSV fixtures in `test/test_data/rescue/`:

1. **`classification_lost_fsm.tsv`**
   - 5 FSM with >1 exon (filter_result="Isoform")
   - 3 FSM with >1 exon (filter_result="Artifact", lost reference)
   - Mix of represented and lost references

2. **`classification_monoexon.tsv`**
   - FSM monoexons (filter_result="Artifact")
   - ISM monoexons (filter_result="Artifact")
   - Some represented, some lost

3. **`classification_no_lost.tsv`**
   - All references represented
   - Mix of FSM/ISM but all filter_result="Isoform"

4. **`classification_multi_artifact.tsv`**
   - Multiple ISM artifacts → same lost reference
   - Test "reintroduced" logic

5. **`reference_mini.gtf`**
   - Small GTF with 3-5 genes, 10-15 transcripts
   - For testing GTF parsing functions

## File Structure

```
test/
├── unit/
│   └── utilities/
│       └── rescue/
│           ├── __init__.py
│           ├── test_automatic_rescue.py      # Unit tests for automatic_rescue.py functions
│           └── test_rescue_helpers.py        # Unit tests for rescue_helpers.py functions
├── process/
│   ├── conftest.py                           # Update with rescue fixtures
│   └── test_rescue_automatic.py              # Process test for run_automatic_rescue
└── test_data/
    └── rescue/
        ├── classification_lost_fsm.tsv
        ├── classification_monoexon.tsv
        ├── classification_no_lost.tsv
        ├── classification_multi_artifact.tsv
        └── reference_mini.gtf
```

## Implementation Order

1. **Phase 1: Set up test data**
   - Create `test/test_data/rescue/` directory
   - Create minimal test TSV and GTF files

2. **Phase 2: Unit tests**
   - Start with `test_automatic_rescue.py`
   - Then `test_rescue_helpers.py`
   - Focus on one function at a time

3. **Phase 3: Fixtures**
   - Add rescue fixtures to `test/process/conftest.py`
   - Create `sample_rescue_classification` fixture
   - Reuse patterns from QC fixtures

4. **Phase 4: Process test**
   - Create `test_rescue_automatic.py`
   - Test full `run_automatic_rescue` workflow
   - Use pytest parametrize for monoexon variations

5. **Phase 5: Expand coverage**
   - Add edge cases discovered during testing
   - Test error handling
   - Add integration tests if needed

## Running Tests

```bash
# Run all rescue unit tests
pytest test/unit/utilities/rescue/ -v

# Run rescue process tests
pytest test/process/test_rescue_automatic.py -v

# Run specific test
pytest test/unit/utilities/rescue/test_automatic_rescue.py::test_get_lost_reference_id_some_lost -v

# Run with coverage
pytest test/unit/utilities/rescue/ --cov=src.utilities.rescue --cov-report=html
```

## Key Testing Patterns

### 1. DataFrame Fixtures
```python
@pytest.fixture
def sample_rescue_classification():
    """Classification DataFrame with lost FSM references."""
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.1.2', 'REF1', 'PB.2.1'],
        'structural_category': ['full-splice_match', 'full-splice_match', 
                               'full-splice_match', 'incomplete-splice_match'],
        'exons': [2, 2, 2, 2],
        'filter_result': ['Isoform', 'Artifact', 'Artifact', 'Artifact'],
        'associated_transcript': ['REF1', 'REF1', 'REF1', 'REF2'],
        'associated_gene': ['GENE1', 'GENE1', 'GENE1', 'GENE2']
    })
```

### 2. Parameterized Tests
```python
@pytest.mark.parametrize("monoexons,expected_count", [
    ("all", 3),
    ("fsm", 2),
    (None, 1)
])
def test_run_automatic_rescue_monoexon_modes(
    sample_rescue_classification, monoexons, expected_count
):
    # Test logic
    pass
```

### 3. Temporary Output Files
```python
def test_save_automatic_rescue(tmp_path, sample_rescue_classification):
    prefix = tmp_path / "test_output"
    result_df = save_automatic_rescue(inclusion_df, class_df, str(prefix))
    
    # Verify output files created
    assert (tmp_path / "test_output_automatic_rescue_table.tsv").exists()
```

## Next Steps

After completing `run_automatic_rescue` tests, expand to:
- `rescue_candidates()` 
- `rescue_targets()`
- `run_candidate_mapping()`
- `run_rules_rescue()` / `run_ML_rescue()`
- `save_rescue_results()`

Each function should follow the same pattern:
1. Unit tests for helper functions
2. Process tests for main workflows
3. Proper fixtures and test data
