# Requantification Module Test Suite

## Overview

This document describes the test suite for the SQANTI3 requantification module, which redistributes counts from filtered artifact transcripts to rescued isoforms.

## Test Structure

### Unit Tests
Location: `test/unit/utilities/rescue/test_requant_helpers.py`

Tests individual helper functions in isolation:

- **`TestGetUnrescuedArtifacts`**: Tests identification of artifacts not in rescue table
- **`TestBuildArtifactTable`**: Tests combining rescued and unrescued artifacts
- **`TestPrepareCountMatrices`**: Tests splitting counts into base (valid) and source (artifact) matrices
- **`TestCalculateDistributionFractions`**: Tests proportional weight calculation for count distribution
- **`TestDistributeIntegerCounts`**: Tests integer count distribution with conservation
- **`TestExportCounts`**: Tests file export functionality
- **`TestCalculateTPM`**: Tests TPM (Transcripts Per Million) calculation

### Process/Integration Tests
Location: `test/process/test_requantification.py`

Tests the complete requantification pipeline:

- **`TestRedistributeCountsVectorized`**: Tests the main count redistribution pipeline
  - Simple one-to-one redistribution
  - Multi-target redistribution (one artifact → multiple isoforms)
  - Count conservation verification
  - TD (Transcript Divergency) isoform creation

- **`TestRunRequant`**: Tests the main entry point with file I/O
  - Full workflow with output file generation
  - Extended table generation

- **`TestRequantificationPipeline`**: Tests end-to-end pipeline including TPM
  - Complete pipeline from input file to all outputs
  - TPM calculation and output

- **`TestEdgeCases`**: Tests edge cases and error conditions
  - Zero-count artifacts
  - No artifacts scenario
  - All targets with zero counts

### Test Fixtures
Location: `test/process/conftest_requantification.py`

Reusable test data fixtures:

- `sample_counts_simple`: Basic count matrix (3 isoforms, 2 samples)
- `sample_counts_complex`: Complex count matrix (6 isoforms, 3 samples)
- `sample_rescue_simple`: Simple rescue mapping (1:1)
- `sample_rescue_multi_target`: Multi-target rescue mapping (1:N)
- `sample_rescue_complex`: Complex rescue scenarios
- `sample_classification_simple`: Basic classification table
- `sample_classification_complex`: Complex classification table
- `sample_counts_file`: Temporary count file for I/O tests
- Various edge case fixtures (zero counts, all artifacts, etc.)

## Running Tests

### Run All Requantification Tests
```bash
# From project root
pytest test/unit/utilities/rescue/test_requant_helpers.py test/process/test_requantification.py -v
```

### Run Only Unit Tests
```bash
pytest test/unit/utilities/rescue/test_requant_helpers.py -v
```

### Run Only Process Tests
```bash
pytest test/process/test_requantification.py -v
```

### Run Specific Test Class
```bash
pytest test/unit/utilities/rescue/test_requant_helpers.py::TestCalculateDistributionFractions -v
```

### Run with Coverage
```bash
pytest test/unit/utilities/rescue/test_requant_helpers.py test/process/test_requantification.py --cov=src.utilities.rescue.requant_helpers --cov=src.utilities.rescue.sq_requant --cov-report=html
```

## Test Scenarios Covered

### Basic Functionality
- ✅ Simple artifact-to-isoform redistribution
- ✅ Multi-target redistribution (proportional split)
- ✅ Count conservation (totals preserved)
- ✅ Integer arithmetic (no fractional counts)
- ✅ File I/O operations

### Edge Cases
- ✅ Zero-count artifacts
- ✅ No artifacts present
- ✅ All targets with zero counts (uniform split)
- ✅ All transcripts are artifacts
- ✅ Missing isoforms in count table

### Complex Scenarios
- ✅ Multiple samples
- ✅ Multiple artifacts to same target
- ✅ One artifact to multiple targets
- ✅ TD (Transcript Divergency) isoform creation
- ✅ Mixed rescued and unrescued artifacts

### Output Verification
- ✅ Reassigned counts file format
- ✅ Extended comparison table (old vs new)
- ✅ TPM calculation and output
- ✅ Proper column naming conventions

## Key Testing Principles

1. **Count Conservation**: All tests verify that total counts are preserved through redistribution
2. **Integer Arithmetic**: No fractional counts are allowed; remainders are handled explicitly
3. **Proportional Distribution**: When splitting to multiple targets, use proportional weights
4. **File Integrity**: All output files are tested for correct format and content
5. **Edge Case Handling**: Comprehensive testing of boundary conditions

## Adding New Tests

When adding new functionality to the requantification module:

1. **Unit Test**: Add test to `test_requant_helpers.py` if testing individual helper function
2. **Process Test**: Add test to `test_requantification.py` if testing integrated workflow
3. **Fixture**: Add reusable test data to `conftest_requantification.py` if needed
4. **Documentation**: Update this README with new test scenarios

### Test Template
```python
def test_new_feature(self, fixture_name):
    """Test description explaining what is being tested."""
    # Arrange: Set up test data
    input_data = ...
    
    # Act: Execute the function
    result = function_under_test(input_data)
    
    # Assert: Verify expected behavior
    assert result.meets_expectation
```

## Continuous Integration

These tests are automatically run in CI/CD pipeline on:
- Pull requests
- Commits to main branch
- Release tags

## Troubleshooting

### Test Failures

If tests fail:

1. Check that all dependencies are installed: `pip install -r requirements.txt`
2. Verify pytest is properly configured: `pytest --version`
3. Run with verbose output: `pytest -vv`
4. Check specific test with debug output: `pytest test_file.py::TestClass::test_method -vv -s`

### Common Issues

- **Import errors**: Ensure `main_path` is correctly added to `sys.path`
- **File not found**: Check that `tmp_path` fixture is being used correctly
- **Assertion failures**: Verify test data matches expected function behavior

## Related Documentation

- [Main SQANTI3 Documentation](../../wiki/Home.md)
- [Rescue Pipeline Documentation](../../wiki/Running-SQANTI3-rescue.md)
- [Test Data Documentation](../test_data/README.md)
