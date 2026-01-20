# Quick Start Guide: Running Rescue Tests

## What Has Been Created

### 1. Test Data Files (`test/test_data/rescue/`)
- `classification_lost_fsm.tsv` - Classification with lost FSM references
- `classification_monoexon.tsv` - Classification with monoexonic transcripts
- `classification_no_lost.tsv` - All references represented (no rescue needed)
- `classification_multi_artifact.tsv` - Multiple artifacts per reference
- `reference_mini.gtf` - Minimal GTF for testing GTF parsing functions

### 2. Unit Tests (`test/unit/utilities/rescue/`)
- `test_automatic_rescue.py` - Tests for functions in `automatic_rescue.py`
  - `get_lost_reference_id()` - 4 test cases
  - `rescue_lost_reference()` - 3 test cases
  - `rescue_fsm_monoexons()` - 4 test cases
  - `save_automatic_rescue()` - 4 test cases
  
- `test_rescue_helpers.py` - Tests for functions in `rescue_helpers.py`
  - `read_classification()` - 1 test case
  - `get_rescue_gene_targets()` - 5 test cases
  - `parse_rescue_gtf()` - 4 test cases
  - `get_rescue_reference_targets()` - 5 test cases

### 3. Process Tests (`test/process/`)
- `test_rescue_automatic.py` - Tests for `run_automatic_rescue()` workflow
  - Multi-exon rescue scenarios - 3 tests
  - Monoexon rescue modes - 4 tests
  - Multiple artifacts handling - 1 test
  - Real data files - 4 tests
  - Output format validation - 5 tests

### 4. Documentation
- `README_RESCUE_TESTS.md` - Complete testing strategy and guidelines

## Running the Tests

### Run All Rescue Tests
```bash
# From the sqanti3 root directory
pytest test/unit/utilities/rescue/ test/process/test_rescue_automatic.py -v
```

### Run Unit Tests Only
```bash
# Test automatic_rescue.py functions
pytest test/unit/utilities/rescue/test_automatic_rescue.py -v

# Test rescue_helpers.py functions
pytest test/unit/utilities/rescue/test_rescue_helpers.py -v
```

### Run Process Tests Only
```bash
pytest test/process/test_rescue_automatic.py -v
```

### Run Specific Test Class
```bash
# Example: Test only monoexon rescue
pytest test/process/test_rescue_automatic.py::TestRunAutomaticRescueMonoexon -v
```

### Run Single Test
```bash
pytest test/unit/utilities/rescue/test_automatic_rescue.py::TestGetLostReferenceId::test_some_references_lost -v
```

### Run with Coverage Report
```bash
# Generate coverage for rescue utilities
pytest test/unit/utilities/rescue/ --cov=src.utilities.rescue --cov-report=html

# Open coverage report in browser
xdg-open htmlcov/index.html  # Linux
```

## Expected Initial Results

Some tests may fail initially because:
1. The rescue functions might have minor bugs or edge cases
2. Test expectations might need adjustment based on actual function behavior
3. File I/O operations in `save_automatic_rescue()` might need mock testing

This is **expected and valuable** - failing tests help identify issues!

## Troubleshooting Common Issues

### Import Errors
If you get `ModuleNotFoundError`:
```bash
# Make sure you're running from the sqanti3 root directory
cd /home/pabloati/Programs/sqanti3

# Or add to PYTHONPATH
export PYTHONPATH=/home/pabloati/Programs/sqanti3:$PYTHONPATH
```

### Missing Dependencies
If tests fail due to missing packages (like `gtfparse`):
```bash
# Install test dependencies
pip install gtfparse pytest pytest-cov
```

### File Path Issues
If test data files can't be found:
```bash
# Verify test data exists
ls -la test/test_data/rescue/

# If files are missing, they were created in the previous steps
```

## Debugging Failed Tests

### Run with More Verbose Output
```bash
pytest test/unit/utilities/rescue/test_automatic_rescue.py -vv
```

### Run with Print Statements
```bash
pytest test/unit/utilities/rescue/test_automatic_rescue.py -v -s
```

### Stop at First Failure
```bash
pytest test/unit/utilities/rescue/ -x
```

### Run Only Failed Tests from Last Run
```bash
pytest test/unit/utilities/rescue/ --lf
```

## Next Steps After Running Tests

### 1. Review Test Results
- Note which tests pass/fail
- Understand why tests fail (bug in code or test?)
- Adjust tests or fix code as needed

### 2. Expand Test Coverage
Once `run_automatic_rescue` is fully tested, move to:
- `rescue_candidates()` - Similar pattern
- `rescue_targets()` - Similar pattern
- `run_candidate_mapping()` - More complex (needs minimap2 mocking)
- `run_rules_rescue()` / `run_ML_rescue()` - Complex workflows

### 3. Add Integration Tests
Create end-to-end tests in `test/functional/` that:
- Run complete rescue pipeline
- Use real-world size test data
- Verify final output files

### 4. Add Edge Case Tests
Based on real-world usage:
- Very large classification files
- Empty classification files
- Malformed input data
- Missing columns in classification

## Test Development Workflow

1. **Write failing test** - Define expected behavior
2. **Run test** - Confirm it fails
3. **Fix code or test** - Make it pass
4. **Refactor** - Improve code quality
5. **Repeat** - Add more tests

## Continuous Testing During Development

Use pytest's watch mode (requires pytest-watch):
```bash
pip install pytest-watch
ptw test/unit/utilities/rescue/
```

This re-runs tests automatically when you save files!

## Questions to Answer with These Tests

1. **Does `get_lost_reference_id` correctly identify lost references?**
   - What if all are lost? None are lost? Mix?

2. **Does `rescue_lost_reference` properly check for FSM?**
   - What if reference has multiple FSM? No FSM?

3. **Does monoexon rescue work correctly?**
   - Different parameter values? Edge cases?

4. **Is the rescue table generated properly?**
   - Correct columns? Correct reintroduced logic?

5. **Does the full workflow handle complex scenarios?**
   - Multiple artifacts? No lost refs? All lost?

## Success Metrics

Your test suite is successful when:
- ✅ All tests pass consistently
- ✅ Coverage > 80% for rescue utilities
- ✅ Tests catch real bugs before production
- ✅ Tests document expected behavior clearly
- ✅ Tests run fast (< 5 seconds for unit tests)

Good luck! 🚀
