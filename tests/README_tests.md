# Tests Directory

This directory contains the comprehensive test suite for METAINFORMANT, organized by module and functionality.

## Test Organization

### Module-Based Tests
Tests are organized to mirror the source code structure:

```
tests/
├── test_core_*.py          # Core infrastructure tests
├── test_dna_*.py           # DNA analysis tests
├── test_rna_*.py           # RNA workflow tests
├── test_protein_*.py       # Protein analysis tests
├── test_math_*.py          # Mathematical biology tests
├── test_simulation_*.py    # Simulation tests
├── test_visualization_*.py # Visualization tests
└── test_*.py               # Cross-cutting and integration tests
```

### Test Data (`data/`)
Organized test fixtures and sample data by domain.

## Running Tests

### All Tests
```bash
# Run complete test suite
python -m pytest tests/ -v

# With coverage reporting
python -m pytest tests/ --cov=src/metainformant --cov-report=html
```

### Module-Specific Tests
```bash
# Core infrastructure
python -m pytest tests/test_core_*.py

# DNA analysis
python -m pytest tests/test_dna_*.py

# RNA workflows
python -m pytest tests/test_rna_*.py
```

## Test Development Guidelines

### Writing Tests
```python
def test_function_expected_behavior(tmp_path: Path) -> None:
    """Test description following repository style."""
    # Arrange: Set up test data and preconditions
    input_data = prepare_test_data()

    # Act: Execute the function under test
    result = module_under_test.function(input_data)

    # Assert: Verify expected behavior
    assert result.property == expected_value
    assert validates_result(result)
```

### Test Naming Conventions
- **File**: `test_{domain}_{module}.py` or `test_{domain}_{module}_{aspect}.py`
- **Function**: `test_{function_name}_{expected_behavior}`
- **Class**: `Test{ClassName}` for class-based tests

### Test Best Practices
- **Isolation**: Each test should be independent
- **Clarity**: Clear test names and assertions
- **Coverage**: Test both success and failure cases
- **Documentation**: Include docstrings explaining test purpose
- **Fixtures**: Use pytest fixtures for shared test data

## Related Documentation

- See individual module README files for module-specific testing details
- See `NO_MOCKING_POLICY.md` for testing philosophy and guidelines
- See main project documentation for development workflow
- See `scripts/run_tests.sh` for test execution tooling

This test suite ensures the reliability, correctness, and performance of all METAINFORMANT functionality.
