# Agent Directives: tests

## Role
Test suite agent context for METAINFORMANT's comprehensive testing infrastructure.

## Directory Structure
- `conftest.py` - pytest fixtures and shared test configuration
- `conftest_examples.py` - fixtures for example validation tests
- `utils.py` - shared test utilities and helper functions
- `data/` - test data fixtures organized by module
- `rna/` - RNA-specific test modules

## Test Organization
Tests follow the naming convention `test_{module}_{feature}.py`:
- `test_core_*` - Core infrastructure tests (I/O, config, paths, validation)
- `test_dna_*` - DNA sequence analysis tests
- `test_rna_*` - RNA/amalgkit workflow tests
- `test_gwas_*` - GWAS pipeline tests
- `test_visualization_*` - Plotting and visualization tests
- `test_*_comprehensive.py` - Full coverage tests for modules

## Rules and Constraints

### NO MOCKING POLICY
**CRITICAL**: This project enforces a strict NO MOCKING policy:
- All tests use REAL implementations, REAL API calls, REAL algorithms
- Never use `unittest.mock`, `pytest-mock`, or fake data
- When external services unavailable, use `pytest.skip()` not mocks
- Test fixtures contain REAL sample data, not generated placeholders

### Test Markers
Use these pytest markers appropriately:
- `@pytest.mark.slow` - Tests taking >30 seconds
- `@pytest.mark.network` - Tests requiring real network/API calls
- `@pytest.mark.external_tool` - Tests requiring CLI tools (amalgkit, muscle)
- `@pytest.mark.integration` - Cross-module integration tests

### Test Patterns
```python
# Correct: Real implementation test
def test_real_functionality(tmp_path: Path) -> None:
    result = actual_function(tmp_path / "output.json")
    assert result is not None

# Correct: Network test with graceful skip
@pytest.mark.network
def test_api_call() -> None:
    try:
        result = fetch_from_api()
    except requests.RequestException as e:
        pytest.skip(f"API unavailable: {e}")

# WRONG: Never mock
def test_bad_example(mocker):  # PROHIBITED
    mocker.patch("module.function")  # NEVER DO THIS
```

## Running Tests
```bash
# Fast tests (~15s)
bash scripts/package/test.sh

# Ultra-fast core tests (~5s)
bash scripts/package/test.sh --mode ultra-fast

# Full coverage (~2-5min)
bash scripts/package/test.sh --mode coverage

# Single file
pytest tests/test_core_io.py -v

# Pattern match
bash scripts/package/test.sh --pattern "test_rna_*"
```

## Key Files
- `NO_MOCKING_POLICY.md` - Detailed policy documentation
- `README.md` - Test suite overview
- `README_tests.md` - Additional test guidelines
