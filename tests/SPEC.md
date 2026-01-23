# Specification: tests

## Scope

Comprehensive pytest test suite for METAINFORMANT. Contains real-implementation tests covering all 19 domain modules with strict adherence to the NO MOCKING policy. Tests validate functionality using actual algorithms, real file I/O, and genuine API calls.

## Architecture

- **Dependency Level**: Test
- **Component Type**: Test Suite
- **Framework**: pytest with custom markers and fixtures

### Test Hierarchy
```
tests/
├── conftest.py              # Global fixtures and pytest configuration
├── conftest_examples.py     # Fixtures for example validation
├── data/                    # Real test data fixtures (not mocks)
├── rna/                     # RNA-specific test modules
└── test_{module}_{feature}.py  # Module-specific tests
```

## Data Structures

### Test File Types
- **test_core_*.py**: Core infrastructure tests (I/O, config, paths, logging, validation)
- **test_dna_*.py**: DNA sequence analysis, alignment, population genetics
- **test_rna_*.py**: RNA-seq workflows, amalgkit integration
- **test_gwas_*.py**: GWAS pipelines, association studies
- **test_*_comprehensive.py**: Full coverage tests for modules
- **conftest.py**: pytest fixtures providing real sample data and temp directories

### Test Markers (pytest.mark)
- `slow`: Tests taking >30 seconds
- `network`: Tests requiring real network/API calls
- `external_tool`: Tests requiring CLI tools (amalgkit, muscle, bcftools)
- `integration`: Cross-module integration tests

## Interface

### Running Tests
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

# Specific test
pytest tests/test_core_functionality.py::TestCoreIO::test_json_operations -v
```

### Writing Tests
```python
# Real implementation test pattern
def test_real_functionality(tmp_path: Path) -> None:
    """Test with real implementation."""
    result = actual_function(tmp_path / "output.json")
    assert result is not None
    assert (tmp_path / "output.json").exists()

# Network test with graceful skip
@pytest.mark.network
def test_api_call() -> None:
    try:
        result = fetch_from_api()
        assert result is not None
    except requests.RequestException as e:
        pytest.skip(f"API unavailable: {e}")
```

### Policy Enforcement
- **NO MOCKING**: Never use unittest.mock, pytest-mock, or fake data
- **Real Fixtures**: All test data in data/ contains actual biological samples
- **Graceful Skips**: Use pytest.skip() when external dependencies unavailable
