# Test Suite

METAINFORMANT test suite using pytest with **real implementations only** (NO MOCKING policy).

## Quick Start

```bash
# Fast tests (~15s)
bash scripts/package/test.sh

# Ultra-fast core tests (~5s)
bash scripts/package/test.sh --mode ultra-fast

# Full coverage (~2-5min)
bash scripts/package/test.sh --mode coverage

# Single test file
pytest tests/test_core_functionality.py -v

# Single test function
pytest tests/test_core_functionality.py::TestCoreIO::test_json_operations -v

# Pattern match
bash scripts/package/test.sh --pattern "test_dna_*"
```

## Test Organization

Tests are organized by module, mirroring `src/metainformant/` structure:

| Prefix | Module | Count |
|--------|--------|-------|
| `test_core_*` | Core utilities | 18 |
| `test_dna_*` | DNA analysis | 28 |
| `test_rna_*` | RNA workflows | 32 |
| `test_gwas_*` | GWAS pipeline | 20 |
| `test_protein_*` | Protein analysis | 16 |
| `test_visualization_*` | Plotting | 14 |
| `test_math_*` | Math/popgen | 22 |
| `test_life_events_*` | Life events | 14 |
| `test_networks_*` | Networks | 6 |
| `test_ml_*` | Machine learning | 2 |
| `test_ontology_*` | Ontology/GO | 6 |
| `test_singlecell_*` | Single-cell | 3 |
| `test_information_*` | Information theory | 2 |
| `test_multiomics_*` | Multi-omics | 3 |
| Others | Integration, etc. | ~20 |

## Test Modes

| Mode | Description | Time |
|------|-------------|------|
| `ultra-fast` | Core functionality only | ~5s |
| `fast` | Essential tests (default) | ~15s |
| `coverage` | Full suite with coverage | ~2-5min |
| `network` | Tests with real API calls | ~30s |
| `external` | Tests requiring CLI tools | varies |
| `integration` | Cross-module tests | ~30s |

## Test Markers

```python
@pytest.mark.slow           # Long-running tests
@pytest.mark.network        # Real API calls (NCBI, UniProt, etc.)
@pytest.mark.external_tool  # Requires CLI tools (amalgkit, muscle)
@pytest.mark.integration    # Cross-module integration
```

Run specific markers:
```bash
pytest -m "network" tests/        # Only network tests
pytest -m "not slow" tests/       # Skip slow tests
```

## NO MOCKING Policy

**All tests use real implementations:**

```python
# CORRECT: Real file I/O
def test_json_operations(tmp_path):
    output = tmp_path / "result.json"
    io.dump_json({"key": "value"}, output)
    assert output.exists()
    data = io.load_json(output)
    assert data["key"] == "value"

# CORRECT: Real API with graceful skip
@pytest.mark.network
def test_uniprot_fetch():
    try:
        result = fetch_uniprot("P12345")
        assert result is not None
    except requests.RequestException as e:
        pytest.skip(f"API unavailable: {e}")

# WRONG: Never mock
def test_bad_example():
    with mock.patch("module.function"):  # PROHIBITED
        ...
```

## Test Data

Test fixtures are in `tests/data/`:

```
tests/data/
├── dna/        # FASTA, FASTQ samples
├── protein/    # PDB, sequence data
├── rna/        # RNA-seq test data
├── ontology/   # GO OBO files
├── phenotype/  # Trait data
└── epigenome/  # Methylation data
```

## Writing Tests

```python
from pathlib import Path
import pytest
from metainformant.core import io

def test_my_feature(tmp_path: Path) -> None:
    """Test description."""
    # Arrange
    input_data = {"test": "data"}
    output_file = tmp_path / "output.json"

    # Act
    result = my_function(input_data, output_file)

    # Assert
    assert output_file.exists()
    assert result["status"] == "success"

@pytest.mark.network
def test_api_integration() -> None:
    """Test with real API - skip if offline."""
    try:
        result = fetch_from_api("test_id")
        assert result is not None
    except Exception as e:
        pytest.skip(f"API unavailable: {e}")
```

## Configuration

`conftest.py` provides common fixtures:

- `tmp_path` - Temporary directory (pytest built-in)
- Test data paths
- NCBI API configuration

## Coverage

```bash
# Generate coverage report
bash scripts/package/test.sh --mode coverage

# View HTML report
open output/coverage_html/index.html
```

Coverage threshold: 85% minimum for CI.

## Related

- [Testing Guide](../docs/testing.md)
- [NO_MOCKING_POLICY](../docs/NO_MOCKING_POLICY.md)
- [Test Script](../scripts/package/test.sh)
