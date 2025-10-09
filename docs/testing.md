# Testing

## Testing Policy (STRICTLY NO MOCKS/FAKES)

**ABSOLUTE PROHIBITION**: Never use fake/mocked/stubbed methods, objects, or network shims in code or tests.

**Real Implementation Only**: Tests must exercise real code paths and real external behavior without exceptions:
- **Networked tests**: perform real HTTP requests with short timeouts. If offline, skip gracefully with clear messages.
- **CLI-dependent tests** (e.g., amalgkit): run only when the dependency is available on PATH; otherwise skip with dependency notes.
- **Database tests**: use real database connections or skip when unavailable.
- **API tests**: make real API calls or skip when network/credentials unavailable.

**Anti-Pattern Enforcement**: Mocking is an anti-pattern that creates brittle tests disconnected from reality.

**Quality Assurance**: Real implementations reveal actual bugs, performance issues, and integration problems.

**Environment Setup**: It is acceptable to set environment variables for test setup, but do not monkeypatch or replace functions.

**Test Artifacts**: Tests must write all artifacts only under `output/` directory.

**Reproducibility**: Prefer deterministic seeds and stable filenames for reproducible test runs.

**Clear Documentation**: When external dependencies are unavailable, tests must clearly document what is being skipped and why.

## Running Tests

### Basic Execution

```bash
# Run all tests with minimal output
uv run pytest -q

# Run with coverage
uv run pytest --cov=src/metainformant --cov-report=html

# Run specific test files
uv run pytest tests/test_core_text.py -v
```

### Professional Test Runner

Use the comprehensive test runner script:

```bash
# Standard test run
./scripts/run_tests.sh

# Fast tests only (skip slow/network/external dependencies)
./scripts/run_tests.sh --fast

# Include network tests
./scripts/run_tests.sh --network

# Generate coverage reports without running tests
./scripts/run_tests.sh --report-only

# Run tests matching a pattern
./scripts/run_tests.sh --pattern "core_*"
```

### Test Categories

Tests are organized with custom pytest markers:

- **`slow`**: Tests that take significant time to complete
- **`network`**: Tests requiring internet connectivity
- **`external_tool`**: Tests requiring external CLI tools (e.g., amalgkit)
- **`integration`**: End-to-end integration tests

### Configuration

Test configuration is managed in `pyproject.toml`:

```toml
[tool.pytest.ini_options]
addopts = ["--cov=src/metainformant"]
markers = [
    "slow: marks tests as slow",
    "network: marks tests as requiring network access",
    "external_tool: marks tests as requiring external tools",
    "integration: marks tests as integration tests",
]

[tool.coverage.run]
branch = true
omit = ["*/tests/*", "*/__pycache__/*"]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise AssertionError",
    "raise NotImplementedError",
]
```

## Test Structure and Organization

### Directory Layout

- **`tests/`**: All test files at repository root
- **`tests/conftest.py`**: Shared pytest configuration and fixtures
- **`tests/test_*.py`**: Individual test modules
- **`tests/data/`**: Test input fixtures and data
- **`output/`**: All test outputs and artifacts

### Test Coverage Matrix

The test suite provides comprehensive coverage:

| Module | Primary Tests | Coverage |
|--------|---------------|----------|
| Core Utilities | `test_core_*.py` | Config, I/O, text, logging, parallel, paths, cache, hash, db |
| DNA Analysis | `test_dna_*.py` | Sequences, alignment, MSA, phylogeny, population genetics, FASTQ |
| RNA Analysis | `test_rna_*.py` | amalgkit wrapper, workflow, configs, step runners |
| Single-Cell | `test_singlecell_*.py` | Preprocessing, dimensionality, clustering, trajectory |
| Quality Control | `test_quality_*.py` | FASTQ quality analysis, contamination detection |
| Protein Analysis | `test_protein_*.py` | UniProt, PDB, AlphaFold, InterPro, structure analysis |
| Mathematical Models | `test_math_*.py` | Population genetics, coalescent, epidemiology, selection |
| Simulation | `test_simulation.py` | Sequence generation, RNA counts, agent-based models |
| Visualization | `test_visualization_*.py` | Plots, trees, animations |
| Ontology | `test_ontology_*.py` | Gene Ontology, OBO parsing |
| CLI Interface | `test_cli.py` | Command-line argument parsing and dispatch |

### Test Conventions

#### File Organization
- **One test file per source module**: `src/foo/bar.py` → `tests/test_foo_bar.py`
- **Class-based organization**: Group related tests in classes
- **Descriptive test names**: Use clear, descriptive function names

#### Test Isolation
- **Setup/teardown**: Use `setup_method` and `teardown_method` for clean environments
- **Independent tests**: Each test should be runnable in isolation
- **Output directory**: All test artifacts go to `output/` subdirectories

#### Example Test Structure

```python
import pytest
from metainformant.core.io import ensure_directory
from metainformant.some_module import some_function

class TestSomeFunction:
    def setup_method(self):
        """Setup test environment."""
        self.output_dir = ensure_directory("output/test_some_module")
        self.test_data = {...}

    def teardown_method(self):
        """Cleanup after test."""
        # Cleanup if needed (usually not required)
        pass

    def test_basic_functionality(self):
        """Test basic function behavior."""
        result = some_function(self.test_data)
        assert result is not None

    def test_edge_cases(self):
        """Test edge cases and error conditions."""
        with pytest.raises(ValueError):
            some_function(invalid_input)

    @pytest.mark.slow
    def test_performance_intensive(self):
        """Test that takes significant time."""
        # Test implementation
        pass

    @pytest.mark.network
    def test_network_dependent(self):
        """Test requiring network access."""
        # Skip if offline
        try:
            result = some_network_function()
            assert result is not None
        except ConnectionError:
            pytest.skip("Network not available")
```

## Network and External Tool Testing

### Network Tests

Network-dependent tests use real API calls with graceful offline handling:

```python
import requests
import pytest

def _check_online():
    """Check if network connectivity is available."""
    try:
        response = requests.get('https://httpbin.org/status/200', timeout=5)
        return response.status_code == 200
    except requests.RequestException:
        return False

@pytest.mark.network
@pytest.mark.skipif(not _check_online(), reason="Network not available")
def test_uniprot_api():
    """Test UniProt API with real network call."""
    from metainformant.protein.uniprot import map_ids_uniprot
    
    result = map_ids_uniprot(['P12345'], from_db='UniProtKB_AC-ID', to_db='Gene_Name')
    assert len(result) >= 0  # May be empty if ID not found
```

### External Tool Tests

Tests requiring external CLI tools check for availability:

```python
import shutil
import subprocess
import pytest

@pytest.mark.external_tool
@pytest.mark.skipif(not shutil.which("amalgkit"), reason="amalgkit not available")
def test_amalgkit_execution():
    """Test amalgkit CLI tool integration."""
    result = subprocess.run(['amalgkit', '--version'], capture_output=True, text=True)
    assert result.returncode == 0
```

## Environment Variables

Some tests require environment variables for configuration:

```bash
# For NCBI tests
export NCBI_EMAIL="your.email@example.com"

# For database tests
export TEST_DATABASE_URL="sqlite:///output/test.db"

# Run tests with environment
./scripts/run_tests.sh --network
```

## Test Data Management

### Input Fixtures

- **Small test data**: Include directly in `tests/data/`
- **Generated test data**: Create programmatically in `setup_method`
- **Large test data**: Download during test execution (with caching)

### Output Management

- **Consistent paths**: Use `metainformant.core.io.ensure_directory`
- **Deterministic names**: Use consistent, predictable filenames
- **Cleanup policy**: Generally leave outputs for inspection (they're in `output/`)

### Example Test Data Handling

```python
def setup_method(self):
    """Setup test environment with data."""
    self.output_dir = ensure_directory("output/test_analysis")
    
    # Create test data
    self.test_sequences = [
        "ATCGATCGATCG",
        "GCTAGCTAGCTA",
        "TTTTAAAACCCC"
    ]
    
    # Write test FASTA
    self.test_fasta = f"{self.output_dir}/test.fasta"
    with open(self.test_fasta, 'w') as f:
        for i, seq in enumerate(self.test_sequences):
            f.write(f">seq_{i}\n{seq}\n")
```

## Continuous Integration

The test suite is designed for CI/CD environments:

### CI Configuration

```yaml
# Example GitHub Actions
- name: Run fast tests
  run: ./scripts/run_tests.sh --fast

- name: Run network tests  
  run: ./scripts/run_tests.sh --network
  if: env.ENABLE_NETWORK_TESTS == 'true'

- name: Generate coverage report
  run: ./scripts/run_tests.sh --coverage
```

### Test Selection

- **Default CI**: Run fast tests only
- **Nightly builds**: Include slow and network tests
- **Release testing**: Full test suite including integration tests

## Performance Considerations

### Test Execution Speed

- **Fast tests**: < 1 second per test
- **Slow tests**: Marked with `@pytest.mark.slow`
- **Parallel execution**: Tests designed for parallel execution
- **Resource usage**: Tests avoid excessive memory/CPU usage

### Large Dataset Testing

```python
@pytest.mark.slow
def test_large_dataset_processing(self):
    """Test with large synthetic dataset."""
    # Generate large test data efficiently
    large_data = generate_synthetic_data(n_samples=100000)
    
    # Test processing
    result = process_large_data(large_data)
    assert len(result) == len(large_data)
```

## Development Workflow

### Test-Driven Development

1. **Write failing test**: Implement test for new functionality
2. **Implement feature**: Write minimal code to pass test
3. **Refactor**: Improve implementation while maintaining test passage
4. **Add edge cases**: Test error conditions and boundary cases

### Running Tests During Development

```bash
# Run tests for specific module during development
uv run pytest tests/test_core_text.py -v --tb=short

# Watch mode for continuous testing (with external tool)
# uv pip install pytest-watch
ptw tests/test_core_text.py

# Run with immediate failure reporting
uv run pytest -x tests/test_core_text.py
```

Related: [CLI](./cli.md), [Core](./core.md)
