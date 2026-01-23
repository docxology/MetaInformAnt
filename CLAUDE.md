# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

METAINFORMANT is a comprehensive bioinformatics toolkit for multi-omic analysis across genomics, transcriptomics, proteomics, epigenomics, and systems biology. Built with Python 3.11+, it provides production-ready tools for DNA/RNA analysis, GWAS, population genetics, network analysis, machine learning, and visualization.

**Key Characteristics:**
- **NO MOCKING POLICY**: All code uses real implementations, real API calls, and real algorithms. Never return dummy/placeholder data.
- **UV Package Manager**: All Python package operations use `uv` exclusively (never `pip`)
- **Output Directory Policy**: All program execution outputs go to `output/` directory
- **Multi-domain Integration**: 19 integrated modules with cross-module workflows

## Build, Test & Development Commands

### Environment Setup

```bash
# Setup entire environment (installs dev + scientific dependencies)
bash scripts/package/setup.sh

# Activate virtual environment
source .venv/bin/activate  # Standard filesystem
source /tmp/metainformant_venv/bin/activate  # FAT filesystem (auto-detected)

# Install package in editable mode
uv pip install -e .

# Install with specific extras
uv pip install -e ".[dev,scientific,ml,networks]"

# Sync all dependencies from lock file
uv sync --extra dev --extra scientific
```

### Testing

```bash
# Run fast tests (default, ~15s)
bash scripts/package/test.sh

# Run ultra-fast core tests (~5s)
bash scripts/package/test.sh --mode ultra-fast

# Run full test suite with coverage (~2-5min)
bash scripts/package/test.sh --mode coverage

# Run specific test pattern
bash scripts/package/test.sh --pattern "test_core_*"

# Run integration tests
bash scripts/package/test.sh --mode integration

# Run network tests (real API calls)
bash scripts/package/test.sh --mode network

# Run tests with parallel execution
bash scripts/package/test.sh --mode parallel

# Run single test file
pytest tests/test_dna_comprehensive.py -v

# Run single test function
pytest tests/test_core_functionality.py::TestCoreIO::test_json_operations -v
```

### Code Quality

```bash
# Run all quality checks (format, lint, typecheck)
bash scripts/package/uv_quality.sh

# Format code only
bash scripts/package/uv_quality.sh format

# Lint code only
bash scripts/package/uv_quality.sh lint

# Type check only
bash scripts/package/uv_quality.sh typecheck
```

### Building & Verification

```bash
# Build package
bash scripts/package/build.sh

# Validate build artifacts
bash scripts/package/validate_build.sh

# Verify installation
bash scripts/package/verify.sh
```

## High-Level Architecture

### Core Philosophy

**Real Implementations Only**: METAINFORMANT enforces a strict NO MOCKING policy:
- All functions perform real computations or make genuine API calls
- No dummy data, placeholders, or stubbed methods in source code or tests
- When external dependencies unavailable, raise errors or skip gracefully
- Tests use real implementations with actual file I/O, API calls, and algorithms

### Module Organization

The codebase is organized into 19 domain-specific modules with clear separation of concerns:

**Core Infrastructure** (`src/metainformant/core/`):
- `io.py`: File I/O operations (JSON/JSONL/CSV/TSV, gzip-aware)
- `paths.py`: Path handling, validation, and security (prevents directory traversal)
- `config.py`: Configuration loading with environment variable overrides
- `logging.py`: Centralized logging infrastructure
- `workflow.py`: Workflow orchestration framework
- `cache.py`: JSON-based caching system
- `parallel.py`: Parallel processing utilities
- `validation.py`: Input validation and type checking

**Data Flow Pattern**:
```
Raw Data (data/) → Module Processing → Core Utilities → Output (output/)
                                     ↓
                              Configuration (config/)
```

**Cross-Module Integration**:
- DNA ↔ GWAS: Variant calling and population genetics
- RNA ↔ Single-Cell: Expression data for scRNA-seq
- Protein ↔ Networks: PPI network construction
- Epigenome ↔ Networks: Chromatin interaction networks
- All Modules → Quality: Quality control integration
- All Modules → Visualization: Plotting and visualization

### Configuration System

**Environment Variable Prefixes** (for overriding configs):
- Core: `CORE_` (e.g., `CORE_THREADS=8`)
- DNA: `DNA_` (e.g., `DNA_WORK_DIR=output/dna`)
- RNA: `AK_` (e.g., `AK_THREADS=16`)
- GWAS: `GWAS_` (e.g., `GWAS_WORK_DIR=output/gwas`)
- Protein: `PROT_`, Epigenome: `EPI_`, Ontology: `ONT_`
- Other modules follow `{MODULE}_` pattern

**Config Loading Pattern**:
```python
from metainformant.core.config import load_mapping_from_file, apply_env_overrides

# Load with env overrides
config = load_mapping_from_file("config/domain/workflow.yaml")
config = apply_env_overrides(config, prefix="DOMAIN")
```

### Output Directory Structure

**CRITICAL**: All program execution outputs must go to `output/` directory:
- `output/core/` - Core utility outputs (cache, logs, workflows)
- `output/dna/` - DNA analysis results (phylogeny, variants, alignments)
- `output/rna/` - RNA-seq outputs (amalgkit workflows, quantification)
- `output/gwas/` - GWAS results (association, QC, manhattan plots)
- `output/visualization/` - All generated plots and figures
- And 15+ other domain-specific subdirectories

**Never create documentation, reports, test scripts, or planning documents in `output/`** - these belong in `docs/`, `tests/`, or `scripts/` respectively.

### Temporary File Management

**CRITICAL for External Drives**: System `/tmp` may be full (small tmpfs). Always use repository-local temp:

```python
import os, tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TEMP_DIR = REPO_ROOT / ".tmp" / "python"
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# Override temp directory
os.environ["TMPDIR"] = str(TEMP_DIR)
tempfile.tempdir = str(TEMP_DIR)
```

For shell scripts:
```bash
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
export TMPDIR="$REPO_ROOT/.tmp/bash"
mkdir -p "$TMPDIR"
```

### Package Management with UV

**CRITICAL**: Always use `uv` for all Python package operations:

```bash
# Creating environments
uv venv                          # Create virtual environment
uv venv /tmp/metainformant_venv  # Create in /tmp (FAT filesystems)

# Installing packages
uv pip install -e .              # Editable install
uv pip install package_name      # Install package

# Running commands
uv run pytest                    # Run in project env
uv run python script.py          # Execute script

# Managing dependencies
uv add package_name              # Add dependency
uv remove package_name           # Remove dependency
uv sync                          # Sync from lock file
```

**NEVER use `pip`, `pip3`, `python -m pip` directly** - the project requires `uv` for consistency.

### Testing Architecture

Tests are organized by module with comprehensive integration testing:

**Test Markers**:
- `@pytest.mark.slow` - Long-running tests
- `@pytest.mark.network` - Tests requiring real API calls
- `@pytest.mark.external_tool` - Tests requiring CLI tools (amalgkit, muscle, etc.)
- `@pytest.mark.integration` - Cross-module integration tests

**Test Pattern** (enforcing NO MOCKING):
```python
def test_real_functionality(tmp_path: Path) -> None:
    """Test with real implementation."""
    output_file = tmp_path / "result.json"

    # Test real API call or algorithm
    result = actual_function(output_file)

    # Verify real output
    assert output_file.exists()
    data = io.load_json(output_file)
    assert data["key"] == expected_value

@pytest.mark.network
def test_api_integration() -> None:
    """Test real API - skip if offline."""
    try:
        result = fetch_from_api()
        assert result is not None
    except requests.RequestException as e:
        pytest.skip(f"API unavailable: {e}")
```

### Workflow Execution Patterns

**Python-based workflows** (recommended for programmatic control):
```python
from metainformant.rna.workflow import load_workflow_config, execute_workflow

config = load_workflow_config("config/amalgkit/species.yaml")
results = execute_workflow(config)
```

**CLI-based workflows** (recommended for end-to-end execution):
```bash
# RNA-seq workflow
python3 scripts/rna/run_workflow.py --config config/amalgkit/species.yaml

# GWAS workflow
uv run metainformant gwas run --config config/gwas/gwas_template.yaml

# DNA analysis
uv run metainformant dna align --input data/sequences.fasta --output output/dna/alignment
```

## Important Development Patterns

### Path Handling & Security

```python
from metainformant.core import paths

# Always expand and resolve paths
resolved = paths.expand_and_resolve("~/data/input.txt")

# Validate path containment (prevent directory traversal)
if not paths.is_within(resolved, base_path="/safe/directory"):
    raise ValueError("Path outside safe directory")
```

### File I/O Operations

```python
from metainformant.core import io

# Reading
data = io.load_json("config/example.yaml")
df = io.load_csv("data/table.csv")
records = list(io.read_jsonl("data/records.jsonl"))

# Writing (auto-creates parent directories)
io.dump_json(data, "output/module/result.json")
io.write_csv(df, "output/module/table.csv")

# Gzip-aware operations
with io.open_text_auto("data/file.txt.gz") as f:
    content = f.read()
```

### Logging Pattern

```python
from metainformant.core.logging import get_logger

logger = get_logger(__name__)

def process_data():
    logger.info("Starting operation")
    try:
        result = expensive_computation()
        logger.info(f"Completed: {result}")
        return result
    except Exception as e:
        logger.error(f"Operation failed: {e}", exc_info=True)
        raise
```

### Error Handling (No Dummy Data Fallbacks)

```python
# ✅ CORRECT: Raise error when dependency unavailable
def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]:
    try:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        logger.error(f"Failed to fetch UniProt record: {e}")
        raise  # Propagate error

# ❌ INCORRECT: Never return dummy data
def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]:
    try:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
        return response.json()
    except Exception:
        return {'accession': uniprot_id, 'sequence': 'DUMMY'}  # PROHIBITED
```

### Optional Dependencies

```python
# Lazy import pattern for optional dependencies
def requires_networkx():
    try:
        import networkx as nx
        return nx
    except ImportError:
        raise ImportError(
            "networkx required for this operation. "
            "Install with: uv pip install networkx"
        )

def graph_operation():
    nx = requires_networkx()
    # Use nx for graph operations
```

## Module-Specific Notes

### RNA Module (Amalgkit Integration)
- External CLI tool: `amalgkit` (install via `uv pip install amalgkit`)
- Workflow orchestration in `src/metainformant/rna/workflow.py`
- Config prefix: `AK_`
- Output: `output/amalgkit/<species>/<step>/`

### GWAS Module
- Requires external tools: bcftools, samtools, bwa (optional)
- Real VCF parsing (no dummy genotype matrices)
- Config prefix: `GWAS_`
- Comprehensive visualization suite in `gwas/visualization_*.py`

### DNA Module
- Population genetics: real Tajima's D, nucleotide diversity calculations
- Phylogenetic analysis with real tree construction
- Alignment algorithms: Needleman-Wunsch, Smith-Waterman implementations

### Visualization Module
- 14 specialized plotting modules with 57+ unique plot types
- All functions return real matplotlib Figure objects (never None)
- Output to `output/visualization/` with informative filenames

## Code Style & Standards

- **Python Version**: 3.11+ minimum
- **Line Length**: 120 characters
- **Formatter**: Black (run via `uv run black`)
- **Linter**: Flake8 (run via `uv run flake8`)
- **Type Checker**: mypy (run via `uv run mypy`)
- **Import Order**: `__future__` → stdlib → third-party → local → optional (with try/except)

## Critical Policies

1. **NO MOCKING**: Never use mock/fake/stub/placeholder implementations. See `docs/NO_MOCKING_POLICY.md`
2. **UV Only**: Never use `pip` - always use `uv` for package management
3. **Output Directory**: All execution outputs to `output/` (never docs/reports in `output/`)
4. **Temp Files**: Use `.tmp/` in repo root (not system `/tmp`)
5. **Real APIs**: All network tests make real API calls or skip gracefully
6. **Path Security**: Always validate paths with `paths.is_within()`

## Common Workflows

```bash
# Setup development environment
bash scripts/package/setup.sh

# Run fast tests before committing
bash scripts/package/test.sh --mode ultra-fast

# Check code quality
bash scripts/package/uv_quality.sh

# Run comprehensive workflow demo
python scripts/core/run_demo.py

# Generate documentation
bash scripts/package/uv_docs.sh
```

## External Tool Dependencies

Optional tools for specific workflows:
- **amalgkit**: RNA-seq workflow orchestration (`uv pip install amalgkit`)
- **SRA Toolkit**: Sequencing data download (`apt-get install sra-toolkit`)
- **samtools/bcftools/bwa**: GWAS variant calling (`apt-get install samtools bcftools bwa`)
- **muscle/clustalo**: Multiple sequence alignment (`apt-get install muscle clustalo`)

## Documentation Structure

- `README.md`: Project overview and quick start
- `QUICKSTART.md`: Fast setup commands
- `docs/`: Domain-organized documentation (`docs/<domain>/<topic>.md`)
- `cursorrules/`: Module-specific coding guidelines
- Module READMEs: `src/metainformant/<module>/README.md`

**Never create new documentation in repository root** - update existing docs or add to `docs/` subdirectories.
