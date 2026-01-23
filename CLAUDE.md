# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

METAINFORMANT is a bioinformatics toolkit for multi-omic analysis (genomics, transcriptomics, proteomics, epigenomics, systems biology). Python 3.11+.

**Critical Rules:**
- **NO MOCKING**: Real implementations only. Never return dummy/placeholder data. See `docs/NO_MOCKING_POLICY.md`
- **UV Only**: Use `uv` for all package operations (never `pip`)
- **Output to `output/`**: All execution outputs go to `output/` directory
- **Temp files in `.tmp/`**: Use repository-local `.tmp/` (not system `/tmp`)

## Commands

### Setup & Install

```bash
bash scripts/package/setup.sh              # Full environment setup
source .venv/bin/activate                  # Activate venv
uv pip install -e ".[dev,scientific]"      # Editable install with extras
uv sync --extra dev --extra scientific     # Sync from lock file
```

### Testing

```bash
bash scripts/package/test.sh                          # Fast tests (default)
bash scripts/package/test.sh --mode ultra-fast        # Core only (~5s)
bash scripts/package/test.sh --mode coverage          # Full coverage
bash scripts/package/test.sh --mode network           # Real API tests
bash scripts/package/test.sh --pattern "test_core_*"  # Pattern match
pytest tests/test_dna_comprehensive.py -v             # Single file
pytest tests/test_core_functionality.py::TestCoreIO::test_json_operations -v  # Single test
```

### Code Quality

```bash
bash scripts/package/uv_quality.sh           # All checks (format, lint, typecheck)
bash scripts/package/uv_quality.sh format    # Black formatting
bash scripts/package/uv_quality.sh lint      # Linting
bash scripts/package/uv_quality.sh typecheck # mypy
```

### CLI

```bash
uv run metainformant --help                  # Show all commands
uv run metainformant dna align --input data/sequences.fasta --output output/dna/
uv run metainformant rna run-config --config config/amalgkit/species.yaml
uv run metainformant gwas run --config config/gwas/gwas_template.yaml
```

## Architecture

### Module Structure

```
src/metainformant/
├── core/           # Shared utilities (io, paths, config, logging, cache, parallel)
├── dna/            # DNA: sequences, alignment, phylogeny, population genetics
├── rna/            # RNA-seq: amalgkit integration, workflow orchestration
├── protein/        # Protein: sequences, structures, AlphaFold, UniProt
├── gwas/           # GWAS: association, QC, visualization
├── epigenome/      # Methylation, ChIP-seq, ATAC-seq
├── networks/       # Biological networks, community detection
├── multiomics/     # Multi-omic integration
├── singlecell/     # scRNA-seq analysis
├── visualization/  # 57+ plot types
├── quality/        # QC metrics
├── ml/             # Machine learning pipelines
├── math/           # Population genetics theory, coalescent
├── information/    # Information theory (entropy, MI)
├── ontology/       # GO analysis, semantic similarity
├── phenotype/      # Trait analysis
├── ecology/        # Community diversity
├── simulation/     # Synthetic data generation
└── life_events/    # Event sequence analysis
```

### Data Flow

```
data/ (inputs) → Module Processing → output/ (results)
                       ↑
                 config/ (settings)
```

### Core Utilities

```python
from metainformant.core import io, paths, config
from metainformant.core.logging import get_logger

# File I/O (auto-creates dirs, gzip-aware)
data = io.load_json("config/example.yaml")
io.dump_json(result, "output/module/result.json")

# Path handling with security validation
resolved = paths.expand_and_resolve("~/data/input.txt")
if not paths.is_within(resolved, base_path="/safe/dir"):
    raise ValueError("Path outside safe directory")

# Config with env overrides
cfg = config.load_mapping_from_file("config/workflow.yaml")
cfg = config.apply_env_overrides(cfg, prefix="DOMAIN")

# Logging
logger = get_logger(__name__)
```

### Config Environment Prefixes

- Core: `CORE_`, DNA: `DNA_`, RNA: `AK_`, GWAS: `GWAS_`
- Protein: `PROT_`, Epigenome: `EPI_`, Ontology: `ONT_`
- Others: `{MODULE}_` pattern (e.g., `NET_`, `ML_`, `VIZ_`)

## Testing Requirements

**All tests use real implementations** (NO mocking policy):

```python
def test_real_functionality(tmp_path: Path) -> None:
    output_file = tmp_path / "result.json"
    result = actual_function(output_file)
    assert output_file.exists()

@pytest.mark.network
def test_api_integration() -> None:
    try:
        result = fetch_from_api()
        assert result is not None
    except requests.RequestException as e:
        pytest.skip(f"API unavailable: {e}")
```

**Test Markers:**
- `@pytest.mark.slow` - Long-running
- `@pytest.mark.network` - Real API calls
- `@pytest.mark.external_tool` - CLI tools (amalgkit, muscle)
- `@pytest.mark.integration` - Cross-module

## Code Style

- Python 3.11+, 120 char lines
- Black formatting, mypy type checking
- Import order: `__future__` → stdlib → third-party → local → optional (try/except)
- Always use type hints

## External Dependencies

Optional tools for workflows:
- **amalgkit**: `uv pip install amalgkit`
- **SRA Toolkit**: `apt-get install sra-toolkit`
- **samtools/bcftools/bwa**: `apt-get install samtools bcftools bwa`
- **muscle/clustalo**: `apt-get install muscle clustalo`

## Documentation

- `docs/`: Domain documentation (`docs/<domain>/<topic>.md`)
- `cursorrules/`: Module-specific coding guidelines
- `.cursorrules`: Detailed project rules

**Never create docs in repo root or `output/`** - use `docs/` subdirectories.
