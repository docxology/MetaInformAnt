# Scripts

Workflow orchestration scripts and utilities for METAINFORMANT.

## Design Philosophy

Scripts are **thin wrappers** around `src/metainformant/` modules:
- Business logic resides in source modules
- Scripts handle CLI, configuration, and orchestration
- All outputs go to `output/` directory

## Directory Structure

```
scripts/
 cloud/           # Cloud deployment and sync utilities
 core/            # Core utility scripts
 dna/             # DNA analysis workflows
 ecology/         # Ecology analysis scripts
 epigenome/       # Epigenome analysis scripts
 eqtl/            # eQTL analysis scripts
 gwas/            # GWAS pipeline scripts
 life_events/     # Life events analysis workflows
 maintenance/     # Refactoring, import fixes, code generation
 math/            # Mathematical biology scripts
 menu/            # Interactive menu system
 ml/              # Machine learning scripts
 multiomics/      # Multi-omics integration
 networks/        # Network analysis scripts
 ontology/        # Ontology analysis
 package/         # Package management (setup, test, build)
 phenotype/       # Phenotype analysis
 popgen/          # Population genetics scripts
 protein/         # Protein analysis scripts
 quality/         # Auditing, linting, export checks
 reorganize/      # Import rewriting and migration
 rna/             # RNA-seq workflow orchestration
 simulation/      # Simulation scripts for all modules
 singlecell/      # Single-cell analysis
 test_examples/   # Example testing, validation, benchmarking
 validate/        # Project completeness audits
 visualization/   # Batch plotting utilities
```

## Key Scripts

### Package Management (`package/`)

```bash
# Setup environment
bash scripts/package/setup.sh

# Run tests
bash scripts/package/test.sh
bash scripts/package/test.sh --mode coverage
bash scripts/package/test.sh --mode network

# Code quality
bash scripts/package/uv_quality.sh         # All checks
bash scripts/package/uv_quality.sh format  # Black only
bash scripts/package/uv_quality.sh lint    # Linting only

# Build
bash scripts/package/build.sh
bash scripts/package/validate_build.sh
```

### RNA Workflows (`rna/`)

```bash
# List available amalgkit workflow configs
uv run python scripts/rna/run_workflow.py --list-configs

# Check workflow status
uv run python scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Discover new species
uv run python scripts/rna/discover_species.py --species "Apis_mellifera"

# Validate amalgkit species configs
uv run python scripts/rna/validate_configs.py
```

### GWAS Workflows (`gwas/`)

```bash
# Run GWAS pipeline
uv run python scripts/gwas/run_amellifera_gwas.py --config config/gwas/gwas_amellifera.yaml --output output/gwas/amellifera

# Generate all visualizations
uv run python scripts/gwas/visualization/visualizations.py --help
```

### Quality & Auditing (`quality/`)

```bash
# Audit documentation presence
uv run python scripts/quality/audit_docs.py

# Check export completeness
uv run python scripts/quality/check_exports.py

# Audit docstring coverage
uv run python scripts/quality/audit_docstrings.py
```

### Example Testing (`test_examples/`)

```bash
# Validate examples
uv run python scripts/test_examples/validate_examples.py --fast

# Fast validation (pre-commit)
uv run python scripts/test_examples/validate_examples.py --fast

# Benchmark examples
uv run python scripts/test_examples/benchmark_examples.py
```

### Visualization (`visualization/`)

```bash
# Generate batch plots
uv run python scripts/visualization/run_visualization.py --input data/matrix.csv --plot-type heatmap --output output/plots --dry-run
```

## Script Pattern

All scripts follow this pattern:

```python
#!/usr/bin/env python3
"""Script description.

Usage:
    python3 scripts/<domain>/script.py --config config.yaml
"""
import argparse
from pathlib import Path

from metainformant.<domain> import module

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--output", default="output/", help="Output directory")
    args = parser.parse_args()

    # Load config
    config = module.load_config(args.config)

    # Run workflow
    results = module.run_workflow(config)

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    module.save_results(results, output_dir)

if __name__ == "__main__":
    main()
```

## Common Utilities (`_common.sh`)

Shared bash utilities for scripts:

```bash
source scripts/package/_common.sh

# Available functions
setup_environment    # Setup Python/UV environment
check_uv            # Verify UV is installed
sync_dependencies   # Sync from pyproject.toml
print_status        # Colored status messages
```

## Environment Variables

Scripts respect these environment variables:

| Variable | Purpose | Default |
|----------|---------|---------|
| `METAINFORMANT_THREADS` | Parallel threads | 8 |
| `METAINFORMANT_OUTPUT` | Output directory | `output/` |
| `TMPDIR` | Temp directory | `.tmp/` |
| `AK_*` | RNA config overrides | - |
| `GWAS_*` | GWAS config overrides | - |

## Output Structure

All scripts write to `output/`:

```
output/
 amalgkit/<species>/ # RNA workflow outputs
 gwas/<species>/ # GWAS outputs
 visualization/ # Generated plots
 logs/ # Script logs
 ...
```

## Related

- [Configuration](../config/README.md)
- [Source Modules](../src/metainformant/)
- [Testing](../tests/README.md)
