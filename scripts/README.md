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
# Run amalgkit workflow
python3 scripts/rna/run_workflow.py --config config/amalgkit/species.yaml

# Check workflow status
python3 scripts/rna/run_workflow.py --config config/amalgkit/species.yaml --status

# Discover new species
python3 scripts/rna/discover_species.py --species "Apis_mellifera"

# Recover missing samples
python3 scripts/rna/recover_missing_parallel.py --config config/amalgkit/species.yaml
```

### GWAS Workflows (`gwas/`)

```bash
# Run GWAS pipeline
python3 scripts/gwas/run_genome_scale_gwas.py --config config/gwas/species.yaml

# Generate all visualizations
python3 scripts/gwas/generate_plots.py --results output/gwas/results.tsv
```

### Quality & Auditing (`quality/`)

```bash
# Audit documentation presence
python3 scripts/quality/audit_docs.py

# Check export completeness
python3 scripts/quality/check_exports.py

# Audit docstring coverage
python3 scripts/quality/audit_docstrings.py
```

### Example Testing (`test_examples/`)

```bash
# Run all examples
python3 scripts/test_examples/test_examples.py --continue-on-error

# Fast validation (pre-commit)
python3 scripts/test_examples/validate_examples.py --fast

# Benchmark examples
python3 scripts/test_examples/benchmark_examples.py
```

### Visualization (`visualization/`)

```bash
# Generate batch plots
python3 scripts/visualization/batch_plots.py --input data/ --output output/plots/
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
