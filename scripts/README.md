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
├── package/        # Package management (setup, test, build)
├── rna/            # RNA-seq workflow orchestration
├── gwas/           # GWAS pipeline scripts
├── dna/            # DNA analysis workflows
├── protein/        # Protein analysis scripts
├── visualization/  # Batch plotting utilities
├── core/           # Core utility scripts
├── menu/           # Interactive menu system
├── test_examples/  # Example validation
└── <domain>/       # Domain-specific scripts
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
├── amalgkit/<species>/    # RNA workflow outputs
├── gwas/<species>/        # GWAS outputs
├── visualization/         # Generated plots
├── logs/                  # Script logs
└── ...
```

## Related

- [Configuration](../config/README.md)
- [Source Modules](../src/metainformant/)
- [Testing](../tests/README.md)
