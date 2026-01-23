# Specification: scripts

## Scope

Executable scripts and utilities for running METAINFORMANT workflows. Contains thin wrapper orchestrators that invoke core library functions, package management tools, and utility scripts. Scripts are organized by domain module.

## Architecture

- **Dependency Level**: Orchestration
- **Component Type**: Executable Scripts
- **Design Pattern**: Thin Wrapper Orchestration

### Directory Structure
```
scripts/
├── package/            # Package management and build scripts
│   ├── setup.sh        # Environment setup
│   ├── test.sh         # Test runner
│   ├── build.sh        # Package builder
│   └── uv_quality.sh   # Code quality checks
├── core/               # Core utility scripts
├── {module}/           # Domain-specific scripts
│   ├── run_workflow.py # Workflow orchestration
│   └── {feature}.py    # Feature-specific scripts
└── *.py                # Root-level utility scripts
```

## Data Structures

### Script Types
- **Workflow Scripts**: Orchestrate multi-step pipelines (run_workflow.py, run_analysis.py)
- **Utility Scripts**: Single-purpose tools (audit_docs.py, check_exports.py)
- **Package Scripts**: Build, test, and quality tools (setup.sh, test.sh, build.sh)
- **Generation Scripts**: Generate code, docs, or configs (generate_example.py, generate_specs.py)

### Module Subdirectories
Each domain has a corresponding scripts subdirectory:
- core/, dna/, rna/, gwas/, protein/, epigenome/
- networks/, multiomics/, singlecell/, visualization/
- quality/, ml/, math/, menu/, ontology/
- phenotype/, ecology/, simulation/, life_events/

### Package Scripts (scripts/package/)
| Script | Purpose |
|--------|---------|
| setup.sh | Install development environment |
| test.sh | Run pytest with various modes |
| build.sh | Build distributable package |
| uv_quality.sh | Run formatters, linters, type checkers |
| uv_docs.sh | Build Sphinx documentation |
| validate_build.sh | Verify build artifacts |

## Interface

### Running Scripts
```bash
# Package management
bash scripts/package/setup.sh        # Setup environment
bash scripts/package/test.sh         # Run tests
bash scripts/package/uv_quality.sh   # Code quality

# Domain workflows
python scripts/rna/run_workflow.py --config config/amalgkit/species.yaml
uv run python scripts/gwas/run_analysis.py --config config/gwas/gwas.yaml

# Utilities
python scripts/generate_example.py --module dna --feature alignment
python scripts/audit_docs.py
```

### Script Conventions
- All scripts executable directly or via `uv run`
- Use `metainformant.core` utilities for I/O and logging
- Write all outputs to `output/` directory
- Support `--help` flag for argument documentation
- Follow NO MOCKING policy (real implementations only)
- Use argparse or click for CLI argument handling
