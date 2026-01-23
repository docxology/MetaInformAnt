# Agent Directives: scripts

## Role
Executable scripts and utilities for running METAINFORMANT workflows.

## Directory Structure
Domain-specific scripts organized by module:
- `core/` - Core infrastructure utilities
- `dna/` - DNA analysis scripts
- `ecology/` - Ecology analysis scripts
- `epigenome/` - Epigenome analysis scripts
- `gwas/` - GWAS pipeline scripts (with subdirectories)
- `life_events/` - Life events analysis workflows
- `math/` - Mathematical biology scripts
- `menu/` - Interactive menu system
- `ml/` - Machine learning scripts
- `multiomics/` - Multi-omics integration
- `networks/` - Network analysis scripts
- `ontology/` - Ontology analysis
- `package/` - Package management and build scripts
- `phenotype/` - Phenotype analysis
- `popgen/` - Population genetics scripts
- `protein/` - Protein analysis
- `quality/` - Quality control scripts
- `rna/` - RNA-seq and amalgkit scripts
- `simulation/` - Simulation scripts for all modules
- `singlecell/` - Single-cell analysis
- `test_examples/` - Example testing scripts
- `visualization/` - Visualization scripts

## Root-Level Utilities
- `analyze_example_coverage.py` - Check example coverage
- `benchmark_examples.py` - Benchmark example scripts
- `generate_example.py` - Generate example files
- `validate_examples.py` - Validate example scripts
- `refactor_imports.py` - Import refactoring utility

## Rules and Constraints
- All scripts should be executable directly or via `uv run`
- Scripts should use `metainformant.core` utilities for I/O
- Output goes to `output/` directory
- Follow the NO MOCKING policy
