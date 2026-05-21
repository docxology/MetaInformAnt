# Agent Directives: scripts

## Role
Executable scripts and utilities for running METAINFORMANT workflows.

## Directory Structure
Domain-specific scripts organized by module:
- `cloud/` - Cloud deployment and sync utilities
- `core/` - Core infrastructure utilities
- `dna/` - DNA analysis scripts
- `ecology/` - Ecology analysis scripts
- `epigenome/` - Epigenome analysis scripts
- `eqtl/` - eQTL analysis scripts
- `gwas/` - GWAS pipeline scripts (with subdirectories)
- `life_events/` - Life events analysis workflows
- `maintenance/` - Refactoring, import fixes, code generation
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
- `quality/` - Auditing, linting, export checks
- `reorganize/` - Import rewriting and migration
- `rna/` - RNA-seq and amalgkit scripts
- `simulation/` - Simulation scripts for all modules
- `singlecell/` - Single-cell analysis
- `test_examples/` - Example testing, validation, benchmarking
- `validate/` - Project completeness audits
- `visualization/` - Visualization scripts

## Rules and Constraints
- All scripts should be executable directly or via `uv run`
- Scripts should use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is acceptable for CLI glue, subprocess wrappers, and narrow parser internals when tested.
- Output goes to `output/` directory
- Follow the NO MOCKING policy
