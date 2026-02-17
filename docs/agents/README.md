# METAINFORMANT Agent Rules

This directory contains module-specific agent rules for the METAINFORMANT project. These files serve as the authoritative source for agent behaviors, coding patterns, and conventions.

## Structure

- **`rules/core.md`**: Core utilities and shared patterns
- **`rules/dna.md`**: DNA sequence analysis and genomics
- **`rules/rna.md`**: RNA transcriptomic analysis
- **`rules/gwas.md`**: Genome-wide association studies
- **`rules/protein.md`**: Protein sequence and structure analysis
- **`rules/singlecell.md`**: Single-cell transcriptomic analysis
- **`rules/networks.md`**: Biological network analysis
- **`rules/ml.md`**: Machine learning for biological data
- **`rules/multiomics.md`**: Multi-omic data integration
- **`rules/information.md`**: Information-theoretic analysis
- **`rules/life_events.md`**: Life course event analysis
- **`rules/math.md`**: Mathematical and theoretical biology
- **`rules/visualization.md`**: Plotting and visualization
- **`rules/quality.md`**: Data quality assessment
- **`rules/ontology.md`**: Functional annotation and ontologies
- **`rules/phenotype.md`**: Phenotypic trait analysis
- **`rules/epigenome.md`**: Epigenetic modification analysis
- **`rules/ecology.md`**: Ecological metadata and community analysis
- **`rules/simulation.md`**: Synthetic data generation
- **`rules/longread.md`**: Long-read sequencing (PacBio/Nanopore)
- **`rules/metagenomics.md`**: Metagenomic analysis (amplicon, shotgun)
- **`rules/structural_variants.md`**: CNV/SV detection and annotation
- **`rules/spatial.md`**: Spatial transcriptomics (Visium, MERFISH, Xenium)
- **`rules/pharmacogenomics.md`**: Clinical pharmacogenomics
- **`rules/menu.md`**: Interactive menu and discovery system

## Usage

These rules provide the "Brain" for our agents. When working on a specific module, consult the corresponding rule file in `rules/`.

### When Working on a Module

1. **Start with the Core Rules (`rules/core.md`)** for:
   - Directory structure and path handling
   - Configuration patterns
   - Testing policy (no mocks)
   - Import patterns and code style
   - Documentation guidelines

2. **Then consult the module-specific rule file** for:
   - Domain-specific patterns and conventions
   - Module-specific dependencies
   - Output directory structures
   - Environment variable prefixes
   - Integration patterns

### Finding Module Rules

- Module name maps directly to filename: `dna` → `rules/dna.md`
- All modules follow the same structure for consistency
- Each file includes purpose, dependencies, patterns, and examples

## Common Rules

All modules follow these core principles:

- **Package Management**: **ALWAYS use `uv`** for all Python package management - `uv venv`, `uv pip install`, `uv run`, `uv sync`, `uv add`, `uv remove`. Never use `pip` directly.
- Write outputs to `output/` by default
- Use `config/` for configuration with env overrides
- **STRICTLY NO MOCKING** in tests (real implementations only - see `rules/core.md` NO_MOCKING_POLICY)
- **CRITICAL**: Always use `metainformant.core.io` for JSON/CSV/TSV operations - never use direct `import json` or `import csv`
- Use `metainformant.core.paths` utilities for all path operations
- Use `metainformant.core` utilities for I/O, logging, paths, config
- Comprehensive type hints (Python 3.11+)
- Update existing docs, never create root-level docs
- All tests use `tmp_path` fixture and write to `output/` only

## Module-Specific Details

Each rule file includes:

- **Purpose**: What the module does
- **Dependencies**: Required and optional dependencies
- **Key Submodules**: Main components and their roles
- **Patterns**: Code patterns with examples (including I/O operations using `core.io`)
- **Configuration**: Config structure and environment variables (with proper prefixes)
- **Output Paths**: Where module writes its outputs
- **Integration**: How it connects with other modules
- **Testing**: Module-specific testing patterns (NO_MOCKING_POLICY enforced)
