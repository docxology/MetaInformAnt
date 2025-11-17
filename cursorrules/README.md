# METAINFORMANT Modular Cursor Rules

This directory contains module-specific cursor rules for the METAINFORMANT project. Each module has its own `.cursorrules` file with domain-specific patterns, conventions, and guidelines.

## Structure

- **`core.cursorrules`**: Core utilities and shared patterns
- **`dna.cursorrules`**: DNA sequence analysis and genomics
- **`rna.cursorrules`**: RNA transcriptomic analysis
- **`gwas.cursorrules`**: Genome-wide association studies
- **`protein.cursorrules`**: Protein sequence and structure analysis
- **`singlecell.cursorrules`**: Single-cell transcriptomic analysis
- **`networks.cursorrules`**: Biological network analysis
- **`ml.cursorrules`**: Machine learning for biological data
- **`multiomics.cursorrules`**: Multi-omic data integration
- **`information.cursorrules`**: Information-theoretic analysis
- **`life_events.cursorrules`**: Life course event analysis
- **`math.cursorrules`**: Mathematical and theoretical biology
- **`visualization.cursorrules`**: Plotting and visualization
- **`quality.cursorrules`**: Data quality assessment
- **`ontology.cursorrules`**: Functional annotation and ontologies
- **`phenotype.cursorrules`**: Phenotypic trait analysis
- **`epigenome.cursorrules`**: Epigenetic modification analysis
- **`ecology.cursorrules`**: Ecological metadata and community analysis
- **`simulation.cursorrules`**: Synthetic data generation

## Usage

The main `.cursorrules` file in the repository root contains common rules applicable to all modules. Module-specific rules in this directory supplement the main rules with domain-specific patterns.

### When Working on a Module

1. **Start with the main `.cursorrules`** for:
   - Directory structure and path handling
   - Configuration patterns
   - Testing policy (no mocks)
   - Import patterns and code style
   - Documentation guidelines

2. **Then consult the module-specific `.cursorrules`** for:
   - Domain-specific patterns and conventions
   - Module-specific dependencies
   - Output directory structures
   - Environment variable prefixes
   - Integration patterns

### Finding Module Rules

- Module name maps directly to filename: `dna` → `dna.cursorrules`
- All modules follow the same structure for consistency
- Each file includes purpose, dependencies, patterns, and examples

## Common Rules

All modules follow these core principles:
- Write outputs to `output/` by default
- Use `config/` for configuration with env overrides
- **STRICTLY NO MOCKING** in tests (real implementations only - see main `.cursorrules` NO_MOCKING_POLICY)
- **CRITICAL**: Always use `metainformant.core.io` for JSON/CSV/TSV operations - never use direct `import json` or `import csv`
- Use `metainformant.core.paths` utilities for all path operations
- Use `metainformant.core` utilities for I/O, logging, paths, config
- Comprehensive type hints (Python 3.11+)
- Update existing docs, never create root-level docs
- All tests use `tmp_path` fixture and write to `output/` only

## Module-Specific Details

Each module cursorrules file includes:
- **Purpose**: What the module does
- **Dependencies**: Required and optional dependencies
- **Key Submodules**: Main components and their roles
- **Patterns**: Code patterns with examples (including I/O operations using `core.io`)
- **Configuration**: Config structure and environment variables (with proper prefixes)
- **Output Paths**: Where module writes its outputs (matching main `.cursorrules`)
- **Integration**: How it connects with other modules
- **Testing**: Module-specific testing patterns (NO_MOCKING_POLICY enforced)

**Key Improvements Made**:
- All files now emphasize using `metainformant.core.io` instead of direct `json`/`csv` imports
- All files include NO_MOCKING_POLICY in testing sections
- All files document environment variable prefixes consistently
- All files reference main `.cursorrules` for common patterns
- I/O pattern examples added to all modules
- Path handling patterns with `paths.expand_and_resolve()` and `paths.is_within()` added to all modules
- Type hints guidance added to all modules (including `from __future__ import annotations`)
- Consistent section ordering: Patterns → Configuration → Integration → Reference → Testing

For example, see `rna.cursorrules` for detailed amalgkit workflow patterns or `gwas.cursorrules` for statistical analysis conventions.

