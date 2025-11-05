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
- No mocks in tests (real implementations only)
- Use `metainformant.core` utilities for I/O, logging, paths
- Comprehensive type hints (Python 3.11+)
- Update existing docs, never create root-level docs

## Module-Specific Details

Each module cursorrules file includes:
- **Purpose**: What the module does
- **Dependencies**: Required and optional dependencies
- **Key Submodules**: Main components and their roles
- **Patterns**: Code patterns with examples
- **Configuration**: Config structure and environment variables
- **Output Paths**: Where module writes its outputs
- **Integration**: How it connects with other modules
- **Testing**: Module-specific testing patterns

For example, see `rna.cursorrules` for detailed amalgkit workflow patterns or `gwas.cursorrules` for statistical analysis conventions.

