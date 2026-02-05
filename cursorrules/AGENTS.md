# Agent Directives: cursorrules

## Role
Module-specific Cursor AI rules for consistent development patterns.

## Contents
Each `.cursorrules` file contains domain-specific guidelines:
- `core.cursorrules` - Core infrastructure patterns
- `dna.cursorrules` - DNA sequence analysis patterns
- `ecology.cursorrules` - Ecology analysis patterns
- `epigenome.cursorrules` - Epigenome analysis patterns
- `gwas.cursorrules` - GWAS pipeline patterns
- `information.cursorrules` - Information theory patterns
- `life_events.cursorrules` - Life events patterns
- `math.cursorrules` - Mathematical biology patterns
- `ml.cursorrules` - Machine learning patterns
- `multiomics.cursorrules` - Multi-omics patterns
- `networks.cursorrules` - Network analysis patterns
- `ontology.cursorrules` - Ontology patterns
- `phenotype.cursorrules` - Phenotype analysis patterns
- `protein.cursorrules` - Protein analysis patterns
- `quality.cursorrules` - Quality control patterns
- `rna.cursorrules` - RNA-seq and amalgkit patterns
- `simulation.cursorrules` - Simulation patterns
- `singlecell.cursorrules` - Single-cell patterns
- `visualization.cursorrules` - Visualization patterns
- `longread.cursorrules` - Long-read sequencing (PacBio/Nanopore) patterns
- `metagenomics.cursorrules` - Metagenomics (amplicon, shotgun) patterns
- `structural_variants.cursorrules` - Structural variant detection patterns
- `spatial.cursorrules` - Spatial transcriptomics patterns
- `pharmacogenomics.cursorrules` - Clinical pharmacogenomics patterns
- `menu.cursorrules` - Interactive menu system patterns

## Usage
These rules are automatically loaded by Cursor AI when working in the corresponding module directories. They ensure:
- Consistent coding patterns across the project
- Proper use of core utilities
- Adherence to NO MOCKING policy
- Correct output directory usage

## Key Patterns Enforced
All cursorrules enforce:
- Use `metainformant.core.io` for file operations
- Use `metainformant.core.logging` for logging
- Write outputs to `output/` directory
- Use absolute imports from `metainformant`
- Never use mocks in tests
- Use `uv` for package management

## See Also
- `README.md` - Full documentation of cursorrules system
- Root `.cursorrules` - Global project rules
