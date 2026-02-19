# Agent Directives: Rules

## Role

Module-specific Agent Rules for consistent development patterns.

## Contents

Each Rule file contains domain-specific guidelines:

- `core.md` - Core infrastructure patterns
- `dna.md` - DNA sequence analysis patterns
- `ecology.md` - Ecology analysis patterns
- `epigenome.md` - Epigenome analysis patterns
- `gwas.md` - GWAS pipeline patterns
- `information.md` - Information theory patterns
- `life_events.md` - Life events patterns
- `math.md` - Mathematical biology patterns
- `ml.md` - Machine learning patterns
- `multiomics.md` - Multi-omics patterns
- `networks.md` - Network analysis patterns
- `ontology.md` - Ontology patterns
- `phenotype.md` - Phenotype analysis patterns
- `protein.md` - Protein analysis patterns
- `quality.md` - Quality control patterns
- `rna.md` - RNA-seq and amalgkit patterns
- `simulation.md` - Simulation patterns
- `singlecell.md` - Single-cell patterns
- `visualization.md` - Visualization patterns
- `longread.md` - Long-read sequencing (PacBio/Nanopore) patterns
- `metagenomics.md` - Metagenomics (amplicon, shotgun) patterns
- `structural_variants.md` - Structural variant detection patterns
- `spatial.md` - Spatial transcriptomics patterns
- `pharmacogenomics.md` - Clinical pharmacogenomics patterns
- `metabolomics.md` - Metabolomics analysis patterns
- `menu.md` - Interactive menu system patterns

## Usage

These rules are automatically loaded by the Agent when working in the corresponding module directories. They ensure:

- Consistent coding patterns across the project
- Proper use of core utilities
- Adherence to NO MOCKING policy
- Correct output directory usage

## Key Patterns Enforced

All rules enforce:

- Use `metainformant.core.io` for file operations
- Use `metainformant.core.utils.logging` for logging
- Write outputs to `output/` directory
- Use absolute imports from `metainformant`
- Never use mocks in tests
- Use `uv` for package management

## See Also

- `../README.md` - Full documentation of the Agent Rules system
- Global Project Rules (root `.cursorrules`)
