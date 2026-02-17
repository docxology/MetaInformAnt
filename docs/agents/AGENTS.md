# Agent Directives: Rules

## Role

Module-specific Agent Rules for consistent development patterns.

## Contents

Each Rule file contains domain-specific guidelines:

- `rules/core.md` - Core infrastructure patterns
- `rules/dna.md` - DNA sequence analysis patterns
- `rules/ecology.md` - Ecology analysis patterns
- `rules/epigenome.md` - Epigenome analysis patterns
- `rules/gwas.md` - GWAS pipeline patterns
- `rules/information.md` - Information theory patterns
- `rules/life_events.md` - Life events patterns
- `rules/math.md` - Mathematical biology patterns
- `rules/ml.md` - Machine learning patterns
- `rules/multiomics.md` - Multi-omics patterns
- `rules/networks.md` - Network analysis patterns
- `rules/ontology.md` - Ontology patterns
- `rules/phenotype.md` - Phenotype analysis patterns
- `rules/protein.md` - Protein analysis patterns
- `rules/quality.md` - Quality control patterns
- `rules/rna.md` - RNA-seq and amalgkit patterns
- `rules/simulation.md` - Simulation patterns
- `rules/singlecell.md` - Single-cell patterns
- `rules/visualization.md` - Visualization patterns
- `rules/longread.md` - Long-read sequencing (PacBio/Nanopore) patterns
- `rules/metagenomics.md` - Metagenomics (amplicon, shotgun) patterns
- `rules/structural_variants.md` - Structural variant detection patterns
- `rules/spatial.md` - Spatial transcriptomics patterns
- `rules/pharmacogenomics.md` - Clinical pharmacogenomics patterns
- `rules/menu.md` - Interactive menu system patterns

## Usage

These rules are automatically loaded by the Agent when working in the corresponding module directories. They ensure:

- Consistent coding patterns across the project
- Proper use of core utilities
- Adherence to NO MOCKING policy
- Correct output directory usage

## Key Patterns Enforced

All rules enforce:

- Use `metainformant.core.io` for file operations
- Use `metainformant.core.logging` for logging
- Write outputs to `output/` directory
- Use absolute imports from `metainformant`
- Never use mocks in tests
- Use `uv` for package management

## See Also

- `README.md` - Full documentation of the Agent Rules system
- Global Project Rules (root `.cursorrules`)
