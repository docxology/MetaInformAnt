# Agent Directives: dna

**Context**: DNA sequence analysis and genomics module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `alignment/` — exports: `distances`, `msa`, `pairwise`
- `annotation/` — exports: `gene_prediction`, `functional`
- `expression/` — exports: `codon`, `transcription`, `translation`
- `external/` — exports: `entrez`, `genomes`, `ncbi`
- `integration/` — exports: `rna`
- `io/` — exports: `fasta`, `fastq`
- `phylogeny/` — exports: `tree`, `tree_analysis`, `tree_construction`
- `population/` — exports: `analysis`, `core`, `visualization`, `visualization_core`, `visualization_stats`
- `sequence/` — exports: `composition`, `consensus`, `core`, `kmer`, `motifs`
- `variation/` — exports: `calling`, `mutations`, `variants`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
