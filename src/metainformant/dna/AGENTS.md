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
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/dna/](../../../docs/dna/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Rna module**: [../rna/AGENTS.md](../rna/AGENTS.md)
- **Gwas module**: [../gwas/AGENTS.md](../gwas/AGENTS.md)
- **Quality module**: [../quality/AGENTS.md](../quality/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
