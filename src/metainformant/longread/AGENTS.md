# Agent Directives: longread

**Context**: PacBio and Oxford Nanopore long-read analysis: signal I/O, quality assessment, assembly, methylation calling, haplotype phasing, and structural variant detection.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `modified_bases`, `structural`, `phasing`
- `assembly/` — exports: `overlap`, `consensus`, `hybrid`
- `io/` — exports: `fast5`, `bam`, `formats`
- `methylation/` — exports: `calling`
- `phasing/` — exports: `haplotyping`
- `quality/` — exports: `metrics`, `filtering`
- `utils/` — exports: `batch`, `summary`
- `visualization/` — exports: `plots`
- `workflow/` — exports: `orchestrator`, `orchestrator_core`, `pipeline_stages`, `pipelines`, `reporting`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/longread/](../../../docs/longread/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Structural_Variants module**: [../structural_variants/AGENTS.md](../structural_variants/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
