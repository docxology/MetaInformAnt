# Agent Directives: structural_variants

**Context**: Structural variant analysis: CNV detection via circular binary segmentation, SV calling from split/discordant reads, annotation, filtering, population genotyping, and visualization.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `annotation/` — exports: `overlap`, `functional_impact`
- `detection/` — exports: `cnv`, `sv_calling`, `breakpoints`
- `filtering/` — exports: `quality_filter`, `merge`
- `population/` — exports: `sv_population`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/structural_variants/](../../../docs/structural_variants/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Longread module**: [../longread/AGENTS.md](../longread/AGENTS.md)
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
