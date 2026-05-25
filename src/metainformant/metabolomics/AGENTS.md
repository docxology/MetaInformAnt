# Agent Directives: metabolomics

**Context**: Metabolomics analysis: mass spectrometry data processing, metabolite identification, pathway mapping, and metabolite-gene integration.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `identification`
- `io/` — exports: `formats`
- `pathways/` — exports: `enrichment`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/metabolomics/](../../../docs/metabolomics/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Multiomics module**: [../multiomics/AGENTS.md](../multiomics/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
