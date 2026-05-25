# Agent Directives: ontology

**Context**: Gene ontology and functional annotation module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `core/` — exports: `go`, `obo`, `types`
- `pathway_enrichment/` — exports: `enrichment`
- `query/` — exports: `query`, `serialize`
- `visualization/` — exports: `visualization`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/ontology/](../../../docs/ontology/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Epigenome module**: [../epigenome/AGENTS.md](../epigenome/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
