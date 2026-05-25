# Agent Directives: protein

**Context**: Protein sequence and structure analysis module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `database/` — exports: `interpro`, `uniprot`
- `domains/` — exports: `classification`, `detection`
- `function/` — exports: `prediction`
- `sequence/` — exports: `alignment`, `proteomes`, `sequences`
- `structure/` — exports: `alphafold`, `analysis`, `contacts`, `general`, `io`
- `visualization/` — exports: `general`
- `workflow/` — exports: `orchestration`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/protein/](../../../docs/protein/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
