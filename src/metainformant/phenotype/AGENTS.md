# Agent Directives: phenotype

**Context**: Phenotype module for MetaInformAnt.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `life_course`
- `behavior/` — exports: `ethogram`, `sequence`
- `chemical/` — exports: `compound`, `profile`
- `data/` — exports: `antwiki`, `scraper`
- `electronic/` — exports: `tracking`
- `gwas_integration/` — exports: `phewas`
- `integration/` — exports: `cross_omic`
- `morphological/` — exports: `measurement`, `profile`
- `sonic/` — exports: `signal`
- `visualization/` — exports: `visualization`
- `workflow/` — exports: `pipeline`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/phenotype/](../../../docs/phenotype/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Gwas module**: [../gwas/AGENTS.md](../gwas/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
