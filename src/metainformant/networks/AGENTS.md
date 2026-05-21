# Agent Directives: networks

**Context**: Network analysis module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `community`, `graph`, `pathway`
- `config/` — exports: `config`
- `interaction/` — exports: `ppi`, `regulatory`
- `regulatory/` — exports: `grn_inference`, `motif_analysis`
- `workflow/` — exports: `workflow`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/networks/](../../../docs/networks/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Gwas module**: [../gwas/AGENTS.md](../gwas/AGENTS.md)
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Ml module**: [../ml/AGENTS.md](../ml/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
