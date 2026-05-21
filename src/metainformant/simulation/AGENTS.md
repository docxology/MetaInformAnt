# Agent Directives: simulation

**Context**: Evolutionary and Population Genetics Simulation module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `benchmark/` — exports: `generators`
- `methylation/`
- `models/` — exports: `agents`, `popgen`, `rna`, `sequences`
- `visualization/` — exports: `visualization`
- `workflow/` — exports: `workflow`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/simulation/](../../../docs/simulation/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Math module**: [../math/AGENTS.md](../math/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
