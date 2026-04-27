# Agent Directives: core

**Context**: Core utilities for METAINFORMANT bioinformatics toolkit.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `data/` — exports: `db`, `validation`
- `engine/` — exports: `workflow_manager`
- `execution/` — exports: `discovery`, `parallel`, `workflow`
- `io/` — exports: `atomic`, `cache`, `checksums`, `disk`, `download`
- `output/`
- `ui/` — exports: `tui`
- `utils/` — exports: `config`, `errors`, `hash`, `logging`, `optional_deps`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/core/](../../../docs/core/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **All modules**: [../../../docs/index.md](../../../docs/index.md)
