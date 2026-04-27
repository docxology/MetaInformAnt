# Agent Directives: singlecell

**Context**: Single-cell analysis module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `clustering`, `dimensionality`, `nonlinear_methods`, `pca_methods`, `trajectory`
- `celltyping/` — exports: `annotation`
- `data/` — exports: `integration`, `preprocessing`
- `differential/` — exports: `expression`
- `doublet/`
- `velocity/` — exports: `rna_velocity`
- `visualization/` — exports: `visualization`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/singlecell/](../../../docs/singlecell/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Rna module**: [../rna/AGENTS.md](../rna/AGENTS.md)
- **Spatial module**: [../spatial/AGENTS.md](../spatial/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
