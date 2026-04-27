# Agent Directives: epigenome

**Context**: Epigenome analysis module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `tracks`
- `assays/` — exports: `atacseq`, `chipseq`, `methylation`
- `chromatin_state/` — exports: `state_learning`
- `peak_calling/` — exports: `peak_detection`
- `visualization/` — exports: `visualization`
- `workflow/` — exports: `workflow`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/epigenome/](../../../docs/epigenome/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Protein module**: [../protein/AGENTS.md](../protein/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
