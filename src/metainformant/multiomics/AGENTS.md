# Agent Directives: multiomics

**Context**: Multi-omics integration module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `integration`
- `methods/` — exports: `factorization`, `clustering`
- `pathways/` — exports: `enrichment`
- `survival/` — exports: `analysis`
- `visualization/` — exports: `visualization`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/multiomics/](../../../docs/multiomics/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Rna module**: [../rna/AGENTS.md](../rna/AGENTS.md)
- **Metabolomics module**: [../metabolomics/AGENTS.md](../metabolomics/AGENTS.md)
- **Protein module**: [../protein/AGENTS.md](../protein/AGENTS.md)
