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
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
