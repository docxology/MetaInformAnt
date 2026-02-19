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
