# Agent Directives: metabolomics



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `identification`
- `io/` — exports: `formats`
- `pathways/` — exports: `enrichment`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
