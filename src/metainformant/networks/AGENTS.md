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
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
