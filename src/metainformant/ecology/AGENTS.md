# Agent Directives: ecology

**Context**: Ecology and community analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `community`, `functional`, `indicators`, `macroecology`, `ordination`
- `phylogenetic/` — exports: `diversity`
- `traits/`
- `visualization/` — exports: `visualization`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
