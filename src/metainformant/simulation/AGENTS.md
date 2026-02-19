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
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
