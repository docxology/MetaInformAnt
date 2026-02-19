# Agent Directives: protein

**Context**: Protein sequence and structure analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `database/` — exports: `interpro`, `uniprot`
- `domains/` — exports: `classification`, `detection`
- `function/` — exports: `prediction`
- `sequence/` — exports: `alignment`, `proteomes`, `sequences`
- `structure/` — exports: `alphafold`, `analysis`, `contacts`, `general`, `io`
- `visualization/` — exports: `general`
- `workflow/` — exports: `orchestration`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
