# Agent Directives: structural_variants



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `annotation/` — exports: `overlap`, `functional_impact`
- `detection/` — exports: `cnv`, `sv_calling`, `breakpoints`
- `filtering/` — exports: `quality_filter`, `merge`
- `population/` — exports: `sv_population`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
