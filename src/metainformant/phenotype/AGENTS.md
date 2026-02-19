# Agent Directives: phenotype

**Context**: Phenotype module for MetaInformAnt.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `life_course`
- `behavior/` — exports: `ethogram`, `sequence`
- `chemical/` — exports: `compound`, `profile`
- `data/` — exports: `antwiki`, `scraper`
- `electronic/` — exports: `tracking`
- `gwas_integration/` — exports: `phewas`
- `integration/` — exports: `cross_omic`
- `morphological/` — exports: `measurement`, `profile`
- `sonic/` — exports: `signal`
- `visualization/` — exports: `visualization`
- `workflow/` — exports: `pipeline`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
