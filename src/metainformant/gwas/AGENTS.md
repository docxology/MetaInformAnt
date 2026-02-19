# Agent Directives: gwas

**Context**: Genome-Wide Association Studies (GWAS) module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `annotation`, `association`, `calling`, `correction`, `heritability`
- `data/` — exports: `config`, `download`, `genome`, `metadata`, `sra_download`
- `finemapping/` — exports: `colocalization`, `credible_sets`, `eqtl`
- `heritability/` — exports: `estimation`
- `visualization/` — exports: `config`, `eqtl_visualization`, `general`, `genomic`, `interactive`
- `workflow/` — exports: `workflow`, `workflow_config`, `workflow_execution`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
