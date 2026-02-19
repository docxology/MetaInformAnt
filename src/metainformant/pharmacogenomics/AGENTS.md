# Agent Directives: pharmacogenomics



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `alleles/` — exports: `star_allele`, `diplotype`, `phenotype`
- `annotations/` — exports: `cpic`, `pharmgkb`, `drug_labels`
- `clinical/` — exports: `pathogenicity`, `drug_interaction`, `reporting`
- `interaction/` — exports: `drug_interactions`
- `metabolism/` — exports: `metabolizer_status`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
