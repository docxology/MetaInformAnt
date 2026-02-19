# Agent Directives: singlecell

**Context**: Single-cell analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `clustering`, `dimensionality`, `nonlinear_methods`, `pca_methods`, `trajectory`
- `celltyping/` — exports: `annotation`
- `data/` — exports: `integration`, `preprocessing`
- `differential/` — exports: `expression`
- `doublet/`
- `velocity/` — exports: `rna_velocity`
- `visualization/` — exports: `visualization`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
