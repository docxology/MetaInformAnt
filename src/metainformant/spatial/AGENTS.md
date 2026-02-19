# Agent Directives: spatial



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `autocorrelation`, `clustering`, `deconvolution`, `neighborhood`
- `communication/` — exports: `cell_communication`
- `deconvolution/` — exports: `spatial_deconvolution`
- `integration/` — exports: `scrna_mapping`
- `io/` — exports: `merfish`, `visium`, `xenium`
- `niche/`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
