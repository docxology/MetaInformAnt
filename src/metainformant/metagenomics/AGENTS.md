# Agent Directives: metagenomics



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `amplicon/` — exports: `otu_clustering`, `asv_denoising`, `taxonomy`
- `comparative/` — exports: `differential_abundance`
- `diversity/` — exports: `metrics`
- `functional/` — exports: `annotation`, `pathways`
- `shotgun/` — exports: `assembly`, `binning`, `profiling`
- `visualization/` — exports: `plots`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
