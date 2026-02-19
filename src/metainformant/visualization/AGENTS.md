# Agent Directives: visualization

**Context**: Visualization and plotting utilities module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `dimred`, `information`, `quality`, `quality_assessment`, `quality_omics`
- `config/` — exports: `palettes`, `themes`
- `dashboards/` — exports: `composite`, `interactive`
- `genomics/` — exports: `expression`, `genomics`, `networks`, `trees`
- `interactive_dashboards/` — exports: `dashboards`
- `plots/` — exports: `animations`, `basic`, `general`, `multidim`, `specialized`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
