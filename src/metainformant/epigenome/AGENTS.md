# Agent Directives: epigenome

**Context**: Epigenome analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `tracks`
- `assays/` — exports: `atacseq`, `chipseq`, `methylation`
- `chromatin_state/` — exports: `state_learning`
- `peak_calling/` — exports: `peak_detection`
- `visualization/` — exports: `visualization`
- `workflow/` — exports: `workflow`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
