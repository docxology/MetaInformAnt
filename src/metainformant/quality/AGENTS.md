# Agent Directives: quality

**Context**: Quality control analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `contamination`, `metrics`
- `batch/`
- `io/` — exports: `fastq`
- `reporting/` — exports: `multiqc_integration`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
