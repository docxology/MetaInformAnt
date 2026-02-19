# Agent Directives: longread



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `modified_bases`, `structural`, `phasing`
- `assembly/` — exports: `overlap`, `consensus`, `hybrid`
- `io/` — exports: `fast5`, `bam`, `formats`
- `methylation/` — exports: `calling`
- `phasing/` — exports: `haplotyping`
- `quality/` — exports: `metrics`, `filtering`
- `utils/` — exports: `batch`, `summary`
- `visualization/` — exports: `plots`
- `workflow/` — exports: `orchestrator`, `orchestrator_core`, `pipeline_stages`, `pipelines`, `reporting`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
