# Agent Directives: pharmacogenomics

**Context**: Clinical pharmacogenomic analysis: star allele calling, metabolizer phenotyping, CPIC guideline lookups, drug interaction prediction, and report generation.

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

## Related Documentation

- **Module guide**: [../../../docs/pharmacogenomics/](../../../docs/pharmacogenomics/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Protein module**: [../protein/AGENTS.md](../protein/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
