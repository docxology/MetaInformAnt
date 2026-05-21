# Agent Directives: metagenomics

**Context**: Microbiome and metagenomic analysis: amplicon profiling (16S/ITS), shotgun metagenomics, community diversity, functional annotation, and differential abundance testing.

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
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/metagenomics/](../../../docs/metagenomics/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Ecology module**: [../ecology/AGENTS.md](../ecology/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
