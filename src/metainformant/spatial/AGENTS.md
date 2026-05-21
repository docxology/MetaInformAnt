# Agent Directives: spatial

**Context**: Spatial transcriptomics analysis: platform I/O (Visium, MERFISH, Xenium), spatial statistics, cell-cell communication, deconvolution, and scRNA-seq integration.

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
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/spatial/](../../../docs/spatial/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Singlecell module**: [../singlecell/AGENTS.md](../singlecell/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
