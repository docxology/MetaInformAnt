# Personal AI Infrastructure (PAI) - docs/spatial

## Context & Intent

- **Path**: `docs/spatial` (relative to repo root)
- **Purpose**: User documentation for spatial transcriptomics: coordinate I/O,
  neighborhood analysis, autocorrelation, niche analysis, communication,
  deconvolution, integration, visualization.
- **Domain**: `docs` for `src/metainformant/spatial/` (plain arrays and
  dataclasses, no AnnData requirement; see `ARCHITECTURE.md`).

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/` (doc-tree entry point: `docs/index.md`)
- **Children**: none. Flat directory of 20 markdown files: topic guides
  (index, GETTING_STARTED, CONFIGURATION, ARCHITECTURE, CAPABILITIES,
  PERFORMANCE, TROUBLESHOOTING), method pages (autocorrelation, clustering,
  communication, deconvolution, io, neighborhood, niche, visualization,
  INTEGRATION), and AGENTS/PAI/README/SPEC — all `.md`.

## Maintenance Notes

Quoted from the sibling `AGENTS.md`:

- "Document current interfaces under `src/metainformant/spatial/`.";
  "Keep examples runnable or clearly marked as pseudocode."; examples use
  `metainformant.core.io` for data I/O; no unexercised platform claims; output under `output/`.

- `ARCHITECTURE.md` invariants: the sub-packages hardly import each other (three links:
  `analysis.spatially_variable_genes` uses `analysis.autocorrelation`,
  `integration.scrna_mapping.map_scrna_to_spatial` imports `analysis.deconvolution`
  inside the function, and `communication` and `deconvolution` share the private
  `spatial/_gene_names.py`); no analysis module imports `visualization`; loaders return
  `SpatialDataset`, `XeniumDataset` or `MERFISHDataset`, each declaring a
  `coordinate_order`; results are returned as dataclasses or dictionaries, not written
  onto a shared object.

## AI Workflows

- **Onboard**: read `index.md` -> `GETTING_STARTED.md` -> `ARCHITECTURE.md`.
- **Document a new algorithm**: update the matching method page and
  `CAPABILITIES.md`; claim platform support only with a source/test path.
  Keep imports canonical (`metainformant.spatial`) and mark pseudocode.
- **Validation**: repo docs checks + link validation; no dir-local test.
