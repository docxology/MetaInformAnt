# Personal AI Infrastructure (PAI) - docs/spatial

## Context & Intent

- **Path**: `docs/spatial/` (relative to repo root)
- **Purpose**: User documentation for spatial transcriptomics: coordinate I/O,
  neighborhood analysis, autocorrelation, niche analysis, communication,
  deconvolution, integration, visualization.
- **Domain**: `docs` for `src/metainformant/spatial/` (AnnData-based;
  see `ARCHITECTURE.md`).

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/` (doc-tree entry point: `docs/index.md`)
- **Children**: none. Flat directory of 21 markdown files: topic guides
  (index, GETTING_STARTED, CONFIGURATION, ARCHITECTURE, CAPABILITIES,
  PERFORMANCE, TROUBLESHOOTING), method pages (autocorrelation, clustering,
  communication, deconvolution, integration, io, neighborhood, niche,
  visualization, INTEGRATION), and AGENTS/PAI/README/SPEC — all `.md`.

## Maintenance Notes

Quoted from the sibling `AGENTS.md`:

- "Document current interfaces under `src/metainformant/spatial/`.";
  "Keep examples runnable or clearly marked as pseudocode."; examples use
  `metainformant.core.io` for data I/O; no unexercised platform claims; output under `output/`.

- `ARCHITECTURE.md` invariants: dependencies flow upward only; no analysis
  module imports `visualization`; loaders build AnnData with
  `obsm['spatial']`/`uns['platform']`; results land in `adata.obs`/`adata.uns`.

## AI Workflows

- **Onboard**: read `index.md` -> `GETTING_STARTED.md` -> `ARCHITECTURE.md`.
- **Document a new algorithm**: update the matching method page and
  `CAPABILITIES.md`; claim platform support only with a source/test path.
  Keep imports canonical (`metainformant.spatial`) and mark pseudocode.
- **Validation**: repo docs checks + link validation; no dir-local test.
