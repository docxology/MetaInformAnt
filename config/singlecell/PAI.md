# PAI - config/singlecell

## Context & Intent

`config/singlecell/` configures single-cell RNA-seq analysis pipelines covering preprocessing, clustering, and trajectory work. Per `AGENTS.md`: "Single-cell RNA-seq analysis pipeline configuration for preprocessing, clustering, and trajectory inference." It holds one file, `singlecell_template.yaml`, the template for input data, QC, normalization, dimensionality reduction, clustering, and outputs, consumed by `src/metainformant/singlecell/`.

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer; per root `SPEC.md` YAML values can be overridden via environment variables with domain prefixes).
- Sibling template dirs: `config/networks/`, `config/multiomics/`, `config/life_events/` — same conventions.
- Downstream consumers: `src/metainformant/singlecell/` (link in `README.md`); outputs to `output/singlecell/`.

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the single-cell namespace documented by the module."

## AI Workflows

- **New config**: start from `singlecell_template.yaml`. Sections: `input` (h5ad/csv/tsv), `qc` (`min_genes: 200`, `min_cells: 3`, `max_mito_percent: 5.0`, `max_ribo_percent: 50.0`), `normalization` (`method: log`, `target_sum: 10000`), `dimension_reduction` (50 PCA components, 15 UMAP neighbors), `clustering` (`leiden` or `louvain`, `resolution: 0.5`, `find_markers: true`), `output`, `performance`.
- **Verification**: load the YAML with `metainformant.core.utils.config` helpers under `uv run python`; targeted singlecell tests only, never repo-wide (shared-instruction USB constraint).
- **Scope discipline**: YAML here only; no edits to `src/metainformant/singlecell/` from config tasks.
