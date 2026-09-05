# PAI - config/networks

## Context & Intent

`config/networks/` configures biological network analysis pipelines for graph construction and community detection. Per `AGENTS.md`: "Biological network analysis pipeline configuration for graph construction and community detection." It holds one file, `networks_template.yaml`, the template covering input data, analysis toggles, network construction, and outputs, consumed by `src/metainformant/networks/`.

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer; per root `SPEC.md` YAML values can be overridden via environment variables with domain prefixes).
- Sibling template dirs: `config/singlecell/`, `config/multiomics/`, `config/life_events/` — same template conventions.
- Downstream consumers: `src/metainformant/networks/` (link in `README.md`); outputs to `output/networks/`.

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the networks namespace documented by the module."

## AI Workflows

- **New config**: start from `networks_template.yaml`. Sections: `input` (tsv/csv interactions, `has_weights`, `directed`), `analysis` (metrics/communities/centrality true by default; `pathways`/`regulatory` false — pathways require a pathway database), `network` (`min_weight: 0.0`, `normalize_weights`), `output` (json/csv, `include_visualizations: true`), `performance`.
- **Verification**: load the YAML with `metainformant.core.utils.config` helpers under `uv run python`; targeted networks tests only, never repo-wide (shared-instruction USB constraint).
- **Scope discipline**: YAML here only; graph algorithm code in `src/metainformant/networks/` is off-limits for config tasks.
