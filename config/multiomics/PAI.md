# PAI - config/multiomics

## Context & Intent

`config/multiomics/` configures multi-omics integration pipelines for cross-platform data harmonization. Per `AGENTS.md`: "Multi-omics integration pipeline configuration for cross-platform data harmonization." It holds one file, `multiomics_template.yaml`, the documented template for input sources, integration methods, sample mapping, and outputs, consumed by `src/metainformant/multiomics/`.

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer; per root `SPEC.md`, YAML values can be overridden via environment variables with domain prefixes).
- Sibling template dirs: `config/singlecell/`, `config/networks/`, `config/life_events/` — same conventions.
- Downstream consumers: `src/metainformant/multiomics/` (link in `README.md`); outputs to `output/multiomics/`.

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the multi-omics namespace documented by the module."

## AI Workflows

- **New config**: start from `multiomics_template.yaml`. Sections: `inputs` (genomics/transcriptomics/proteomics/metabolomics, each with `file` and `required` flag), `integration` (`joint_pca`, `joint_nmf`, `canonical_correlation`, `n_components: 10`), `samples` (`id_column: sample_id`, `require_all_omics`), `output`, `performance`.
- **Layer discipline**: omics layers default to `required: false`; set `require_all_omics: true` only when the analysis design demands complete cases.
- **Verification**: load the YAML with `metainformant.core.utils.config` helpers under `uv run python`; targeted multiomics tests only, never repo-wide (shared-instruction USB constraint).
- **Scope discipline**: YAML here only; no edits to `src/metainformant/multiomics/` from config tasks.
