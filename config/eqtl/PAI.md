# PAI - config/eqtl

## Context & Intent

`config/eqtl/` configures eQTL integration pipelines that link GWAS variants with RNA-seq expression data. Per `AGENTS.md`, its role is "eQTL integration pipeline configuration linking GWAS variants with RNA-seq expression data". It currently holds one real config, `eqtl_amellifera.yaml`, an *Apis mellifera* RNA-SNP pipeline config used by `scripts/eqtl/rna_snp_pipeline.py`.

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer; environment overrides use domain prefixes per root `SPEC.md`).
- Upstream inputs referenced by `eqtl_amellifera.yaml`: `output/amalgkit/amellifera/...` (expression) and GWAS variant outputs.
- Sibling config domains: `config/gwas/`, `config/amalgkit/` (explicitly cross-linked in `README.md`).
- Downstream consumers: `scripts/eqtl/rna_snp_pipeline.py` and `src/metainformant/eqtl/`.

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the eQTL namespace documented by the module."

## AI Workflows

- **New config**: copy `eqtl_amellifera.yaml` as a base; keep sections `samples`, `alignment`, `variant_calling`, `filtering`, `output`, `eqtl_scan`; validate against the eqtl schema before committing.
- **Smoke-check a config**: run the `load_mapping_from_file` snippet in `README.md` via `uv run python` and confirm top-level keys load.
- **Sample discovery**: default `samples.mode: auto` with `max_samples: 5` auto-discovers completed Amalgkit quant runs; switch to `explicit` with SRR ids for deterministic runs.
- **Scope discipline**: never edit `src/` or `scripts/` from this directory's tasks; config data only.
