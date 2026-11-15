# AGENTS.md — `MetaInformAnt/config/`

YAML/TOML/JSON configuration for all METAINFORMANT domain modules plus the
amalgkit campaign configs (verified 2026-08-29 entries: AGENTS.md, PAI.md, README.md, SPEC.md, amalgkit, config_base, eqtl, gwas, life_events, longread, multiomics, ncbi, networks, phenotype, quality, singlecell).
Each subfolder carries its own `AGENTS.md` and `README.md`.

## Gotchas

- Amalgkit configs under `amalgkit/` (species YAMLs plus template/test/cross-species) are the repository reference
  set; the campaign runs `projects/hymenoptera_amalgkit/config/amalgkit/` against the frozen cohort. Keep both
  aligned with `projects/hymenoptera_amalgkit/` evidence rules.

Repo-wide policy: see the repository-root `AGENTS.md`.