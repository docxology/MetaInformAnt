# AGENTS.md — `MetaInformAnt/config/`

YAML/TOML/JSON configuration for all METAINFORMANT domain modules plus the
amalgkit campaign configs (verified 2026-08-29 entries: AGENTS.md, PAI.md, README.md, SPEC.md, amalgkit, config_base, eqtl, gwas, life_events, longread, multiomics, ncbi, networks, phenotype, quality, singlecell).
Each subfolder already carries its own `AGENTS.md`/`README.md` except
`quality/` (documented during the 2026-08-29 doc pass).

## Gotchas
- Amalgkit configs under `amalgkit/` (28 species YAMLs plus template/test/cross-species) drive the campaign; keep them
  aligned with `projects/hymenoptera_amalgkit/` evidence rules.
Repo-wide policy: see the repository-root `AGENTS.md`.