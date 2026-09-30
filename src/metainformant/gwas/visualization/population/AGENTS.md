# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/visualization/population`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `geography.py` — Geographic visualization for GWAS sample and allele frequency data.
- `population.py` — Population structure visualization for GWAS.
- `population_admixture.py` — Admixture and kinship population visualization for GWAS.
- `population_pca.py` — PCA-related population structure visualization for GWAS.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-visualization-population` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
