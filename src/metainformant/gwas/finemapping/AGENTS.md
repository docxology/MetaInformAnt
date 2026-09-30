# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/finemapping`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `colocalization.py` — Multi-trait colocalization analysis for GWAS fine-mapping.
- `credible_sets.py` — Statistical fine-mapping: credible sets, SuSiE, Bayes factors, and colocalization.
- `eqtl.py` — eQTL analysis methods for gene expression-variant associations.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-finemapping` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
