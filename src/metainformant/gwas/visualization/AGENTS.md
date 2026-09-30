# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/visualization`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `_general_impl.py` — GWAS visualization utilities.
- `config.py` — GWAS visualization configuration and theming.
- `eqtl_visualization.py` — eQTL visualization module for expression-variant analysis.
- `general.py` — GWAS visualization utilities.
- `strain_plots.py` — Strain-aware GWAS visualization functions.
- `utils.py` — Shared utility functions for GWAS visualization.

Subpackages: `genomic`, `interactive`, `population`, `statistical`.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-visualization` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
