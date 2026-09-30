# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/visualization/genomic`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `genome.py` — Genome-wide visualization functions for GWAS.
- `ld.py` — Linkage disequilibrium (LD) visualization for GWAS.
- `regional.py` — Regional visualization functions for GWAS.
- `variants.py` — Variant-level visualization functions for GWAS.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-visualization-genomic` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
