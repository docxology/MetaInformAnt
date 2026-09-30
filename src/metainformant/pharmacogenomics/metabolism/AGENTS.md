# AGENTS.md — `MetaInformAnt/src/metainformant/pharmacogenomics/metabolism`

Source module under the METAINFORMANT bioinformatics toolkit (`pharmacogenomics` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `metabolizer_status.py` — Metabolizer status prediction from pharmacogenomic genotype data.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-pharmacogenomics-metabolism` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
