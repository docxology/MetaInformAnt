# AGENTS.md — `MetaInformAnt/src/metainformant/dna/population`

Source module under the METAINFORMANT bioinformatics toolkit (`dna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `analysis.py` — Advanced population genetics analysis utilities.
- `core.py` — Population genetics analysis utilities.
- `visualization.py` — Population genetics visualization utilities.
- `visualization_core.py` — Population genetics visualization utilities - core plots.
- `visualization_stats.py` — Population genetics visualization utilities - statistical comparison plots.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-dna-population` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
