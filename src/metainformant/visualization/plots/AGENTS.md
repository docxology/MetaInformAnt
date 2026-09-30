# AGENTS.md — `MetaInformAnt/src/metainformant/visualization/plots`

Source module under the METAINFORMANT bioinformatics toolkit (`visualization` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `animations.py` — Animated visualization functions for dynamic data exploration.
- `basic.py` — Basic plotting functions for fundamental chart types.
- `cross_species.py` — Cross-species visualization module.
- `general.py` — Enhanced plotting functions for biological data visualization.
- `multidim.py` — Multi-dimensional data visualization functions.
- `specialized.py` — Specialized visualization functions for advanced plotting types.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-visualization-plots` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
