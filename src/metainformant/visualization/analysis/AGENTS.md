# AGENTS.md — `MetaInformAnt/src/metainformant/visualization/analysis`

Source module under the METAINFORMANT bioinformatics toolkit (`visualization` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `dimred.py` — Dimensionality reduction visualization functions.
- `information.py` — Information theory visualization functions.
- `quality.py` — Quality control data visualization functions.
- `quality_assessment.py` — Quality assessment visualization for coverage, errors, and batch effects.
- `quality_omics.py` — Multi-omics quality control visualization functions.
- `quality_sequencing.py` — Sequencing quality control visualization functions.
- `statistical.py` — Statistical plotting functions for data analysis and visualization.
- `timeseries.py` — Time series analysis visualization functions.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-visualization-analysis` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
