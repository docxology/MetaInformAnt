# AGENTS.md — `MetaInformAnt/src/metainformant/ecology/analysis`

Source module under the METAINFORMANT bioinformatics toolkit (`ecology` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `community.py` — Community ecology analysis and biodiversity metrics.
- `functional.py` — Functional ecology analysis: trait-based diversity and community metrics.
- `indicators.py` — Ecological indicator and multivariate community analysis.
- `macroecology.py` — Macroecological analysis methods for species abundance distributions and scaling laws.
- `ordination.py` — Ordination methods for ecological community analysis.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-ecology-analysis` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
