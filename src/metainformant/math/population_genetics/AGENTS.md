# AGENTS.md — `MetaInformAnt/src/metainformant/math/population_genetics`

Source module under the METAINFORMANT bioinformatics toolkit (`math` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `coalescent.py` — Coalescent theory mathematical models.
- `core.py` — Population genetics mathematical models and calculations.
- `demography.py` — Population demography and growth models.
- `effective_size.py` — Effective population size calculations.
- `fst.py` — F-statistics and population differentiation functions.
- `ld.py` — Linkage disequilibrium functions.
- `selection.py` — Selection theory functions.
- `statistics.py` — Population genetics statistical functions.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-math-population_genetics` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
