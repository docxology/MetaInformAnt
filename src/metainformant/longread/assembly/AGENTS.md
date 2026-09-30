# AGENTS.md — `MetaInformAnt/src/metainformant/longread/assembly`

Source module under the METAINFORMANT bioinformatics toolkit (`longread` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `consensus.py` — Consensus sequence generation for long-read assembly.
- `hybrid.py` — Hybrid assembly combining long and short reads.
- `overlap.py` — Overlap computation for long-read assembly using minimizer sketching.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-longread-assembly` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
