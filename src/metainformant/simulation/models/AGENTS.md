# AGENTS.md — `MetaInformAnt/src/metainformant/simulation/models`

Source module under the METAINFORMANT bioinformatics toolkit (`simulation` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `agents.py` — Agent-based ecosystem simulation utilities.
- `popgen.py` — Population genetics simulation utilities for generating synthetic genomic data.
- `rna.py` — RNA expression simulation utilities for generating synthetic transcriptomic data.
- `sequences.py` — Sequence simulation utilities for generating synthetic biological sequences.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-simulation-models` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
