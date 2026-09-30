# AGENTS.md — `MetaInformAnt/src/metainformant/core/execution`

Source module under the METAINFORMANT bioinformatics toolkit (`core` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `discovery.py` — Symbolic mapping and context discovery utilities for METAINFORMANT.
- `parallel.py` — (no docstring)
- `workflow.py` — Config-driven workflow execution and orchestration for METAINFORMANT.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-core-execution` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
