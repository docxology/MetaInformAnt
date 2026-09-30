# AGENTS.md — `MetaInformAnt/src/metainformant/core/utils`

Source module under the METAINFORMANT bioinformatics toolkit (`core` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `config.py` — (no docstring)
- `errors.py` — Error handling and resilience utilities for METAINFORMANT.
- `hash.py` — (no docstring)
- `logging.py` — Logging utilities for METAINFORMANT.
- `optional_deps.py` — Centralized optional dependency warning management.
- `progress.py` — Progress tracking utilities for METAINFORMANT.
- `symbols.py` — Symbol indexing and cross-referencing utilities for METAINFORMANT.
- `text.py` — (no docstring)
- `timing.py` — Timing and performance utilities for METAINFORMANT."""
- `watchdog.py` — Process Watchdog Utility.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-core-utils` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
