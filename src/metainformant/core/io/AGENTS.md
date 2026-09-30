# AGENTS.md — `MetaInformAnt/src/metainformant/core/io`

Source module under the METAINFORMANT bioinformatics toolkit (`core` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `atomic.py` — Atomic file operations for METAINFORMANT."""
- `cache.py` — JSON-based caching with TTL support for METAINFORMANT.
- `checksums.py` — File checksum and integrity verification for METAINFORMANT."""
- `disk.py` — Disk space monitoring and file system utilities for METAINFORMANT.
- `download.py` — Robust, modular download utilities with progress + heartbeat.
- `download_manager.py` — Download Manager module.
- `download_robust.py` — Robust file download utilities for MetaInformAnt.
- `errors.py` — Core I/O exceptions."""
- `io.py` — (no docstring)
- `paths.py` — (no docstring)
- `sra_environment.py` — Internal campaign-local environments for SRA command-line tools.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-core-io` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
