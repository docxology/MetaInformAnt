# AGENTS.md — `MetaInformAnt/tests/core/`

tests for `metainformant.core` utilities (io, paths, config, validation, execution).
Files (verified 2026-08-29): __init__.py, test_core_atomic.py, test_core_cache.py, test_core_checksums.py, test_core_compatibility.py, test_core_comprehensive.py, test_core_config.py, test_core_db.py … (33 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/core` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.