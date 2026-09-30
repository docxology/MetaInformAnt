# AGENTS.md — `MetaInformAnt/tests/core`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `core` domain of METAINFORMANT. Tests import from `src/metainformant/core` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_core_atomic.py`
- `test_core_cache.py`
- `test_core_checksums.py`
- `test_core_compatibility.py`
- `test_core_comprehensive.py`
- `test_core_config.py`
- `test_core_db.py`
- `test_core_discovery.py`
- `test_core_discovery_cache.py`
- `test_core_disk.py`
- `test_core_download.py`
- `test_core_errors.py`
- `test_core_functionality.py`
- … (+19 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/core/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
