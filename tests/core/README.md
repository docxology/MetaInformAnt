# core tests

pytest suite for the `core` domain of METAINFORMANT. Tests import from `src/metainformant/core` and follow the repo's real-implementation policy.

## Files

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

Run from the repo root: `pytest tests/core/ -v`. Tests follow the real-implementation policy (no mocks).
