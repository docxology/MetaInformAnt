# AGENTS.md — `MetaInformAnt/tests/other`

Verified against disk 2026-08-30 (doc-realization fleet pass). Cross-cutting tests that don't belong to one domain: domain-module imports, examples validation, import verification, orchestrators, ortholog generation, tissue config.

## Layout

- `__init__.py`
- `test_domain_modules.py`
- `test_examples.py`
- `test_import_verification.py`
- `test_orchestrators.py`
- `test_ortholog_generation.py`
- `test_tissue_config.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/other/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
