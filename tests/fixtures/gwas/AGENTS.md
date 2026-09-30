# AGENTS.md — `MetaInformAnt/tests/fixtures/gwas`

Verified against disk 2026-08-30 (doc-realization fleet pass). GWAS test-fixture generator (synthetic genotype/phenotype data for the GWAS test suite).

## Layout

- `__init__.py`
- `generate_test_data.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/fixtures/gwas/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- `generate_test_data.py` is executed (not imported) to create fixture files; keep outputs deterministic.
