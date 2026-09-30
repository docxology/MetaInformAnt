# AGENTS.md — `MetaInformAnt/tests/integration`

Verified against disk 2026-08-30 (doc-realization fleet pass). Cross-domain integration tests: eQTL integration/scripts, RNA→GWAS handoff, comprehensive integration.

## Layout

- `__init__.py`
- `test_eqtl_integration.py`
- `test_eqtl_scripts.py`
- `test_integration_comprehensive.py`
- `test_rna_gwas_handoff.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/integration/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- These exercise real handoffs between domains; run them after domain suites.
