# AGENTS.md — `MetaInformAnt/tests/data/gwas`

Verified against disk 2026-08-30 (doc-realization fleet pass). GWAS fixture data directory for the GWAS test suite (contents on disk as of this pass: see sibling fixtures in `tests/fixtures/gwas`).

## Layout


## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/data/gwas/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
