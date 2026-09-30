# AGENTS.md — `MetaInformAnt/tests/cstmm_fixture`

Verified against disk 2026-08-30 (doc-realization fleet pass). Shared test fixture tree for the CSTMM cross-species expression-merge pipeline (`cstmm`/`csca` steps of the hymenoptera_amalgkit workflow).

## Layout

- Subdirectories: `csca`, `cstmm`, `curate`, `merge`, `metadata`, `orthogroups`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/cstmm_fixture/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Never edit fixture contents to make a failing test pass; regenerate via the fixture generator scripts instead (unverified — locate in `scripts/`).
