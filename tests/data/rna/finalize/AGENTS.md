# AGENTS.md — `MetaInformAnt/tests/data/rna/finalize`

Verified against disk 2026-08-30 (doc-realization fleet pass). RNA finalize-stage fixture data (Apis_mellifera tables) for RNA retrieval/finalize tests.

## Layout

- `PAI.md`
- Subdirectories: `Apis_mellifera`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/data/rna/finalize/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Fixture only; do not treat as analysis output.
