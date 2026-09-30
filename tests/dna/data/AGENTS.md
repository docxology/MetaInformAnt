# AGENTS.md — `MetaInformAnt/tests/dna/data`

Verified against disk 2026-08-30 (doc-realization fleet pass). Data subfolder for DNA tests.

## Layout

- Subdirectories: `dna`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/dna/data/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
