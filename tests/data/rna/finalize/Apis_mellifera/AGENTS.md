# AGENTS.md — `MetaInformAnt/tests/data/rna/finalize/Apis_mellifera`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `Apis_mellifera` domain of METAINFORMANT. Tests import from `src/metainformant/Apis_mellifera` and follow the repo's real-implementation policy.

## Layout

- `PAI.md`
- Subdirectories: `tables`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/data/rna/finalize/Apis_mellifera/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
