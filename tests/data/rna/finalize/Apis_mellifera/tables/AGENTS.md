# AGENTS.md — `MetaInformAnt/tests/data/rna/finalize/Apis_mellifera/tables`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `tables` domain of METAINFORMANT. Tests import from `src/metainformant/tables` and follow the repo's real-implementation policy.

## Layout

- `Apis_mellifera.metadata.tsv`
- `Apis_mellifera.uncorrected.tc.tsv`
- `PAI.md`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/data/rna/finalize/Apis_mellifera/tables/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
