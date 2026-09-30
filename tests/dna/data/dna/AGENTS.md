# AGENTS.md — `MetaInformAnt/tests/dna/data/dna`

Verified against disk 2026-08-30 (doc-realization fleet pass). Small DNA test data (toy FASTA) used by DNA-domain tests.

## Layout

- `toy.fasta`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/dna/data/dna/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Keep the file minimal and deterministic.
