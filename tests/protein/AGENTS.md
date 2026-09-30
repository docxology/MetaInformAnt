# AGENTS.md — `MetaInformAnt/tests/protein`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `protein` domain of METAINFORMANT. Tests import from `src/metainformant/protein` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_protein_alignment_algorithms.py`
- `test_protein_alphafold_fetch.py`
- `test_protein_cli.py`
- `test_protein_cli_comp.py`
- `test_protein_cli_structure.py`
- `test_protein_comprehensive.py`
- `test_protein_contacts.py`
- `test_protein_enhancements.py`
- `test_protein_identity_alignment.py`
- `test_protein_interpro.py`
- `test_protein_proteomes.py`
- `test_protein_proteomes_api.py`
- `test_protein_sequences.py`
- … (+4 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/protein/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
