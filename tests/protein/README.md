# protein tests

pytest suite for the `protein` domain of METAINFORMANT. Tests import from `src/metainformant/protein` and follow the repo's real-implementation policy.

## Files

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

Run from the repo root: `pytest tests/protein/ -v`. Tests follow the real-implementation policy (no mocks).
