# AGENTS.md — `MetaInformAnt/tests/protein/`

tests for the protein module.
Files (verified 2026-08-29): __init__.py, test_protein_alignment_algorithms.py, test_protein_alphafold_fetch.py, test_protein_cli.py, test_protein_cli_comp.py, test_protein_cli_structure.py, test_protein_comprehensive.py, test_protein_contacts.py … (18 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/protein` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.