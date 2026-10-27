# AGENTS.md — `MetaInformAnt/tests/protein/`

tests for the protein module.
Files (verified 2026-10-27): __init__.py, test_protein_alignment_algorithms.py, test_protein_alphafold_fetch.py, test_protein_alphafold_local_files.py, test_protein_alphafold_per_residue.py, test_protein_cli.py, test_protein_cli_comp.py, test_protein_cli_structure.py, test_protein_comprehensive.py, test_protein_contact_types.py, test_protein_contacts.py, test_protein_domains_regression.py, test_protein_enhancements.py, test_protein_identity_alignment.py, test_protein_interpro.py, test_protein_proteomes.py, test_protein_proteomes_api.py, test_protein_sequence_logo.py, test_protein_sequences.py, test_protein_structure_example_randomness.py, test_protein_structure_io_rmsd.py, test_protein_structure_motifs_and_pockets.py, test_protein_structure_rmsd.py, test_protein_structure_secondary.py, test_protein_uniprot_pdb.py (24 test files).

## Contract tests

- `test_protein_alphafold_local_files.py`: `get_alphafold_metadata` reads pLDDT and PAE from local files only (`None` for a file not given), `parse_alphafold_pae` reads both PAE layouts, and `find_alphafold_models_by_sequence`, `search_alphafold_by_keyword` and `get_alphafold_coverage` raise `NotImplementedError`.
- `test_protein_alphafold_per_residue.py`: confidence figures are per residue, not per atom: `parse_alphafold_confidence` returns one mean pLDDT per residue (first model, HETATM ignored), `get_alphafold_structure_quality` and `validate_alphafold_structure` count residues, and `validate_alphafold_structure` fills `n_residues`.
- `test_protein_structure_motifs_and_pockets.py`: `analyze_structural_motifs` finds helices, extended segments and zinc sites from coordinates (`helix_turn_helix` raises `NotImplementedError`), and `identify_ligand_binding_sites` finds buried pockets with the documented parameters.
- `test_protein_structure_example_randomness.py`: the protein-structure modules do not draw from the global NumPy or Python random state (their docstring examples use a seeded generator).

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/protein`
Repo-wide policy: see the repository-root `AGENTS.md`.
