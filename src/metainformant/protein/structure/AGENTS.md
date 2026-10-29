# AGENTS.md — `src/metainformant/protein/structure`

Sub-package `structure` of the `protein` domain module
(verified 2026-08-29: 10 files, 0 subdirs).
Part of the METAINFORMANT package; keep public API in sync with docs and
README. `__pycache__/` is generated.
Repo-wide policy: see the repository-root `AGENTS.md`.

## Rules

- AlphaFold figures are per residue, never per atom: a residue's pLDDT is the mean B-factor of its ATOM records, and residues are identified by chain, number and insertion code.
- `get_alphafold_metadata()` reads local files only. Do not add online lookups to it or invent values for facts the files do not hold; an absent file means `None`.
- A function that cannot do what its name says raises `NotImplementedError` (`find_alphafold_models_by_sequence`, `search_alphafold_by_keyword`, `get_alphafold_coverage`, the `helix_turn_helix` motif) instead of returning an empty or fixed result.
- `analyze_structural_motifs()` and `identify_ligand_binding_sites()` use ATOM records only; HETATM water, ions and ligands are not part of the C-alpha trace or the pocket grid (zinc ions are read as ligands of `zinc_site` motifs).
- `identify_domains()` (spectral bisection of the contact graph), the SASA functions (Shrake-Rupley), `identify_secondary_structure_elements()` (P-SEA) and the similarity functions (Kabsch RMSD, TM-score) are the standard algorithms; keep their parameters (`structure/README.md`) in step with the code.
- `predict_secondary_structure()` and `identify_hydrogen_bonds()` stay simple approximations (propensities; distance only); `psipred` and `jpred` raise `NotImplementedError` and no function of this package may substitute another method for a requested one or call an invented service URL.
- A file reader or writer never fills in a field the data lack (chain, residue number, coordinates, element, occupancy): it skips or rejects the record, or marks the value missing (NaN, `None`).
