# Specification: structure

## 🎯 Scope
Protein structure analysis subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 8 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `alphafold.py`
- `analysis.py`
- `contacts.py`
- `general.py`
- `io.py`
- `pdb.py`
- `secondary.py`

## Behavioral Contracts
- `alphafold.get_alphafold_metadata(uniprot_acc, *, version=4, pdb_path=None, pae_path=None)` looks nothing up online: `url` comes from the accession and version, and `confidence_score`, `n_residues`, `plddt_bands`, `mean_pae` and `max_pae` come from the local files given and are `None` otherwise. `FileNotFoundError` for a missing file, `ValueError` for an unparseable file or a PAE matrix and model with different residue counts.
- `alphafold.parse_alphafold_pae(pae_path)` reads both AlphaFold PAE layouts (full matrix, sparse residue pairs) and returns `pae_matrix`, `n_residues`, `mean_pae`, `max_pae`, `reported_max_pae` and `format`.
- `alphafold.parse_alphafold_confidence(pdb_path)` returns one pLDDT per residue (the mean B-factor of its ATOM records, first model, HETATM ignored) and `[]` for a missing file. `get_alphafold_structure_quality()` and `validate_alphafold_structure()` work per residue; `validate_alphafold_structure()` fills `n_residues`.
- `alphafold.find_alphafold_models_by_sequence()`, `search_alphafold_by_keyword()` and `get_alphafold_coverage()` raise `NotImplementedError`.
- `analysis.analyze_structural_motifs(structure, motif_types=None)` finds `"helix"` and `"extended"` segments on the C-alpha trace and `"zinc_site"` coordination sites; `"helix_turn_helix"` raises `NotImplementedError` and an unknown type raises `ValueError`.
- `analysis.identify_ligand_binding_sites(structure, *, grid_spacing=1.0, probe_radius=1.4, min_buriedness=5, ray_length=8.0, min_volume=30.0, lining_distance=4.0)` returns buried pockets from a LIGSITE-style grid search, largest first; a grid of more than 30 million cells raises `ValueError`.
- `general.identify_secondary_structure_elements(coords, ca_indices=None)` and `assign_secondary_structure_psea(ca_coords)` apply P-SEA to the C-alpha trace: elements are `{"type": "helix" | "strand", "start_residue", "end_residue", "length"}` (0-based, inclusive), states are `"H"`, `"E"`, `"C"`.
- `general.calculate_atomic_sasa(coords, radii=None, probe_radius=1.4, n_points=960)` returns the Shrake-Rupley area of every atom; `calculate_solvent_accessible_surface_area()` and `analysis.calculate_surface_area()` return the total; `analysis.calculate_structure_sasa(structure, *, probe_radius=1.4, n_points=960, include_hetero=False)` returns `total`, `per_atom`, `atom_indices`, `per_residue`, `n_atoms`, `probe_radius` and `n_points`.
- `analysis.identify_domains(structure, *, contact_distance=8.0, min_domain_size=40, max_normalized_cut=0.10)` returns domains with `domain_id`, `chain_id`, `residues`, `residue_ids`, `start`, `end`, `size`, `segments` and `n_contacts`; atoms without coordinates raise `ValueError`.
- `general.calculate_tm_score(coords_ref, coords_mobile, n_reference=None)` returns `tm_score`, `d0`, `n_aligned`, `n_reference`, `rotation` and `translation`; `tm_score_d0(L)` is `1.24 (L-15)^(1/3) - 1.8`, 0.5 for `L <= 21`. `analysis.calculate_structural_similarity(s1, s2, *, pair_by="residue_id")` returns `rmsd`, `tm_score`, `d0`, `n_aligned`, `n_residues_1`, `n_residues_2`; `calculate_structural_alignment_quality()` adds `quality_score` (= `tm_score`). Structures that cannot be compared raise `ValueError`.
- `secondary.predict_secondary_structure(sequence, method="window")` accepts `"window"` and `"simple"`; `"psipred"` and `"jpred"` raise `NotImplementedError`. `io.parse_mmcif_file()` reads the `_atom_site` loop by column name (first model) and raises `ValueError` where it cannot; `io.write_pdb_file()` raises `ValueError` for an atom without name, residue name or number, or coordinates.
- Simple approximations remain (see README.md): sequence-propensity secondary structure, distance-only hydrogen bonds, hydropathy transmembrane regions.
