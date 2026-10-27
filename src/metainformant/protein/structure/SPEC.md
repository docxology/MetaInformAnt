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
- Secondary-structure, SASA, domain and similarity functions are simplified approximations (see README.md).
