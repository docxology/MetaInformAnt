# Protein Structure

3D protein structure analysis: PDB/mmCIF file I/O, RMSD calculation via Kabsch algorithm, contact maps, secondary structure prediction, AlphaFold integration, and RCSB PDB access.

## Contents

| File | Purpose |
|------|---------|
| `io.py` | PDB/mmCIF file parsing, writing, chain extraction, format conversion |
| `pdb.py` | RCSB PDB download, atom parsing, backbone extraction, contact finding |
| `general.py` | Kabsch RMSD, radius of gyration, center of mass, inertia tensor, simplified SASA and secondary-structure helpers |
| `analysis.py` | Contact maps, simplified domain identification and surface area, structural motifs, ligand binding pockets |
| `contacts.py` | Residue contacts, hydrogen bonds, salt bridges, disulfide bonds |
| `secondary.py` | Secondary structure assignment (`simple` propensity rules; `psipred` and `jpred` with a local fallback), DSSP parsing |
| `alphafold.py` | AlphaFold DB URL building, model download, per-residue pLDDT and PAE parsing from local files |

## Key Functions

| Function | Description |
|----------|-------------|
| `parse_pdb_file()` | Parse PDB file into structured atom/residue data |
| `write_pdb_file()` | Write structure data back to PDB format |
| `fetch_pdb_structure()` | Download PDB/CIF from RCSB by ID |
| `compute_rmsd_kabsch()` | Optimal RMSD between two coordinate sets |
| `calculate_contact_map()` | Binary residue contact map from coordinates |
| `calculate_residue_contacts()` | Residue-level contact map with atom ranges |
| `identify_hydrogen_bonds()` | Detect hydrogen bond donor-acceptor pairs |
| `predict_secondary_structure()` | Predict H/E/C assignments for a sequence (simplified, see below) |
| `analyze_structural_motifs()` | Helices and extended segments from the C-alpha trace, and zinc coordination sites |
| `identify_ligand_binding_sites()` | Buried pockets found by a LIGSITE-style grid search |
| `build_alphafold_url()` | Construct AlphaFold DB download URL |
| `fetch_alphafold_model()` | Download AlphaFold predicted structure |
| `parse_alphafold_confidence()` | One pLDDT per residue (the mean over the residue's atoms) from a model file |
| `parse_alphafold_pae()` | Predicted aligned error matrix of a PAE JSON file, with its mean and maximum |
| `get_alphafold_metadata()` | URL plus pLDDT, residue count, pLDDT bands and PAE read from local model files; no lookup online |

## Usage

```python
from metainformant.protein.structure.general import compute_rmsd_kabsch
from metainformant.protein.structure.alphafold import fetch_alphafold_model

rmsd = compute_rmsd_kabsch(ref_coords, mobile_coords)
path = fetch_alphafold_model("P12345", out_dir=Path("structures/"))
```

## AlphaFold Files

`parse_alphafold_confidence()` returns one score per residue: the mean B-factor of the residue's ATOM records of the first model (AlphaFold stores a residue's pLDDT in every atom), with HETATM records ignored, an empty list for a missing file, and residues told apart by chain, number and insertion code. `get_alphafold_metadata(uniprot_acc, version=4, pdb_path=None, pae_path=None)` reads the pLDDT (`confidence_score`, `plddt_bands`, `n_residues`) from `pdb_path` and `mean_pae` and `max_pae` from `pae_path` (parsed by `parse_alphafold_pae()`); a value is `None` when its file is not given, and facts that only the AlphaFold service holds (organism, gene, release date) are not returned. `get_alphafold_structure_quality()` fills the per-residue counts and the confidence distribution, and `validate_alphafold_structure()` fills `n_residues` (distinct chain, residue number and insertion code) and `avg_confidence` (over residues). `find_alphafold_models_by_sequence()`, `search_alphafold_by_keyword()` and `get_alphafold_coverage()` raise `NotImplementedError`: they would need a search or statistics service, and an empty result would claim that nothing exists.

## Structural Motifs and Binding Pockets

`analyze_structural_motifs(structure, motif_types=None)` searches the atoms of a structure dictionary for three motif types (`SUPPORTED_MOTIF_TYPES`):

- `"helix"`: right-handed alpha-helices on the C-alpha trace. A window of five residues counts when CA(i)-CA(i+3) and CA(i+1)-CA(i+4) are 4.5-5.7 A, CA(i)-CA(i+4) is 5.5-6.9 A and the virtual torsions of CA(i..i+3) and CA(i+1..i+4) are 20-80 degrees; overlapping windows merge and a helix has at least 6 residues. A left-handed helix is not reported.
- `"extended"`: strand-like segments of at least 4 residues, from windows of four residues with CA(i)-CA(i+2) of 6.0-7.5 A and CA(i)-CA(i+3) of 8.8-11.0 A. Pairing into sheets is not tested, so an extended loop is reported too.
- `"zinc_site"`: zinc ions with the Cys SG, His ND1/NE2 and Asp/Glu carboxylate oxygen atoms within 2.8 A, giving a `coordination_pattern` such as `"C2H2"` and `zinc_finger_like` for four Cys/His ligands.

Chains are split wherever consecutive C-alpha atoms are not 2.9-4.6 A apart, and only ATOM records supply the trace and the ligands. `"helix_turn_helix"` is not detected: asking for it raises `NotImplementedError`, and an unknown type raises `ValueError`.

`identify_ligand_binding_sites(structure, *, grid_spacing=1.0, probe_radius=1.4, min_buriedness=5, ray_length=8.0, min_volume=30.0, lining_distance=4.0)` follows the LIGSITE approach. The protein is placed on a grid; a cell is occupied when it lies within the van der Waals radius of an atom and accessible when a probe of `probe_radius` centred on it touches no atom. Each cell is scanned along the three axes and the four body diagonals, and a direction counts as buried when occupied cells are met on both sides within `ray_length`. Accessible cells buried in at least `min_buriedness` of the seven directions form pockets (26-neighbourhood), ranked by volume; a pocket smaller than `min_volume` (A^3) is dropped. Each pocket has `rank`, `center`, `volume`, `n_grid_points`, `mean_buriedness`, `lining_atoms` and `lining_residues` (atoms within `lining_distance`). Only ATOM records that are not water are used, so a pocket that holds a bound ligand is still found. Volumes are cell counts times `grid_spacing**3`, so they are estimates, and a grid of more than 30 million cells raises `ValueError` (increase `grid_spacing`).

## Simplified Methods

These functions are deliberately simple approximations, not substitutes for the established tools:

| Function | What it does |
|----------|--------------|
| `predict_secondary_structure()` | `"simple"` assigns H, E or C from the residue type alone; `"psipred"` and `"jpred"` submit the sequence to a web service and fall back to a windowed propensity predictor when the service gives no prediction. It is not DSSP |
| `identify_secondary_structure_elements()` (`general.py`) | Reports a five-atom helix element wherever the distance between backbone atoms i and i+3 lies between 5 and 7 A |
| `calculate_surface_area()` and `calculate_solvent_accessible_surface_area()` | Convex-hull area of the coordinates (`calculate_surface_area()` scales it by `1 + probe_radius / 10` and treats fewer than 4 points as a sphere), not a rolling-probe solvent-accessible surface |
| `identify_domains()` | Starts a new domain at a chain change or a jump of more than 10 in the residue number; it uses residue numbers only, not the coordinates |
| `calculate_structural_similarity()` and `calculate_structural_alignment_quality()` | Position-wise RMSD without superposition mapped to a score (`1 / (1 + rmsd / 10)` and `max(0, 1 - rmsd / 20)`) |
