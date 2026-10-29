# Protein Structure

3D protein structure analysis: PDB/mmCIF file I/O, Kabsch RMSD and TM-score, Shrake-Rupley solvent accessible surface area, P-SEA secondary structure assignment, contact-graph domain identification, contact maps, secondary structure prediction, AlphaFold integration, and RCSB PDB access.

## Contents

| File | Purpose |
|------|---------|
| `io.py` | PDB/mmCIF file parsing (the mmCIF `_atom_site` loop is read by column name), writing, chain extraction, format conversion |
| `pdb.py` | RCSB PDB download, atom parsing, backbone extraction, contact finding |
| `general.py` | Kabsch RMSD, TM-score, Shrake-Rupley SASA, P-SEA secondary structure elements, radius of gyration, center of mass, inertia tensor |
| `analysis.py` | Contact maps, spectral domain identification, structure SASA, RMSD/TM-score similarity, structural motifs, ligand binding pockets |
| `contacts.py` | Residue contacts, hydrogen bond candidates (distance only), salt bridges, disulfide bonds |
| `secondary.py` | Chou-Fasman propensity predictors (`window`, `simple`; `psipred` and `jpred` raise `NotImplementedError`), DSSP parsing |
| `alphafold.py` | AlphaFold DB URL building, model download, per-residue pLDDT and PAE parsing from local files |

## Key Functions

| Function | Description |
|----------|-------------|
| `parse_pdb_file()` | Parse PDB file into structured atom/residue data |
| `write_pdb_file()` | Write structure data back to PDB format |
| `fetch_pdb_structure()` | Download PDB/CIF from RCSB by ID |
| `compute_rmsd_kabsch()` | Optimal RMSD between two coordinate sets |
| `calculate_tm_score()` | TM-score (Zhang & Skolnick 2004) of residue-paired coordinates, with its superposition |
| `calculate_structural_similarity()` | Kabsch RMSD and TM-score of the C-alpha atoms of two structures |
| `calculate_atomic_sasa()` / `calculate_structure_sasa()` | Shrake-Rupley solvent accessible surface area per atom, and for a structure per atom and residue |
| `assign_secondary_structure_psea()` / `identify_secondary_structure_elements()` | P-SEA helices and strands from a C-alpha trace |
| `identify_domains()` | Domains by spectral bisection of the residue contact graph |
| `calculate_contact_map()` | Binary residue contact map from coordinates |
| `calculate_residue_contacts()` | Residue-level contact map with atom ranges |
| `identify_hydrogen_bonds()` | Candidate N-H...O hydrogen bonds from the donor-acceptor distance |
| `predict_secondary_structure()` | Predict H/E/C assignments for a sequence from Chou-Fasman propensities (see below) |
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

## Standard Algorithms

These replace earlier approximations (a convex-hull "surface area", numbering gaps as "domains", a single distance window as "secondary structure", and an arbitrary rescaling of an unsuperposed RMSD as "similarity").

| Function | Algorithm and parameters |
|----------|--------------------------|
| `identify_secondary_structure_elements()`, `assign_secondary_structure_psea()` (`general.py`) | P-SEA (Labesse et al. 1997) on the C-alpha trace. The distances CA(i)-CA(i+2), CA(i)-CA(i+3), CA(i)-CA(i+4) and the virtual bond angle and torsion are compared with the windows of an alpha-helix (5.5 +/- 0.5, 5.3 +/- 0.5, 6.4 +/- 0.6 A; 89 +/- 12 and 50 +/- 20 degrees) and of a beta-strand (6.7 +/- 0.6, 9.9 +/- 0.9, 12.4 +/- 1.1 A; 124 +/- 14 and -170 +/- 45 degrees) (`PSEA_HELIX_WINDOWS`, `PSEA_STRAND_WINDOWS`). A window of five residues fits a conformation when its three distances do, or when the angles and torsions at its two inner residues do; the residues of fitting windows form the elements (helix before strand, at least five residues), assigned per chain segment (chain breaks at C-alpha distances outside 2.9-4.6 A). The five-residue window rule is this implementation's reading of the method; the windows are the P-SEA parameter set (see the note at `PSEA_HELIX_WINDOWS` in `general.py` on how they were transcribed). Handedness is not tested by the distance criterion, so a mirror-image helix is a helix. The second argument is now `ca_indices` (the C-alpha atoms in `coords`), no `confidence` is reported |
| `calculate_solvent_accessible_surface_area()`, `calculate_surface_area()`, `calculate_atomic_sasa()`, `calculate_structure_sasa()` | Shrake-Rupley (1973): every atom is a sphere of its Bondi van der Waals radius (`VAN_DER_WAALS_RADII`) plus the probe (1.4 A); `n_points` (default 960) golden-spiral test points per sphere are checked against the other spheres, and the area is the accessible fraction of the sphere areas. Bare coordinates without `radii` use the carbon radius; `calculate_structure_sasa()` takes the radius from each atom's element, uses ATOM records only unless `include_hetero=True`, and refuses an element without a radius. An isolated atom is exact, two overlapping spheres match the spherical-cap formula to about 1/`n_points` |
| `identify_domains()` | Per chain, the residues are the nodes of a graph whose edges join C-alpha atoms within `contact_distance` (8 A). A connected set is ordered by the Fiedler vector of its normalised Laplacian and cut where the normalised cut `cut/vol(A) + cut/vol(B)` is smallest among cuts that leave both parts at least `min_domain_size` (40) residues; the cut is accepted when it costs at most `max_normalized_cut` (0.10), and both parts are treated the same way. Unconnected pieces are separate; a piece smaller than `min_domain_size` joins the domain nearest to it. Domains may be discontinuous (`segments`). The parameters are not fitted to a domain database: planted compact globules of 90-250 residues cost 0.18-0.32 and stay whole, two compact groups joined by a linker cost about 0.005 and split; elongated single-domain structures (long helices, coiled coils, solenoids) are split too |
| `calculate_structural_similarity()`, `calculate_structural_alignment_quality()`, `calculate_tm_score()`, `compare_structures()` (workflow) | The C-alpha atoms are paired by chain, residue number and insertion code (`pair_by="residue_id"`, or by position with `"order"`); `rmsd` is the RMSD after the optimal rigid superposition (SVD Kabsch with a reflection check); `tm_score` is the maximum over superpositions of `sum 1/(1+(d_i/d0)^2) / L` with `d0 = 1.24 (L-15)^(1/3) - 1.8` (0.5 A for `L <= 21`) and `L` the C-alpha atoms of the first structure, searched as the TM-score program does (fragment seeds refined by iterated Kabsch superposition), so it is the exact score of a real superposition and a lower bound of the maximum. No structural alignment is searched for |

## Simplified Methods

These functions are deliberately simple approximations, not substitutes for the established tools:

| Function | What it does |
|----------|--------------|
| `predict_secondary_structure()` | `"window"` (default) assigns H, E or C from the mean Chou-Fasman helix and sheet parameters over 7 residues; `"simple"` from the largest parameter of each residue. Accuracy is that of single-sequence propensity methods, not DSSP. `"psipred"` and `"jpred"` raise `NotImplementedError`: they need an external service or program that this package does not call, and nothing is substituted for them |
| `identify_hydrogen_bonds()` | Pairs of a nitrogen donor and an oxygen acceptor within `distance_threshold` (3.5 A), without the donor-hydrogen-acceptor angle (hydrogens are rarely present). Covalently close pairs (same residue, N(i+1)-O(i) across the peptide bond) and the proline nitrogen are left out; oxygen donors and nitrogen acceptors are not considered |
| `predict_transmembrane_regions()` | A Kyte-Doolittle hydropathy plot (window 19, mean above 1.5), with overlapping windows merged into one region. It is not a topology predictor |
