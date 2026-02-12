# Protein Structure

3D protein structure analysis: PDB/mmCIF file I/O, RMSD calculation via Kabsch algorithm, contact maps, secondary structure prediction, AlphaFold integration, and RCSB PDB access.

## Contents

| File | Purpose |
|------|---------|
| `io.py` | PDB/mmCIF file parsing, writing, chain extraction, format conversion |
| `pdb.py` | RCSB PDB download, atom parsing, backbone extraction, contact finding |
| `general.py` | Kabsch RMSD, radius of gyration, center of mass, inertia tensor, SASA |
| `analysis.py` | Contact maps, domain identification, surface area, ligand binding sites |
| `contacts.py` | Residue contacts, hydrogen bonds, salt bridges, disulfide bonds |
| `secondary.py` | Secondary structure prediction (PSIPRED, JPred, simple), DSSP parsing |
| `alphafold.py` | AlphaFold DB URL building, model download, confidence parsing |

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
| `predict_secondary_structure()` | Predict H/E/C assignments for a sequence |
| `build_alphafold_url()` | Construct AlphaFold DB download URL |
| `fetch_alphafold_model()` | Download AlphaFold predicted structure |
| `parse_alphafold_confidence()` | Extract per-residue pLDDT scores |

## Usage

```python
from metainformant.protein.structure.general import compute_rmsd_kabsch
from metainformant.protein.structure.alphafold import fetch_alphafold_model

rmsd = compute_rmsd_kabsch(ref_coords, mobile_coords)
path = fetch_alphafold_model("P12345", out_dir=Path("structures/"))
```
