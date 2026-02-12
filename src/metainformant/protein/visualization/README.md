# Protein Visualization

Comprehensive plotting for protein data: sequence logos, domain architecture diagrams, contact maps, Ramachandran plots, structure quality metrics, and interactive 3D viewers.

## Contents

| File | Purpose |
|------|---------|
| `general.py` | All protein visualization functions (matplotlib, seaborn, plotly) |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_sequence_logo()` | Sequence logo from aligned protein sequences |
| `plot_domain_architecture()` | Domain layout diagram along protein length |
| `plot_secondary_structure()` | H/E/C secondary structure bar plot |
| `plot_contact_map()` | Residue-residue contact map heatmap |
| `plot_ramachandran_plot()` | Phi/psi backbone dihedral angle scatter |
| `plot_alignment_quality()` | Per-position alignment quality scores |
| `plot_conservation_scores()` | Residue conservation across an MSA |
| `plot_structure_superposition()` | Overlay two aligned structures |
| `plot_protein_properties()` | Multi-panel physicochemical property tracks |
| `plot_pdb_structure_quality()` | B-factor and resolution quality metrics |
| `create_interactive_structure_viewer()` | Plotly-based 3D structure viewer |
| `plot_helical_wheel()` | Helical wheel projection for amphipathic helices |
| `plot_msa_heatmap()` | Heatmap visualization of multiple sequence alignment |

## Usage

```python
from metainformant.protein.visualization.general import plot_contact_map, plot_ramachandran_plot

fig = plot_contact_map(contact_matrix, output_path="output/contacts.png")
fig = plot_ramachandran_plot(phi_angles, psi_angles, output_path="output/rama.png")
```
