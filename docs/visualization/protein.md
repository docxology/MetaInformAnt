# Protein Visualization

This document provides comprehensive documentation for protein structure and sequence visualization capabilities in METAINFORMANT.

## Overview

Protein visualization includes specialized plots for protein sequences, structures, domains, and functional analysis. These tools integrate with protein analysis workflows to create publication-quality figures.

All plotting functions live in `metainformant.protein.visualization.general`; the examples import that module as `prot_viz`. They take in-memory data (strings, NumPy arrays, lists of dicts) rather than file paths, and return the `matplotlib.axes.Axes` (the interactive viewer returns a Plotly figure). Each plotting function accepts `ax=` to draw on existing axes, `output_path=` to save the figure, and `figsize=` to size a figure it creates (the interactive viewer takes only `output_path`). Inputs are type-checked: a wrong type raises `metainformant.core.utils.errors.ValidationError` and an invalid value raises `ValueError`.

## Module Functions

### Sequence Visualization

#### Sequence Logo Plots

`plot_sequence_logo(sequences, *, ax=None, output_path=None, figsize=(10, 4), **kwargs)` draws one stacked bar per alignment position, with one segment per amino acid colored by residue (see the legend), rather than stacked letters. `sequences` is a list or tuple of aligned strings that must all have the same length.

```python
from metainformant.protein.visualization import general as prot_viz

# Create a sequence logo from equal-length aligned sequences
sequences = [
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
    "MVLTPADKSNVKAAWGKVGAHSGEYGAEAL",
    "MGLSPADKTNVKAAWGKIGGHAGEYGAEAL",
]

ax = prot_viz.plot_sequence_logo(sequences, figsize=(12, 4))
```

#### Domain Architecture

`plot_domain_architecture(domains, protein_length, *, ax=None, output_path=None, figsize=(12, 3), **kwargs)` draws the protein as a line with one rectangle per domain. `protein_length` is required and must be a Python `int`. Each domain dict needs `start` and `end` with `0 <= start < end <= protein_length`; `name` (default `'Unknown'`) and `type` (default `'domain'`) are optional, rectangles are colored by `type`, and a name is only written on domains wider than 50 residues.

```python
from metainformant.protein.visualization import general as prot_viz

# Plot protein domain architecture
domains = [
    {'start': 10, 'end': 50, 'name': 'Signal', 'type': 'signal'},
    {'start': 60, 'end': 120, 'name': 'Globular', 'type': 'domain'},
    {'start': 140, 'end': 180, 'name': 'Transmembrane', 'type': 'tm'}
]

ax = prot_viz.plot_domain_architecture(domains, protein_length=200)
```

#### Secondary Structure

`plot_secondary_structure(sequence, secondary_structure, *, ax=None, output_path=None, figsize=(12, 4), **kwargs)` draws one point per residue colored by state. `secondary_structure` has one character per residue: `H` (helix), `E` (sheet), `C` (coil) or `-` (unknown). The two strings must have the same length, and any other character is not drawn.

```python
from metainformant.protein.visualization import general as prot_viz

# Visualize secondary structure along a sequence (one state per residue)
sequence =            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFD"
secondary_structure = "CCHHHHHHHHHHHHCCEEEEECCCHHHHHHHHHHHHHHHHHHHHCCCC"

ax = prot_viz.plot_secondary_structure(sequence, secondary_structure)
```

#### Helical Wheel

`plot_helical_wheel(sequence, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` projects a segment (typically 18-25 residues) onto a circle at 100 degrees per residue, as seen down an alpha helix, and colors each residue by Kyte-Doolittle hydrophobicity.

```python
from metainformant.protein.visualization import general as prot_viz

ax = prot_viz.plot_helical_wheel("LKKLLKLLKKLLKLLKKL")
```

#### Alignment Heatmap

`plot_msa_heatmap(aligned_sequences, labels=None, *, ax=None, output_path=None, figsize=(14, 6), **kwargs)` colors each alignment cell by residue class: hydrophobic (`AILMFWVP`), charged (`DEKRH`), polar (`STNQYC`) or special (`G`). Gaps (`-`) are white and any other letter is drawn as polar. `aligned_sequences` is a list of equal-length strings and `labels` an optional list of row names.

```python
from metainformant.protein.visualization import general as prot_viz

aligned = [
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
    "MVLTPADKSNVKAAWGKVGAHSGEYGAEAL",
    "MG-SPADKTNVKAAWGKIGGHAGEYGAEAL",
]
ax = prot_viz.plot_msa_heatmap(aligned, labels=["human", "horse", "fish"])
```

### Structure Visualization

#### Contact Maps

`plot_contact_map(contacts, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` takes a square 2-D NumPy matrix (residues x residues) and draws it with the `Blues` colormap; a 0/1 contact matrix works as well as probabilities.

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

# Plot a protein contact map
contacts = np.random.default_rng(0).random((100, 100))  # 100x100 contact matrix
contacts = (contacts + contacts.T) / 2
ax = prot_viz.plot_contact_map(contacts, figsize=(8, 8))
```

#### Ramachandran Plots

`plot_ramachandran_plot(phi_angles, psi_angles, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` takes two NumPy arrays of equal length (dihedral angles in degrees), shades their density, and marks the approximate alpha-helix and beta-sheet regions.

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

# Create a Ramachandran plot
rng = np.random.default_rng(0)
phi_angles = rng.normal(-60, 30, 1000)  # Phi dihedral angles
psi_angles = rng.normal(-45, 30, 1000)  # Psi dihedral angles

ax = prot_viz.plot_ramachandran_plot(phi_angles, psi_angles)
```

### Functional Analysis

#### Alignment Quality

`plot_alignment_quality(alignment_scores, positions=None, *, ax=None, output_path=None, figsize=(10, 4), **kwargs)` draws a filled line of per-position scores from a NumPy array; `positions` (a NumPy array) defaults to the indices.

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

# Plot sequence alignment quality scores
quality_scores = np.random.default_rng(0).random(150)  # Quality scores along alignment
ax = prot_viz.plot_alignment_quality(quality_scores)
```

#### Conservation Scores

`plot_conservation_scores(conservation_scores, positions=None, *, ax=None, output_path=None, figsize=(10, 4), **kwargs)` has the same inputs as `plot_alignment_quality` and draws a filled red line.

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

# Visualize sequence conservation
conservation = np.random.default_rng(0).random(150)  # Conservation scores
positions = np.arange(150)

ax = prot_viz.plot_conservation_scores(conservation, positions)
```

#### Protein Properties

`plot_protein_properties(sequence, *, properties=None, ax=None, output_path=None, figsize=(12, 6), **kwargs)` plots per-residue scales along a sequence. Two properties are available: `'hydrophobicity'` (Kyte-Doolittle) and `'charge'` (+1 for K, R and H, -1 for D and E, 0 otherwise). With `properties=None` both are drawn, and any other name is ignored.

```python
from metainformant.protein.visualization import general as prot_viz

# Plot physicochemical properties
sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"

ax = prot_viz.plot_protein_properties(sequence, properties=['hydrophobicity', 'charge'])
```

#### Property Distributions

`plot_property_distribution(properties, *, ax=None, output_path=None, figsize=(10, 6), kind='box', **kwargs)` compares several properties across proteins. `properties` maps a property name to a list of values, and `kind` is `'box'` or `'violin'` (the violin option needs seaborn and otherwise falls back to a box plot).

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

rng = np.random.default_rng(0)
properties = {
    'hydrophobicity': rng.normal(-0.4, 0.5, 20).tolist(),
    'net charge': rng.normal(0, 3, 20).tolist(),
}
ax = prot_viz.plot_property_distribution(properties, kind='violin')
```

### Quality Control

#### PDB Structure Quality

`plot_pdb_structure_quality(b_factors, positions=None, *, ax=None, output_path=None, figsize=(10, 4), **kwargs)` draws a filled line of B-factors from a NumPy array; `positions` (a NumPy array) defaults to the indices.

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

# Plot PDB B-factors along the residues
b_factors = np.random.default_rng(0).normal(20, 10, 200)  # B-factors for residues
ax = prot_viz.plot_pdb_structure_quality(b_factors)
```

### Advanced Features

#### Interactive Structure Viewer

`create_interactive_structure_viewer(pdb_data, *, output_path=None, **kwargs)` builds a Plotly 3-D scatter of atoms and returns the Plotly figure. It requires `plotly` and raises `ImportError` without it. `pdb_data['atoms']` is a non-empty list of dicts with `x`, `y` and `z`; an optional `b_factor` colors the markers and an optional `residue_id` appears in the hover text. Extra keyword arguments go to `fig.update_layout`.

```python
from metainformant.protein.visualization import general as prot_viz

# Create an interactive 3D protein structure (requires plotly)
pdb_data = {
    'atoms': [
        {'x': 0.0, 'y': 0.0, 'z': 0.0, 'residue_id': 1, 'b_factor': 15.0},
        {'x': 1.5, 'y': 0.5, 'z': 0.2, 'residue_id': 2, 'b_factor': 18.0},
        {'x': 3.0, 'y': 1.0, 'z': 0.1, 'residue_id': 3, 'b_factor': 22.0},
    ]
}

fig = prot_viz.create_interactive_structure_viewer(pdb_data)
```

#### Structure Superposition

`plot_structure_superposition(structures, labels=None, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` draws several coordinate arrays in one 3-D scatter, one color per structure. `structures` must be a `list` of `(n_atoms, 3)` NumPy arrays, and an `ax` you pass in must be 3-D axes. The coordinates are plotted as given (as small fixed-size points) and are not aligned; superimpose them first, for example with `metainformant.protein.structure.general.align_structures_kabsch`, which returns `(aligned_coords, rotation_matrix, rmsd)`.

```python
import numpy as np
from metainformant.protein.visualization import general as prot_viz

# Visualize multiple protein structures
rng = np.random.default_rng(0)
structures = [
    rng.random((100, 3)),  # Structure 1 coordinates
    rng.random((100, 3)),  # Structure 2 coordinates
]

ax = prot_viz.plot_structure_superposition(structures, labels=['WT', 'Mutant'])
```

## Integration with Protein Module

### With Structure Analysis
```python
import numpy as np
from metainformant.protein.structure.analysis import calculate_contact_map
from metainformant.protein.structure.io import read_pdb_ca_coordinates
from metainformant.protein.visualization import general as prot_viz

# Contact map from C-alpha coordinates (residues closer than 8 angstroms are in contact)
ca_coords = np.array(read_pdb_ca_coordinates("protein.pdb"))
contacts = calculate_contact_map(ca_coords, threshold=8.0)

# Visualize results
ax = prot_viz.plot_contact_map(contacts)
```

To view a parsed structure interactively, map the atom keys that `parse_pdb_file` returns (`temp_factor`, `res_seq`) onto the keys the viewer reads:

```python
from metainformant.protein.structure.io import parse_pdb_file
from metainformant.protein.visualization import general as prot_viz

structure = parse_pdb_file("protein.pdb")
pdb_data = {
    "atoms": [
        {"x": a["x"], "y": a["y"], "z": a["z"], "b_factor": a["temp_factor"], "residue_id": a["res_seq"]}
        for a in structure["atoms"]
    ]
}
fig = prot_viz.create_interactive_structure_viewer(pdb_data)  # requires plotly
```

### With Sequence Analysis
```python
from metainformant.protein.sequence import sequences
from metainformant.protein.visualization import general as prot_viz

# Analyze protein sequence
seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
window_scores = sequences.hydropathy_score(seq, window_size=19)  # sliding-window hydropathy
print(f"GRAVY: {sequences.gravy(seq):.2f}, windows: {len(window_scores)}")

# Visualize the per-residue property scales
ax = prot_viz.plot_protein_properties(seq, properties=['hydrophobicity'])
```

### With Functional Annotation
```python
from metainformant.protein.database import interpro
from metainformant.protein.visualization import general as prot_viz

# Get domain annotations over the network (an empty list is returned if the request fails)
domains = interpro.fetch_interpro_domains("P68871")  # Hemoglobin subunit beta (147 residues)

# Visualize domain architecture
domain_list = [{'start': d['start'], 'end': d['end'],
               'name': d['name'], 'type': d['type']} for d in domains]
ax = prot_viz.plot_domain_architecture(domain_list, protein_length=147)
```

## Output Options

`output_path` saves the current pyplot figure (it calls `plt.savefig`) at 300 dpi and creates missing parent directories. That is the figure just drawn unless you created other figures in between; when you draw on your own `ax=`, save with `fig.savefig(...)` instead.

```python
from metainformant.protein.visualization import general as prot_viz

sequences = ["MVLSPADKTNVKAAWGKVGAHAGEYGAEAL", "MVLTPADKSNVKAAWGKVGAHSGEYGAEAL"]

# Save to file
ax = prot_viz.plot_sequence_logo(sequences, output_path="sequence_logo.png")
```

The interactive viewer writes an HTML file instead (it requires plotly, and the file suffix is forced to `.html`):

```python
from metainformant.protein.visualization import general as prot_viz

pdb_data = {'atoms': [{'x': 0.0, 'y': 0.0, 'z': 0.0}, {'x': 1.5, 'y': 0.5, 'z': 0.2}]}
fig = prot_viz.create_interactive_structure_viewer(pdb_data, output_path="structure.html")
```

## Performance Considerations

- **Large Structures**: For proteins >1000 residues, consider subsampling for contact maps
- **Multiple Sequences**: Sequence logos work best with 10-100 aligned sequences
- **Interactive Plots**: Require Plotly installation (the `visualization` extra) for web-based viewers
- **Memory Usage**: 3D structure visualization can be memory-intensive

## Data Sources

The plotting functions do not read files; load data with the protein module first.

- **PDB Files**: `metainformant.protein.structure.io.parse_pdb_file` returns a dict with `atoms`, `coordinates` (all atoms) and `residues`; `read_pdb_ca_coordinates` returns the C-alpha coordinates
- **FASTA Files**: `metainformant.protein.sequence.sequences.read_fasta` returns `{id: sequence}`; sequences for logos and alignment heatmaps must already be aligned to equal length
- **Domain Annotations**: `metainformant.protein.database.interpro.fetch_interpro_domains` returns dicts with `start`, `end`, `name` and `type` (plus `interpro_id`, `score` and `database`)
- **Secondary Structure**: one character per residue, `H`, `E`, `C` or `-`

## Dependencies

- **Required**: matplotlib, numpy
- **Optional**: seaborn (density shading in Ramachandran plots and the violin option of `plot_property_distribution`; plain scatter and box plots are used without it), plotly (interactive viewer)

## Examples

### Complete Protein Analysis Workflow
```python
import matplotlib.pyplot as plt
import numpy as np
from metainformant.protein.sequence import sequences
from metainformant.protein.structure.analysis import calculate_contact_map
from metainformant.protein.structure.io import read_pdb_ca_coordinates
from metainformant.protein.visualization import general as prot_viz

# Load protein data
sequence = sequences.read_fasta("protein.fasta")['protein1']
ca_coords = np.array(read_pdb_ca_coordinates("protein.pdb"))
contacts = calculate_contact_map(ca_coords, threshold=8.0)

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Sequence logo (needs several aligned sequences, so this panel is left empty)
axes[0, 0].axis("off")

# Domain architecture
domains = [{'start': 10, 'end': 50, 'name': 'Domain1', 'type': 'domain'}]
prot_viz.plot_domain_architecture(domains, len(sequence), ax=axes[0, 1])

# Contact map
prot_viz.plot_contact_map(contacts, ax=axes[1, 0])

# Properties
prot_viz.plot_protein_properties(sequence, ax=axes[1, 1])

plt.tight_layout()
plt.savefig("protein_analysis.png", dpi=300, bbox_inches='tight')
```

## Troubleshooting

### Common Issues

1. **Empty Plots**: Check input data dimensions and formats
2. **Memory Errors**: Reduce resolution or subsample large datasets
3. **Import Errors**: Install plotly (`uv sync --extra visualization`) for the interactive viewer
4. **Color Issues**: Call `matplotlib.pyplot.rcdefaults()` to reset matplotlib's default styling

### Data Format Requirements

- **Sequences**: List of strings, all same length for logos and alignment heatmaps
- **Coordinates**: Nx3 numpy arrays for 3D structures
- **Domains**: List of dicts with `start` and `end` keys (and optionally `name` and `type`)
- **Angles**: Numpy arrays in degrees for Ramachandran plots

## Related Documentation

- **[Protein Analysis](../protein/index.md)**: Core protein analysis functions for sequences and structures
- **[Contact Analysis](../protein/contacts.md)**: Residue-residue contact calculations
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive protein visualization capabilities integrated with METAINFORMANT's protein analysis workflows.
