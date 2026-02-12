# Protein Visualization

This document provides comprehensive documentation for protein structure and sequence visualization capabilities in METAINFORMANT.

## Overview

Protein visualization includes specialized plots for protein sequences, structures, domains, and functional analysis. These tools integrate with protein analysis workflows to create publication-quality figures.

## Module Functions

### Sequence Visualization

#### Sequence Logo Plots
```python
from metainformant.protein import visualization as prot_viz
import numpy as np

# Create sequence logo from aligned sequences
sequences = [
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    # ... more aligned sequences
]

ax = prot_viz.plot_sequence_logo(sequences, figsize=(12, 4))
```

#### Domain Architecture
```python
# Plot protein domain architecture
domains = [
    {'start': 10, 'end': 50, 'name': 'Signal', 'type': 'signal'},
    {'start': 60, 'end': 120, 'name': 'Globular', 'type': 'domain'},
    {'start': 140, 'end': 180, 'name': 'Transmembrane', 'type': 'tm'}
]

ax = prot_viz.plot_domain_architecture(domains, protein_length=200)
```

#### Secondary Structure
```python
# Visualize secondary structure along sequence
sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
secondary_structure = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"

ax = prot_viz.plot_secondary_structure(sequence, secondary_structure)
```

### Structure Visualization

#### Contact Maps
```python
# Plot protein contact map
contacts = np.random.rand(100, 100)  # 100x100 contact matrix
ax = prot_viz.plot_contact_map(contacts, figsize=(8, 8))
```

#### Ramachandran Plots
```python
# Create Ramachandran plot
phi_angles = np.random.normal(0, 30, 1000)  # Phi dihedral angles
psi_angles = np.random.normal(0, 30, 1000)  # Psi dihedral angles

ax = prot_viz.plot_ramachandran_plot(phi_angles, psi_angles)
```

### Functional Analysis

#### Alignment Quality
```python
# Plot sequence alignment quality scores
quality_scores = np.random.rand(150)  # Quality scores along alignment
ax = prot_viz.plot_alignment_quality(quality_scores)
```

#### Conservation Scores
```python
# Visualize sequence conservation
conservation = np.random.rand(150)  # Conservation scores
positions = np.arange(150)

ax = prot_viz.plot_conservation_scores(conservation, positions)
```

#### Protein Properties
```python
# Plot physicochemical properties
sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"

ax = prot_viz.plot_protein_properties(sequence, properties=['hydrophobicity', 'charge'])
```

### Quality Control

#### PDB Structure Quality
```python
# Plot PDB B-factor distribution
b_factors = np.random.normal(20, 10, 200)  # B-factors for residues
ax = prot_viz.plot_pdb_structure_quality({'b_factors': b_factors})
```

### Advanced Features

#### Interactive Structure Viewer
```python
# Create interactive 3D protein structure
pdb_data = {
    'atoms': [
        {'x': 0.0, 'y': 0.0, 'z': 0.0, 'residue_id': 1, 'b_factor': 15.0},
        # ... more atom data
    ]
}

fig = prot_viz.create_interactive_structure_viewer(pdb_data)
```

#### Structure Superposition
```python
# Visualize multiple protein structures
structures = [
    np.random.rand(100, 3),  # Structure 1 coordinates
    np.random.rand(100, 3),  # Structure 2 coordinates
]

ax = prot_viz.plot_structure_superposition(structures, labels=['WT', 'Mutant'])
```

## Integration with Protein Module

### With Structure Analysis
```python
from metainformant.protein import structure, visualization as prot_viz

# Load and analyze protein structure
pdb_data = structure.parse_pdb_file("protein.pdb")
contacts = structure.calculate_contact_map(pdb_data['coordinates'])

# Visualize results
ax = prot_viz.plot_contact_map(contacts)
```

### With Sequence Analysis
```python
from metainformant.protein import sequences, visualization as prot_viz

# Analyze protein sequence
seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
properties = sequences.calculate_hydropathy_score(seq)

# Visualize properties
ax = prot_viz.plot_protein_properties(seq, properties=['hydrophobicity'])
```

### With Functional Annotation
```python
from metainformant.protein import interpro, visualization as prot_viz

# Get domain annotations
domains = interpro.fetch_interpro_domains("P68871")  # Hemoglobin

# Visualize domain architecture
domain_list = [{'start': d['start'], 'end': d['end'],
               'name': d['name'], 'type': d['type']} for d in domains]
ax = prot_viz.plot_domain_architecture(domain_list, protein_length=200)
```

## Output Options

All visualization functions support:
```python
# Save to file
ax = prot_viz.plot_sequence_logo(sequences, output_path="sequence_logo.png")

# Interactive web viewer
fig = prot_viz.create_interactive_structure_viewer(pdb_data, output_path="structure.html")
```

## Performance Considerations

- **Large Structures**: For proteins >1000 residues, consider subsampling for contact maps
- **Multiple Sequences**: Sequence logos work best with 10-100 aligned sequences
- **Interactive Plots**: Require Plotly installation for web-based viewers
- **Memory Usage**: 3D structure visualization can be memory-intensive

## File Format Support

- **PDB Files**: Standard Protein Data Bank format
- **FASTA Files**: Protein sequence alignments
- **Domain Annotations**: InterPro, Pfam, and custom formats
- **Secondary Structure**: DSSP and prediction formats

## Dependencies

- **Required**: matplotlib, numpy
- **Optional**: seaborn (enhanced styling), plotly (interactive plots)
- **Integration**: BioPython (sequence handling), NetworkX (structure analysis)

## Examples

### Complete Protein Analysis Workflow
```python
from metainformant.protein import sequences, structure, visualization as prot_viz
import numpy as np

# Load protein data
sequence = sequences.read_fasta("protein.fasta")['protein1']
structure_data = structure.parse_pdb_file("protein.pdb")

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Sequence logo (if multiple sequences available)
# ax = prot_viz.plot_sequence_logo(aligned_sequences, ax=axes[0,0])

# Domain architecture
domains = [{'start': 10, 'end': 50, 'name': 'Domain1', 'type': 'domain'}]
prot_viz.plot_domain_architecture(domains, len(sequence), ax=axes[0,1])

# Contact map
coords = structure_data['coordinates']
contacts = structure.calculate_contact_map(coords)
prot_viz.plot_contact_map(contacts, ax=axes[1,0])

# Properties
prot_viz.plot_protein_properties(sequence, ax=axes[1,1])

plt.tight_layout()
plt.savefig("protein_analysis.png", dpi=300, bbox_inches='tight')
```

## Troubleshooting

### Common Issues

1. **Empty Plots**: Check input data dimensions and formats
2. **Memory Errors**: Reduce resolution or subsample large datasets
3. **Import Errors**: Ensure optional dependencies are installed for advanced features
4. **Color Issues**: Use `prot_viz.style.reset_style()` to reset matplotlib defaults

### Data Format Requirements

- **Sequences**: List of strings, all same length for logos
- **Coordinates**: Nx3 numpy arrays for 3D structures
- **Domains**: List of dicts with 'start', 'end', 'name', 'type' keys
- **Angles**: Numpy arrays in degrees for Ramachandran plots

## Related Documentation

- **[Protein Analysis](../protein/)**: Core protein analysis functions
- **[Structure Analysis](../protein/index.md)**: 3D structure algorithms
- **[Sequence Analysis](../protein/proteomes.md)**: Sequence processing utilities
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive protein visualization capabilities integrated with METAINFORMANT's protein analysis workflows.

