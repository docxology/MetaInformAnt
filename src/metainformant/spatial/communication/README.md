# Communication

Cell-cell communication analysis for spatial transcriptomics, computing spatially-informed ligand-receptor interaction scores, communication networks, and communication pattern discovery.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports cell_communication submodule |
| `cell_communication.py` | L-R scoring, communication networks, NMF pattern analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `cell_communication.compute_ligand_receptor_interactions()` | Score L-R interactions between cell types |
| `cell_communication.spatial_interaction_score()` | Compute spatially-weighted interaction scores |
| `cell_communication.build_communication_network()` | Build cell-cell communication network graph |
| `cell_communication.default_lr_database()` | Get built-in ligand-receptor pair database |
| `cell_communication.communication_pattern_analysis()` | Discover communication patterns via NMF |

## Usage

```python
from metainformant.spatial.communication import cell_communication

interactions = cell_communication.compute_ligand_receptor_interactions(
    expression, cell_types
)
scores = cell_communication.spatial_interaction_score(
    expression, coordinates, cell_types
)
network = cell_communication.build_communication_network(interactions)
```
