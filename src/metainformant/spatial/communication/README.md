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
| `cell_communication.compute_ligand_receptor_interactions()` | Score L-R interactions between cell types with a 100-permutation test (one-sided p-value `(b + 1) / 101`, never 0). Ligand and receptor symbols are matched to the columns of `expression` by exact, case-sensitive name, taken from `gene_names=` or from a dataset or DataFrame that carries its own names. A pair with a gene that is absent is skipped and counted in `summary` (`n_pairs_in_database`, `n_pairs_scored`, `n_pairs_skipped`); `ValueError` is raised for missing or wrong-length names, when no pair has both genes present, or when a scored gene names two columns. `seed=None` (the default) draws the permutations from a fresh generator seeded from operating-system entropy, so the p-values differ between calls; pass an integer for reproducible p-values |
| `cell_communication.spatial_interaction_score()` | Compute spatially-weighted interaction scores for gene-index pairs (`ligand_idx`, `receptor_idx`) |
| `cell_communication.build_communication_network()` | Build cell-cell communication network graph |
| `cell_communication.default_lr_database()` | Get built-in ligand-receptor pair database (95 pairs, human HGNC symbols) |
| `cell_communication.communication_pattern_analysis()` | Discover communication patterns via NMF |

The built-in database uses human HGNC symbols, so data with another gene-name convention needs an `lr_database` written in its own names. `analysis.neighborhood.ligand_receptor_spatial` scores `(ligand, receptor)` pairs on spots with their coordinates and also matches genes by name (`gene_names`).

## Usage

```python
from metainformant.spatial.communication import cell_communication

interactions = cell_communication.compute_ligand_receptor_interactions(
    expression, cell_types, gene_names=gene_names, seed=42
)
print(interactions["summary"]["n_pairs_scored"], "of", interactions["summary"]["n_pairs_in_database"])
network = cell_communication.build_communication_network(interactions["interactions"])
scores = cell_communication.spatial_interaction_score(
    expression, coordinates, [{"ligand_idx": 0, "receptor_idx": 1}], max_distance=100.0
)
```
