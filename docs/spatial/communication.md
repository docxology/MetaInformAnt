# Cell–Cell Communication

Infer signalling activity from spatial proximity of ligand‑expressing and
receptor‑expressing spots / cells.

## Core functions

```python-snippet
from metainformant.spatial.analysis.communication import (
    ligand_receptor_score,           # one pair
    spatial_ligand_receptor_pairs,   # many pairs, multiple testing
    permutation_significance,        # Monte Carlo p‑value
)
```

## Single pair scoring

Given a spatial AnnData (spots / cells) and a ligand and receptor gene:

```python
score = ligand_receptor_score(
    adata,
    ligand='EGF',
    receptor='EGFR',
    radius=200.0,          # µm, interaction neighbourhood
    method='exp_mean',     # mean expression of receptor among ligand‑expressing
)
print(f'Score = {score.value:.3f} (p = {score.pvalue:.3g})')
```

Methods:

| name | meaning |
|------|---------|
| `exp_mean` | average receptor expression conditional on ligand > 0 |
| `exp_fraction` | what fraction of ligand (+) neighbours express receptor |
| `coexpression` | Pearson correlation between ligand and receptor across spots |

## Database‑driven pair scan

Built‑in LR database (CellChat / CellPhoneDB v2) accessed via `load_lr_database()`:

```python-snippet
from metainformant.spatial.analysis.communication import load_lr_database, spatial_ligand_receptor_pairs
lr_db = load_lr_database('human')   # ~2 300 pairs

results = spatial_ligand_receptor_pairs(
    adata,
    lr_db=lr_db,
    radius=150.0,
    n_permutations=999,
    n_jobs=4,
)
# DataFrame: ligand, receptor, score, pvalue, qvalue, cell_type_pair
sig = results[results['qvalue'] < 0.05]
print(f'Significant pairs: {len(sig)} / {len(results)}')
```

Multiple‑testing correction uses Benjamini–Hochberg within each ligand or within
each cell‑type pair; `qvalue` column is FDR.

## Visualisation

```python-snippet
from metainformant.spatial.visualization import lr_network_graph
lr_network_graph(adata, results, min_score=0.3, edge_weight_scale=2.0)
```

Creates a bipartite graph (ligands left, receptors right) with edge thickness ∝
score, colour by significance.

## Custom LR pair

```python-snippet
from metainformant.spatial.analysis.communication import register_custom_lr
register_custom_lr(
    ligand='MYOKINE',
    receptor='MUSCLE_RECEPTOR',
    pathway='muscle_growth',
)
# now appears in future scans
```

## Pitfalls

- **Radius too small** — misses true interactions; try 100–300 µm for Visium.
- **Expression threshold** — by default ligand presence requires >0 counts; use
  `ligand_threshold=1` to require at least 1 UMI.
- **Batch effects** — integrate batches first; otherwise false positives from
  technical variation.
- **Cell‑type proportion bias** — abundant cell types dominate; normalise by cell‑type
  frequency using `normalize_by_celltype_freq=True`.
