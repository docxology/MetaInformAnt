# Cell–Cell Communication

Infer signalling activity from spatial proximity of ligand‑expressing and
receptor‑expressing spots / cells.

## Core functions

```python-snippet
from metainformant.spatial.analysis.neighborhood import ligand_receptor_spatial
from metainformant.spatial.communication.cell_communication import default_lr_database
```

- `ligand_receptor_spatial()` scores ligand–receptor pairs from a spots × genes
  matrix, the spot coordinates and a list of `(ligand, receptor)` gene tuples.
- `default_lr_database()` returns the built‑in list of ligand–receptor pairs.

## Single pair scoring

Given a spots × genes expression matrix, spot coordinates and a ligand and receptor gene:

```python
res = ligand_receptor_spatial(
    expression,
    [('EGF', 'EGFR')],        # (ligand, receptor) pairs
    coordinates,
    gene_names=gene_names,
    radius=200.0,             # µm, interaction neighbourhood (k nearest neighbours if omitted)
)
score = res['scores'][('EGF', 'EGFR')]
print(f'Score = {score:.3f}')
```

The score of a pair is the mean over spots of the ligand expression in the spot
times the mean receptor expression of its neighbours; `res['per_spot_scores']`
holds the per‑spot values. No p‑value is computed.

## Database‑driven pair scan

The built‑in LR database is a curated subset of 95 ligand–receptor pairs (chemokine,
growth factor, Wnt, Notch, Hedgehog, TGF‑β, interleukin, adhesion and other
signalling), returned by `default_lr_database()`:

```python-snippet
from metainformant.spatial.analysis.neighborhood import ligand_receptor_spatial
from metainformant.spatial.communication.cell_communication import default_lr_database

lr_db = default_lr_database()   # {"pairs": [{"ligand": ..., "receptor": ...}, ...]}
lr_pairs = [(p['ligand'], p['receptor']) for p in lr_db['pairs']]

results = ligand_receptor_spatial(
    expression,
    lr_pairs,
    coordinates,
    gene_names=gene_names,
    radius=150.0,
)
# results['scores']: {(ligand, receptor): score}; pairs with a gene missing from gene_names are skipped
ranked = sorted(results['scores'].items(), key=lambda kv: kv[1], reverse=True)
print(f"Pairs scored: {results['n_pairs_tested']} / {len(lr_pairs)}")
```

## Visualisation

Per‑spot scores can be mapped back onto the tissue:

```python-snippet
from metainformant.spatial.visualization.plots import plot_spatial_scatter

top_pair = ranked[0][0]
plot_spatial_scatter(
    coordinates,
    results['per_spot_scores'][top_pair],
    'lr_scores.png',
    title=' → '.join(top_pair),
    colorbar_label='LR score',
)
```

## Custom LR pair

`ligand_receptor_spatial()` takes plain `(ligand, receptor)` tuples, so a custom pair
is simply appended to the list:

```python-snippet
lr_pairs.append(('MYOKINE', 'MUSCLE_RECEPTOR'))
results = ligand_receptor_spatial(expression, lr_pairs, coordinates, gene_names=gene_names, radius=150.0)
```

## Pitfalls

- **Radius too small** — misses true interactions; try 100–300 µm for Visium.
- **Batch effects** — integrate batches first; otherwise false positives from
  technical variation.
