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
    radius=200.0,             # interaction neighbourhood, in coordinate units (k nearest neighbours if omitted)
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

## Cell-type level scoring

`compute_ligand_receptor_interactions` in `metainformant.spatial.communication.cell_communication` scores every ligand–receptor pair between every pair of cell types, without using the coordinates:

```python-snippet
from metainformant.spatial.communication.cell_communication import compute_ligand_receptor_interactions

result = compute_ligand_receptor_interactions(
    expression,                  # cells x genes matrix (scipy sparse is fine)
    cell_types,                  # one label per cell
    lr_database=None,            # default: default_lr_database()
    gene_names=gene_names,       # keyword-only; one name per column of expression
    seed=42,                     # keyword-only; None (the default) is not reproducible
)
significant = [ix for ix in result['interactions'] if ix['p_value'] < 0.05]
# each interaction: ligand, receptor, source_type, target_type, score, p_value
print(result['n_significant'], result['summary'])

# A dataset (or a DataFrame) carries its own gene names, so gene_names can be left out
result = compute_ligand_receptor_interactions(xenium_dataset, cell_types, seed=42)
```

The score of a (source, target) pair is the mean ligand expression in the source type times the mean receptor expression in the target type, divided by the product of the two genes' means over all cells. The p-value is the fraction of 100 label permutations (not configurable) whose score is at least as large.

**Seeds.** The permutations come from a local generator and numpy's global random state is never touched. With `seed=None` every interaction draws from a fresh generator seeded from operating-system entropy, so the p-values differ between calls (the scores do not); an integer `seed` makes them reproducible.

**Genes are matched by name.** Ligand and receptor symbols are matched to the columns of `expression` by exact, case-sensitive name; nothing else maps a name to a column, so a gene the data do not contain is never scored.

- The names come from `gene_names` (one per column, in column order). When `expression` is a dataset object with `.expression` and `.gene_names` (`SpatialDataset`, `XeniumDataset`, `MERFISHDataset`) or a pandas DataFrame whose columns are the genes, the names are taken from it; explicit `gene_names` win over the ones an object carries.
- A database pair is scored only when both of its genes are among the names. The other pairs are skipped and counted in `result['summary']`: `n_pairs_in_database`, `n_pairs_scored` and `n_pairs_skipped` (`total_pairs_tested` is the number of interaction rows, that is scored pair × source type × target type, and `unique_ligands` and `unique_receptors` count the scored pairs).
- `ValueError` is raised when no gene names are available, when `gene_names` does not have one entry per column or `cell_types` one label per cell, when no pair of the database has both of its genes among the names, and when a scored gene names more than one column (duplicate column names: make them unique first, for example by collapsing the duplicate columns).
- The built-in database uses human HGNC symbols, so data with another convention (another species, Ensembl identifiers) needs an `lr_database` written in its own gene names. A numeric string is an ordinary name: `'12'` matches only a column called `'12'`.

The result feeds two helpers: `build_communication_network(result['interactions'], min_score=0.0)` sums the scores per source–target cell-type pair into an adjacency matrix, edge list and hub types, and `communication_pattern_analysis(result, n_patterns=5)` takes the whole result dictionary and decomposes the (cell-type pair × ligand–receptor pair) score matrix into patterns with non-negative matrix factorisation. `spatial_interaction_score(expression, coordinates, lr_pairs, max_distance=100.0)` is a coordinate-based variant over gene-index pairs (`{'ligand_idx': ..., 'receptor_idx': ...}`) with an exponential distance decay. See [CAPABILITIES.md](CAPABILITIES.md) for their signatures.

## Pitfalls

- **Radius too small** — a radius below the spot spacing finds no neighbours and misses true interactions. The radius is in the units of `coordinates` (pixels for a Visium dataset from `load_visium`), so measure the spot spacing in those units first.
- **Gene-name convention** — the built-in database holds human HGNC symbols and genes are matched by exact name, so mouse symbols (`Egf`) or Ensembl identifiers match nothing. Check `result['summary']['n_pairs_scored']` against `n_pairs_in_database`.
- **Batch effects** — integrate batches first; otherwise false positives from
  technical variation.
