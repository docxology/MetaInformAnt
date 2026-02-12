# Cell-Cell Communication

The communication module analyzes cell-cell signaling in spatial transcriptomics data using ligand-receptor interaction databases and spatial proximity information. It scores interactions between cell types, builds communication networks, and discovers communication patterns.

## Key Concepts

### Ligand-Receptor Interactions

Cell-cell communication occurs when a ligand expressed by one cell binds to a receptor on a neighboring cell. The module scores these interactions by computing the product of mean ligand expression in the source cell type and mean receptor expression in the target cell type, normalized by background expression.

### Spatial Proximity Constraint

Unlike bulk interaction analyses, spatial communication scoring incorporates physical proximity. Only interactions between spatially adjacent cell types are considered biologically plausible. Distance thresholds or spatial neighborhood graphs filter interactions.

### Statistical Significance

Significance is assessed by permutation testing: cell type labels are shuffled to build a null distribution, and p-values are computed for each ligand-receptor-source-target combination.

### Communication Networks

Interactions are assembled into directed networks where nodes are cell types and edges represent significant ligand-receptor signaling events. Edge weights reflect interaction strength.

### Communication Patterns

Non-negative matrix factorization (NMF) on the interaction matrix reveals communication patterns -- recurring signaling programs that span multiple ligand-receptor pairs and cell type pairs.

## Function Reference

### default_lr_database

```python
def default_lr_database() -> dict[str, list[dict[str, str]]]
```

Return the built-in ligand-receptor pair database. Returns a dictionary with key `"pairs"` containing a list of `{"ligand": ..., "receptor": ...}` entries covering major signaling pathways.

### compute_ligand_receptor_interactions

```python
def compute_ligand_receptor_interactions(
    expression: Any,
    cell_types: list[str],
    lr_database: dict[str, list[dict[str, str]]] | None = None,
) -> dict[str, Any]
```

Compute ligand-receptor interaction scores between all cell type pairs. Returns a dictionary with `interactions` (list of interaction dicts with `ligand`, `receptor`, `source_type`, `target_type`, `score`, `p_value`), `n_significant`, and `summary`.

### spatial_interaction_score

```python
def spatial_interaction_score(
    expression: Any,
    coordinates: Any,
    cell_types: list[str],
    ligand: str,
    receptor: str,
    distance_threshold: float = 100.0,
) -> dict[str, Any]
```

Score a single ligand-receptor interaction with spatial proximity constraint. Only considers cell pairs within the distance threshold. Returns `score`, `p_value`, `n_interacting_pairs`, and `source_target_types`.

### build_communication_network

```python
def build_communication_network(
    interactions: dict[str, Any],
    significance_threshold: float = 0.05,
) -> dict[str, Any]
```

Build a directed communication network from significant interactions. Returns `nodes` (cell types), `edges` (list of source, target, weight, ligand, receptor), `adjacency_matrix`, and network summary statistics.

### communication_pattern_analysis

```python
def communication_pattern_analysis(
    interactions: dict[str, Any],
    n_patterns: int = 5,
) -> dict[str, Any]
```

Discover communication patterns using NMF decomposition of the interaction matrix. Returns `patterns` (list of pattern dictionaries with contributing ligand-receptor pairs and cell types), `pattern_scores`, and reconstruction error.

## Usage Examples

```python
from metainformant.spatial import io, communication

# Load spatial data with cell type annotations
dataset = io.load_visium("path/to/spaceranger_output/")
cell_types = [...]  # Cell type labels per spot

# Compute all ligand-receptor interactions
lr_db = communication.default_lr_database()
interactions = communication.compute_ligand_receptor_interactions(
    dataset.expression, cell_types, lr_database=lr_db
)
print(f"Significant interactions: {interactions['n_significant']}")

# Score a specific interaction with spatial constraint
score = communication.spatial_interaction_score(
    dataset.expression, dataset.coordinates, cell_types,
    ligand="TGFB1", receptor="TGFBR1", distance_threshold=150.0,
)
print(f"Score: {score['score']:.4f} (p={score['p_value']:.4e})")

# Build communication network
network = communication.build_communication_network(
    interactions, significance_threshold=0.05
)
print(f"Network: {len(network['nodes'])} cell types, {len(network['edges'])} edges")

# Discover communication patterns
patterns = communication.communication_pattern_analysis(interactions, n_patterns=5)
for i, pattern in enumerate(patterns["patterns"]):
    print(f"Pattern {i+1}: {len(pattern['lr_pairs'])} L-R pairs")
```

## Configuration

- **Environment prefix**: `SPATIAL_`
- **Optional dependencies**: numpy, scipy (cdist, mannwhitneyu)
- Built-in LR database covers major signaling pathways (TGFB, WNT, NOTCH, chemokines, etc.)
- Permutation test uses 1000 permutations by default
- Distance threshold should be set based on platform resolution (100um for Visium, 10-50um for subcellular platforms)

## Related Modules

- `spatial.analysis.neighborhood` -- Neighborhood enrichment and niche detection
- `spatial.deconvolution` -- Cell type fractions feed into communication analysis
- `spatial.visualization` -- `plot_interaction_heatmap`, `plot_neighborhood_graph`
- `spatial.io` -- Data loading for spatial datasets
