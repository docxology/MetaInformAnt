# Networks: Gene Regulatory Network Analysis

The regulatory networks module provides tools for inferring and analyzing gene regulatory networks from expression data and other biological information.

The classes and functions below live in `metainformant.networks.interaction.regulatory_core` (the `GeneRegulatoryNetwork` class and graph-level helpers) and `metainformant.networks.interaction.regulatory_analysis` (inference, motifs, cascades, validation). The `metainformant.networks.interaction.regulatory` facade re-exports both.

## Network Construction

### GeneRegulatoryNetwork Class

Container for gene regulatory networks with regulatory relationships. Regulations are added one at a time with `add_regulation(tf, target, **attributes)`; the attributes are stored on the edge of the underlying directed `networkx` graph (`network.graph`):

```python
from metainformant.networks.interaction.regulatory_core import GeneRegulatoryNetwork

# Create regulatory network
reg_network = GeneRegulatoryNetwork(name="demo_grn")

# Add regulatory relationships
regulations = [
    ("TF1", "gene1", {"effect": 1.0, "confidence": 0.9, "type": "activation"}),
    ("TF2", "gene2", {"effect": -0.8, "confidence": 0.7, "type": "repression"}),
    ("TF1", "gene3", {"effect": 0.5, "confidence": 0.6, "type": "activation"}),
]
for regulator, target, attributes in regulations:
    reg_network.add_regulation(regulator, target, **attributes)

print(
    f"Network has {len(reg_network.transcription_factors)} regulators, "
    f"{len(reg_network.get_targets())} targets and {len(reg_network)} regulations"
)
```

`len(network)` is the number of regulations, `network.transcription_factors` lists the regulators, `network.genes` lists every gene and `network.regulations` yields `(tf, target, attributes)` tuples. Class methods build a network in one call: `GeneRegulatoryNetwork.from_interactions([(tf, target, confidence), ...])` and `GeneRegulatoryNetwork.from_tf_target_file(path)` (tab-separated `TF<TAB>target<TAB>confidence`).

### infer_grn()

Infer gene regulatory networks from expression data. `expression_data` is a samples x genes matrix (a pandas DataFrame or a NumPy array with one column per gene) and `gene_names` must have one entry per column:

```python
import numpy as np
import pandas as pd
from metainformant.networks.interaction.regulatory_analysis import infer_grn

# Expression matrix with samples as rows and genes as columns
rng = np.random.default_rng(0)
tf1, tf2 = rng.normal(size=40), rng.normal(size=40)
expression_data = pd.DataFrame({
    "TF1": tf1,
    "TF2": tf2,
    "gene1": 0.9 * tf1 + rng.normal(scale=0.3, size=40),
    "gene2": -0.8 * tf2 + rng.normal(scale=0.3, size=40),
    "gene3": rng.normal(size=40),
})
regulators = ["TF1", "TF2"]  # Known transcription factors

inferred_network = infer_grn(
    expression_data,
    list(expression_data.columns),  # gene_names, one per column
    method="correlation",           # Inference method
    threshold=0.7,                  # Minimum absolute correlation
    tf_genes=regulators,            # Candidate regulators
)

print(f"Inferred {len(inferred_network)} regulatory relationships")
```

The result is a `GeneRegulatoryNetwork`. Each inferred edge stores the score as `weight` (the absolute correlation or the mutual information) and has `type="regulates"`, so inferred edges carry no sign. Without `tf_genes`, the `'correlation'` method treats the genes in the top quartile by variance as regulators and `'mutual_info'` treats every gene as one.

**Available inference methods** (any other value raises `ValueError`):
- `'correlation'`: absolute Pearson correlation between each regulator and every other gene
- `'mutual_info'`: mutual information (`sklearn.feature_selection.mutual_info_regression`)

`infer_grn` has no Granger-causality method: a Granger test needs time-ordered observations, which a samples x genes matrix does not guarantee, so `method="granger"` raises `ValueError`. For time series, use `information_flow_network(series, method="granger_causality")` from `metainformant.information.network_info.information_flow`.

`'mutual_info'` is the one randomised method, because scikit-learn adds a tiny random jitter to continuous inputs when it estimates mutual information (the jitter breaks ties between equal values, as in count data). `infer_grn` and `infer_regulatory_network_from_expression` take a `random_state` (an integer seed or a `numpy.random.RandomState`) and share one generator across all the `mutual_info_regression` calls of a run, which are made pair by pair (regulators in the order given, and for each regulator the targets in column or index order). The same seed gives the same weights and edges, and the global NumPy random state is neither seeded nor consumed; `None` draws fresh entropy. `'correlation'` is deterministic and ignores `random_state`.

Other routines take a genes x samples matrix and return a plain dictionary with an `edges` list of `{"source", "target", "weight"}` records: `infer_grn_correlation` (`'pearson'` or `'spearman'`), `infer_grn_mutual_info` (ARACNE-style pruning) and `infer_grn_regression` (`'lasso'` or `'ridge'`) in `metainformant.networks.regulatory.grn_inference`.

`infer_regulatory_network_from_expression(expression_df, regulators, method="correlation", threshold=0.5, random_state=None)` takes a genes x samples DataFrame and returns a `networkx.DiGraph` (correlation edges keep the sign of the correlation). It supports the same two methods. `method="grnboost2"` is not implemented, because GRNBoost2 needs the `arboreto` package, and raises `NotImplementedError`; any other unknown method raises `ValueError`.

### add_transcription_factor()

Add transcription factor metadata (a `GeneRegulatoryNetwork` method):

```python
# Add TF information
tf_info = {
    "TF1": {
        "family": "bHLH",
        "binding_motif": "CANNTG",
        "expression_level": "high",
        "activity_status": "active"
    },
    "TF2": {
        "family": "Homeodomain",
        "binding_motif": "TAAT",
        "expression_level": "medium",
        "activity_status": "repressed"
    }
}

for tf, info in tf_info.items():
    reg_network.add_transcription_factor(tf, **info)

# Query TF information
tf_data = reg_network.gene_metadata["TF1"]
print(f"TF1 family: {tf_data['family']}")
print(f"TF1 binding motif: {tf_data['binding_motif']}")
```

## Network Analysis

### get_targets()

Find target genes for a regulator (a `GeneRegulatoryNetwork` method). Without an argument it returns every target gene in the network:

```python
# Get all targets of a regulator
targets = reg_network.get_targets("TF1")
print(f"TF1 regulates {len(targets)} genes")

# Get targets with regulation strength (stored on the edges)
for target in targets:
    effect = reg_network.graph["TF1"][target].get("effect", 0.0)
    regulation_type = "activation" if effect > 0 else "repression"
    print(f"  {target}: {effect:.2f} ({regulation_type})")
```

### get_regulators()

Find regulators for a target gene (a `GeneRegulatoryNetwork` method):

```python
# Get all regulators of a gene
regulators = reg_network.get_regulators("gene1")
print(f"gene1 is regulated by {len(regulators)} transcription factors")

# Get regulators with confidence scores (stored on the edges)
for regulator in regulators:
    confidence = reg_network.graph[regulator]["gene1"].get("confidence", 1.0)
    print(f"  {regulator}: confidence = {confidence:.2f}")
```

### filter_by_confidence()

Filter regulatory relationships by confidence or type (`GeneRegulatoryNetwork` methods that return a new network):

```python
# Filter high-confidence regulations
high_conf_network = reg_network.filter_by_confidence(0.8)

# Keep only activation relationships (matches the edge's "type" attribute)
activation_network = reg_network.filter_by_regulation_type("activation")
```

## Regulatory Motifs

### regulatory_motifs()

Identify regulatory motifs and patterns. The function takes only the network:

```python
from metainformant.networks.interaction.regulatory_analysis import regulatory_motifs

# Let TF1 regulate TF2 and gene2, and TF2 regulate gene2: the three edges
# TF1 -> TF2, TF2 -> gene2 and TF1 -> gene2 close a feed-forward loop
reg_network.add_regulation("TF1", "TF2", effect=0.9, confidence=0.8, type="activation")
reg_network.add_regulation("TF1", "gene2", effect=0.6, confidence=0.7, type="activation")

# Find regulatory motifs
motifs = regulatory_motifs(reg_network)

print(f"Found {len(motifs)} regulatory motifs")  # 1: the feed-forward loop

for motif in motifs[:5]:  # Show first 5 motifs
    print(f"Motif: {motif['motif_type']}")
    print(f"  Genes: {motif['genes']}")
    print(f"  Confidence: {motif['confidence']:.2f}")
```

**Reported motif types:**
- `'feed_forward_loop'`: TF1 -> TF2, TF2 -> target and the direct TF1 -> target edge, where TF2 is itself a regulator (`genes` is `[TF1, TF2, target]`). All three edges are required: without the direct TF1 -> target edge the chain is a cascade, not a feed-forward loop, and is not reported
- `'auto_regulation'`: a regulator that regulates itself
- `'mutual_regulation'`: two regulators that regulate each other

The `confidence` of a motif is a fixed score per motif type (0.8, 1.0 and 0.9 respectively), not a value derived from the edge confidences. The `GeneRegulatoryNetwork.regulatory_motifs()` method applies the same feed-forward rule and returns a dictionary with one list of gene lists per motif type (`feed_forward`, `auto_regulation`, `mutual_regulation` and a `cascade` entry that stays empty; use `detect_regulatory_cascades` for cascades). Related graph-level functions: `detect_regulatory_cascades(reg_network, max_length=5, min_confidence=0.0)` returns paths of three or more genes, `identify_regulatory_hubs(reg_network.graph, direction="out", threshold=0.1)` returns the top fraction of nodes by out-degree as `(gene, degree)` tuples (`direction` can also be `"in"` or `"both"`), and `analyze_regulatory_motifs(reg_network.graph)` reports feed-forward loops, feedback loops and cascades of a plain `networkx` graph.

### Randomised analyses and seeds

Two analyses draw random numbers and take a `random_state`: the z-scores of `compute_network_motifs` and the noise of `analyze_regulatory_dynamics`. In both, the same seed reproduces the output, and the global `random` and NumPy random states are neither seeded nor consumed; `None` draws fresh entropy.

```python
from metainformant.networks.interaction.regulatory_core import analyze_regulatory_dynamics
from metainformant.networks.regulatory.grn_inference import compute_network_motifs

# Motif counts for an edge list of {"source", "target"} records, with z-scores against 100 randomised
# networks (the targets are shuffled). random_state is an integer seed or a random.Random
edges = [{"source": source, "target": target} for source, target in reg_network.graph.edges()]
report = compute_network_motifs(edges, motif_size=3, random_state=7)
print(report["motif_counts"], report["z_scores"], report["significant_motifs"])

# Noisy simulation of the regulatory graph. random_state is an integer seed or a numpy.random.Generator
simulation = analyze_regulatory_dynamics(reg_network.graph, {"TF1": 1.0}, time_steps=20, random_state=7)
print(simulation["dynamics"]["gene1"]["final_state"])
```

`compute_network_motifs` counts `feed_forward_loop` (A -> B, B -> C and the direct A -> C), `cascade` (A -> B -> C without A -> C) and `co_regulation` for `motif_size=3`, or `mutual_regulation` and `single_regulation` for `motif_size=2`. A `z_scores` entry is 0.0 when the randomised counts do not vary, and `significant_motifs` lists the motif types with an absolute z-score above 2.

## Integration Examples

### With Expression Data

```python-snippet
import pandas as pd
from metainformant.rna.analysis import expression
from metainformant.networks.interaction.regulatory_analysis import infer_grn, regulatory_motifs

# counts: raw count matrix with genes as rows and samples as columns
counts = pd.read_csv("counts.csv", index_col=0)

# Normalize, then transpose: infer_grn expects samples as rows and genes as columns
normalized = expression.normalize_counts(counts, method="log2cpm")
samples_by_genes = normalized.T

# Infer regulatory network
reg_network = infer_grn(
    samples_by_genes,
    list(samples_by_genes.columns),
    method="correlation",
    threshold=0.7
)

# Find regulatory motifs
motifs = regulatory_motifs(reg_network)

# Analyze motif expression patterns
for motif in motifs:
    motif_expression = normalized.loc[normalized.index.isin(motif['genes'])]
    # Analyze coordinated expression in motif
```

### With Sequence Data

```python-snippet
from metainformant.dna.sequence import motifs as dna_motifs

# Promoter sequences of the target genes (gene -> DNA string)
promoters = {"gene1": "TTCACGTGAA", "gene2": "GGGGCCCC", "gene3": "ACCAATTGGT"}

# Scan each promoter for the TF's binding motif (IUPAC code: N matches any base)
reg_network.add_transcription_factor("TF1", binding_motif="CANNTG")
for gene, promoter in promoters.items():
    positions = dna_motifs.find_motif_positions(promoter, "CANNTG")
    if positions:
        # Record the binding sites on the TF1 -> gene regulation
        reg_network.add_regulation("TF1", gene, type="regulates", binding_sites=positions)
```

## Visualization

### Network Visualization

```python
import matplotlib.pyplot as plt
import networkx as nx

# The underlying directed graph and the set of regulators
graph = reg_network.graph
regulators = set(reg_network.transcription_factors)

# Create visualization
fig, ax = plt.subplots(figsize=(14, 10))

# Color nodes by type (regulators vs targets)
node_colors = []
node_sizes = []

for node in graph.nodes():
    if node in regulators:
        node_colors.append('red')    # Regulators in red
        node_sizes.append(300)
    else:
        node_colors.append('blue')   # Targets in blue
        node_sizes.append(200)

# Draw network
pos = nx.spring_layout(graph, k=2, iterations=50)
nx.draw(
    graph,
    pos=pos,
    ax=ax,
    node_color=node_colors,
    node_size=node_sizes,
    edge_color='gray',
    alpha=0.7,
    with_labels=True,
    font_size=8
)

# Add edge labels for regulation type
edge_labels = {}
for source, target, data in graph.edges(data=True):
    regulation_type = data.get("type", "regulates")
    edge_labels[(source, target)] = regulation_type[:3]  # First 3 letters

nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_size=6)

plt.title('Gene Regulatory Network')
plt.show()
```

### Regulation Strength Visualization

```python
# Visualize regulation strength
fig, ax = plt.subplots(figsize=(12, 8))

# Get regulation strengths (the "effect" attribute set when each regulation was added)
edge_weights = [data.get("effect", 0.0) for _, _, data in graph.edges(data=True)]

# Color edges by regulation sign
edge_colors = ['red' if w > 0 else 'blue' for w in edge_weights]

nx.draw(
    graph,
    pos=pos,
    ax=ax,
    node_color=node_colors,
    node_size=node_sizes,
    edge_color=edge_colors,
    width=[abs(w)*3 for w in edge_weights],  # Edge width by strength
    with_labels=True,
    font_size=8
)

plt.title('Regulatory Network (Red: Activation, Blue: Repression)')
plt.show()
```

## Testing

Regulatory network functionality is tested in `tests/networks/test_networks_regulatory.py`:

```bash
# Run regulatory network tests
uv run pytest tests/networks/test_networks_regulatory.py -v

# Test specific classes and functions
uv run pytest tests/networks/test_networks_regulatory.py::TestInferGRN::test_infer_grn_basic -v
uv run pytest tests/networks/test_networks_regulatory.py::TestRegulatoryMotifs::test_regulatory_motifs_basic -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection in regulatory networks
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
