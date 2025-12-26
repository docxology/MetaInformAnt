# Networks: Gene Regulatory Network Analysis

The regulatory networks module provides tools for inferring and analyzing gene regulatory networks from expression data and other biological information.

## Network Construction

### GeneRegulatoryNetwork Class

Container for gene regulatory networks with regulatory relationships:

```python
from metainformant.networks import GeneRegulatoryNetwork

# Create regulatory network
reg_network = GeneRegulatoryNetwork()

# Add regulatory relationships
regulations = [
    {"regulator": "TF1", "target": "gene1", "effect": 1.0, "confidence": 0.9, "type": "activation"},
    {"regulator": "TF2", "target": "gene2", "effect": -0.8, "confidence": 0.7, "type": "repression"},
    {"regulator": "TF1", "target": "gene3", "effect": 0.5, "confidence": 0.6, "type": "activation"}
]

reg_network.add_regulations(regulations)
print(f"Network has {reg_network.n_regulators} regulators and {reg_network.n_targets} targets")
```

### infer_grn()

Infer gene regulatory networks from expression data:

```python
from metainformant.networks import infer_grn

# Infer network from expression data
expression_data = pd.read_csv("expression_matrix.csv", index_col=0)
regulators = ["TF1", "TF2", "TF3", "TF4"]  # Known transcription factors

reg_network = infer_grn(
    expression_data,
    regulators=regulators,
    method="correlation",      # Inference method
    threshold=0.7,            # Minimum correlation/regression coefficient
    max_targets=50            # Maximum targets per regulator
)

print(f"Inferred {reg_network.n_regulations} regulatory relationships")
```

**Available inference methods:**
- `'correlation'`: Pearson/Spearman correlation
- `'regression'`: Linear regression coefficients
- `'mutual_info'`: Mutual information
- `'granger'`: Granger causality (requires time-series data)
- `'ensemble'`: Ensemble of multiple methods

### add_transcription_factor()

Add transcription factor metadata:

```python
from metainformant.networks import add_transcription_factor

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

reg_network.add_transcription_factors(tf_info)

# Query TF information
tf_data = reg_network.get_tf_info("TF1")
print(f"TF1 family: {tf_data['family']}")
print(f"TF1 binding motif: {tf_data['binding_motif']}")
```

## Network Analysis

### get_targets()

Find target genes for a regulator:

```python
from metainformant.networks import get_targets

# Get all targets of a regulator
targets = get_targets(reg_network, "TF1")
print(f"TF1 regulates {len(targets)} genes")

# Get targets with regulation strength
targets_with_effect = get_targets(reg_network, "TF1", include_effect=True)
for target, effect in targets_with_effect.items():
    regulation_type = "activation" if effect > 0 else "repression"
    print(f"  {target}: {effect:.2f} ({regulation_type})")
```

### get_regulators()

Find regulators for a target gene:

```python
from metainformant.networks import get_regulators

# Get all regulators of a gene
regulators = get_regulators(reg_network, "gene1")
print(f"gene1 is regulated by {len(regulators)} transcription factors")

# Get regulators with confidence scores
regulators_with_confidence = get_regulators(reg_network, "gene1", include_confidence=True)
for regulator, confidence in regulators_with_confidence.items():
    print(f"  {regulator}: confidence = {confidence:.2f}")
```

### filter_by_confidence()

Filter regulatory relationships by confidence:

```python
from metainformant.networks import filter_by_confidence

# Filter high-confidence regulations
high_conf_network = filter_by_confidence(
    reg_network,
    min_confidence=0.8,
    regulation_types=["activation"]  # Only activation relationships
)

# Filter by effect size
strong_regulations = filter_by_confidence(
    reg_network,
    min_effect_size=0.7
)
```

## Regulatory Motifs

### regulatory_motifs()

Identify regulatory motifs and patterns:

```python
from metainformant.networks import regulatory_motifs

# Find regulatory motifs
motifs = regulatory_motifs(
    reg_network,
    motif_type="feed_forward",  # Type of motif to find
    min_confidence=0.7
)

print(f"Found {len(motifs)} regulatory motifs")

for motif in motifs[:5]:  # Show first 5 motifs
    print(f"Motif: {motif['type']}")
    print(f"  Nodes: {motif['nodes']}")
    print(f"  Confidence: {motif['confidence']:.2f}")
```

**Available motif types:**
- `'feed_forward'`: Feed-forward loops
- `'feedback'`: Feedback loops
- `'cascade'`: Regulatory cascades
- `'hub'`: Hub regulators with many targets
- `'module'`: Co-regulated gene modules

## Integration Examples

### With Expression Data

```python
from metainformant.rna import workflow
from metainformant.networks import infer_grn, regulatory_motifs

# Load expression data
expression_data = workflow.extract_expression_patterns(rna_data)

# Infer regulatory network
reg_network = infer_grn(
    expression_data,
    method="correlation",
    threshold=0.7
)

# Find regulatory motifs
motifs = regulatory_motifs(reg_network, motif_type="feed_forward")

# Analyze motif expression patterns
for motif in motifs:
    motif_genes = motif['nodes']
    motif_expression = expression_data[expression_data.index.isin(motif_genes)]
    # Analyze coordinated expression in motif
```

### With Sequence Data

```python
from metainformant.dna import motifs as dna_motifs
from metainformant.networks import add_transcription_factor

# Load DNA motifs
binding_sites = dna_motifs.find_motif_positions(sequences, "CANNTG")

# Add TF binding information to regulatory network
for tf, motif in binding_sites.items():
    tf_info = {
        tf: {
            "binding_motif": "CANNTG",
            "binding_sites": motif,
            "motif_strength": 0.9
        }
    }
    reg_network.add_transcription_factors(tf_info)
```

## Visualization

### Network Visualization

```python
import matplotlib.pyplot as plt
import networkx as nx

# Create visualization
fig, ax = plt.subplots(figsize=(14, 10))

# Color nodes by type (regulators vs targets)
node_colors = []
node_sizes = []

for node in reg_network.nodes():
    if reg_network.is_regulator(node):
        node_colors.append('red')    # Regulators in red
        node_sizes.append(300)
    else:
        node_colors.append('blue')   # Targets in blue
        node_sizes.append(200)

# Draw network
pos = nx.spring_layout(reg_network, k=2, iterations=50)
nx.draw(
    reg_network,
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
for edge in reg_network.edges():
    regulation_type = reg_network.get_regulation_type(edge[0], edge[1])
    edge_labels[edge] = regulation_type[:3]  # First 3 letters

nx.draw_networkx_edge_labels(reg_network, pos, edge_labels=edge_labels, font_size=6)

plt.title('Gene Regulatory Network')
plt.show()
```

### Regulation Strength Visualization

```python
# Visualize regulation strength
fig, ax = plt.subplots(figsize=(12, 8))

# Get regulation strengths
edge_weights = [reg_network.get_effect(u, v) for u, v in reg_network.edges()]

# Color edges by regulation strength
edge_colors = ['red' if w > 0 else 'blue' for w in edge_weights]

nx.draw(
    reg_network,
    pos=pos,
    ax=ax,
    node_color=node_colors,
    node_size=node_sizes,
    edge_color=edge_colors,
    edge_cmap=plt.cm.RdBu,
    width=[abs(w)*3 for w in edge_weights],  # Edge width by strength
    with_labels=True,
    font_size=8
)

plt.title('Regulatory Network (Red: Activation, Blue: Repression)')
plt.show()
```

## Testing

Regulatory network functionality is tested in `tests/test_networks_regulatory.py`:

```bash
# Run regulatory network tests
uv run pytest tests/test_networks_regulatory.py -v

# Test specific functions
uv run pytest tests/test_networks_regulatory.py::test_infer_grn -v
uv run pytest tests/test_networks_regulatory.py::test_regulatory_motifs -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection in regulatory networks
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
