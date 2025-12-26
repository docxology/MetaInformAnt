# Networks: Graph Analysis

The graph analysis module provides fundamental graph operations, metrics, and utilities for biological network analysis.

## Network Creation

### create_network()

Create biological networks from various data sources:

```python
from metainformant.networks import create_network, add_edges_from_interactions, add_edges_from_correlation

# Create from edge list
nodes = ["gene1", "gene2", "gene3"]
network = create_network(nodes, directed=False)

# Add edges from interaction list
interactions = [
    ("gene1", "gene2", 0.8),
    ("gene2", "gene3", 0.6),
    ("gene3", "gene1", 0.9)
]

add_edges_from_interactions(network, interactions)

# Create from correlation matrix
correlation_matrix = pd.DataFrame(...)
correlation_network = create_network(correlation_matrix.index.tolist())

# Add correlation-based edges
add_edges_from_correlation(
    correlation_network,
    correlation_matrix.values,
    correlation_matrix.index.tolist(),
    threshold=0.7           # Minimum correlation
)
```

### add_edges_from_correlation()

Add edges based on correlation between node features:

```python
from metainformant.networks import add_edges_from_correlation

# Add correlation-based edges
add_edges_from_correlation(
    network,
    node_features,           # Feature matrix (nodes × features)
    threshold=0.6,          # Minimum correlation
    method="pearson",       # Correlation method
    max_edges=1000          # Limit number of edges
)
```

### add_edges_from_interactions()

Add edges from known biological interactions:

```python
from metainformant.networks import add_edges_from_interactions

# Add known interactions
interactions = [
    {"source": "TF1", "target": "gene1", "type": "regulation", "confidence": 0.9},
    {"source": "gene1", "target": "gene2", "type": "interaction", "confidence": 0.7}
]

add_edges_from_interactions(network, interactions, confidence_threshold=0.5)
```

## Network Metrics

### network_metrics()

Calculate comprehensive network statistics:

```python
from metainformant.networks import network_metrics

# Get basic network metrics
metrics = network_metrics(network)

print(f"Nodes: {metrics['n_nodes']}")
print(f"Edges: {metrics['n_edges']}")
print(f"Average degree: {metrics['avg_degree']:.2f}")
print(f"Density: {metrics['density']:.3f}")
print(f"Clustering coefficient: {metrics['clustering_coeff']:.3f}")
```

**Available metrics:**
- Number of nodes and edges
- Average degree and degree distribution
- Network density
- Clustering coefficient
- Connected components
- Diameter (if connected)
- Average shortest path length

### centrality_measures()

Calculate node centrality measures:

```python
from metainformant.networks import centrality_measures

# Calculate all centrality measures
centrality = centrality_measures(network)

print("Top 5 nodes by degree centrality:")
for node, score in sorted(centrality['degree'].items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {node}: {score:.3f}")

print("Top 5 nodes by betweenness centrality:")
for node, score in sorted(centrality['betweenness'].items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {node}: {score:.3f}")
```

**Available measures:**
- **Degree centrality**: Number of connections
- **Betweenness centrality**: Control over information flow
- **Closeness centrality**: Distance to other nodes
- **Eigenvector centrality**: Influence in the network
- **PageRank**: Importance based on link structure

## Shortest Paths

### shortest_paths()

Find shortest paths between nodes:

```python
from metainformant.networks import shortest_paths

# Find shortest path between two nodes
path = shortest_paths(network, source="gene1", target="gene5")
print(f"Shortest path: {' -> '.join(path)}")

# Find shortest paths from one node to all others
all_paths = shortest_paths(network, source="gene1")

# Find shortest path lengths between all pairs
path_lengths = shortest_paths(network, source=None, target=None)
```

**Parameters:**
- `source`: Source node (or None for all-pairs)
- `target`: Target node (or None for single-source)
- `weight`: Edge weight attribute for weighted shortest paths

## Network Manipulation

### Network Filtering

Filter networks based on various criteria:

```python
# Filter by edge weight/confidence
filtered_network = network.copy()
edges_to_remove = [edge for edge, data in filtered_network.edges(data=True) if data.get('weight', 0) < 0.5]
filtered_network.remove_edges_from(edges_to_remove)

# Filter by degree
high_degree_nodes = [node for node, degree in filtered_network.degree() if degree > 5]
subgraph = filtered_network.subgraph(high_degree_nodes)
```

### Network Components

Work with connected components:

```python
# Find connected components
components = list(nx.connected_components(network))
print(f"Number of components: {len(components)}")

# Largest component
largest_component = max(components, key=len)
largest_subgraph = network.subgraph(largest_component)
print(f"Largest component size: {len(largest_component)} nodes")
```

## Integration Examples

### With Expression Data

```python
from metainformant.networks import create_network, network_metrics
import pandas as pd

# Create correlation network from expression data
expression_data = pd.read_csv("expression_matrix.csv", index_col=0)

# Calculate correlation matrix
correlation_matrix = expression_data.T.corr()

# Create network
correlation_network = create_network(
    correlation_matrix,
    method="correlation",
    threshold=0.8,
    directed=False
)

# Analyze network properties
metrics = network_metrics(correlation_network)
print(f"Correlated genes form {metrics['n_components']} modules")
```

### With Genomic Data

```python
from metainformant.networks import create_network, centrality_measures

# Create interaction network from genomic coordinates
genomic_interactions = [
    ("chr1:1000-2000", "chr1:3000-4000", 0.9),
    ("chr2:5000-6000", "chr2:7000-8000", 0.7)
]

genomic_network = create_network(genomic_interactions, directed=False)

# Find hub regions
centrality = centrality_measures(genomic_network)
hub_regions = [node for node, score in centrality['degree'].items() if score > centrality['degree'].mean() + centrality['degree'].std()]
```

## Performance Considerations

### Large Network Optimization

```python
# Use sparse representations for large networks
import networkx as nx

# Convert to sparse format for memory efficiency
if network.number_of_nodes() > 10000:
    # Use adjacency matrix representation
    adjacency_matrix = nx.adjacency_matrix(network)
    # Process with sparse matrix operations
```

### Parallel Processing

```python
# Parallel shortest path calculation
from concurrent.futures import ProcessPoolExecutor

def calculate_path_length(source_target):
    source, target = source_target
    try:
        return nx.shortest_path_length(network, source, target)
    except nx.NetworkXNoPath:
        return float('inf')

# Calculate all pairwise distances in parallel
node_pairs = [(s, t) for s in network.nodes() for t in network.nodes() if s != t]

with ProcessPoolExecutor(max_workers=4) as executor:
    distances = list(executor.map(calculate_path_length, node_pairs))
```

## Testing

Graph analysis functionality is tested in `tests/test_networks_graph.py`:

```bash
# Run graph analysis tests
uv run pytest tests/test_networks_graph.py -v

# Test specific functions
uv run pytest tests/test_networks_graph.py::test_create_network -v
uv run pytest tests/test_networks_graph.py::test_centrality_measures -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection algorithms
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
- [Regulatory Networks](./regulatory.md): Gene regulatory network analysis
