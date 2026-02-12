# Networks: Community Detection

The community detection module provides algorithms for identifying functional modules and communities within biological networks.

## Community Detection Algorithms

### detect_communities()

Detect communities in biological networks using various algorithms:

```python
from metainformant.networks import detect_communities

# Louvain algorithm (default)
communities = detect_communities(
    network,
    method='louvain',
    resolution=1.0,        # Higher = more communities
    random_state=42
)

# Leiden algorithm
communities = detect_communities(
    network,
    method='leiden',
    resolution=1.0,
    n_iterations=2
)

# Hierarchical clustering
communities = detect_communities(
    network,
    method='hierarchical',
    n_clusters=10,
    linkage='ward'
)
```

**Supported Methods:**
- `'louvain'`: Louvain algorithm (fast, good for large networks)
- `'leiden'`: Leiden algorithm (improved Louvain)
- `'hierarchical'`: Agglomerative hierarchical clustering
- `'spectral'`: Spectral clustering

### modularity()

Calculate modularity score for community assignments:

```python
from metainformant.networks import modularity

# Calculate modularity
mod_score = modularity(network, communities)
print(f"Modularity: {mod_score:.3f}")

# Compare different community assignments
communities1 = detect_communities(network, resolution=0.5)
communities2 = detect_communities(network, resolution=1.0)

mod1 = modularity(network, communities1)
mod2 = modularity(network, communities2)
print(f"Resolution 0.5: {mod1:.3f}")
print(f"Resolution 1.0: {mod2:.3f}")
```

### community_metrics()

Calculate comprehensive community metrics:

```python
from metainformant.networks import community_metrics

# Get detailed community statistics
metrics = community_metrics(network, communities)

print(f"Number of communities: {metrics['n_communities']}")
print(f"Community size range: {metrics['size_range']}")
print(f"Largest community: {metrics['largest_community']}")
print(f"Modularity: {metrics['modularity']:.3f}")
```

**Metrics include:**
- Number of communities
- Community size distribution
- Modularity score
- Community density
- Internal vs external connectivity

## Functional Analysis

### Community Enrichment

Analyze functional enrichment within communities:

```python
# Get community assignments
communities = detect_communities(network, method='louvain')

# Functional enrichment analysis
for community_id in set(communities.values()):
    community_nodes = [node for node, comm in communities.items() if comm == community_id]

    # Analyze enrichment for this community
    enrichment = analyze_community_enrichment(community_nodes, annotation_db)
    print(f"Community {community_id}: {len(community_nodes)} nodes")
    print(f"Top enriched terms: {enrichment[:5]}")
```

## Network Visualization

### Community Visualization

Visualize communities on network layouts:

```python
import matplotlib.pyplot as plt

# Plot network colored by communities
fig, ax = plt.subplots(figsize=(12, 8))

# Color nodes by community
colors = [communities.get(node, 0) for node in network.nodes()]
nx.draw(
    network,
    ax=ax,
    node_color=colors,
    node_size=100,
    with_labels=False,
    cmap='tab20',
    alpha=0.7
)

plt.title('Network Communities')
plt.show()
```

## Performance Optimization

### Large Network Processing

For networks with >10k nodes:

```python
# Use efficient algorithms for large networks
communities = detect_communities(
    network,
    method='louvain',
    resolution=1.0,
    use_multiprocessing=True,
    max_workers=4
)
```

## Integration Examples

### With Protein Interaction Networks

```python
from metainformant.networks import ppi, community
from metainformant.ontology.core import go

# Load PPI network
protein_network = ppi.load_string_network()

# Detect functional modules
modules = community.detect_communities(protein_network, method='leiden')

# Functional enrichment of modules
for module_id, proteins in modules.items():
    enriched_terms = go.enrich_proteins(proteins)
    print(f"Module {module_id}: {enriched_terms[:3]}")
```

### With Gene Expression Data

```python
from metainformant.networks import regulatory, community

# Build regulatory network
reg_network = regulatory.infer_network(expression_data)

# Find regulatory modules
reg_modules = community.detect_communities(reg_network)

# Analyze module expression patterns
for module_id, genes in reg_modules.items():
    module_expression = expression_data[expression_data.index.isin(genes)]
    # Analyze coordinated expression patterns
```

## Testing

Community detection functionality is tested in `tests/test_networks_community.py`:

```bash
# Run community detection tests
uv run pytest tests/test_networks_community.py -v

# Test specific algorithms
uv run pytest tests/test_networks_community.py::test_detect_communities_louvain -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
- [Regulatory Networks](./regulatory.md): Gene regulatory network analysis
