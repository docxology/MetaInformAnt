# Networks: Community Detection

The community detection module provides algorithms for identifying functional modules and communities within biological networks.

The examples below use a small toy network of three clusters; `detect_communities` accepts a `networkx` graph, a `BiologicalNetwork`, or any wrapper (`ProteinNetwork`, `GeneRegulatoryNetwork`) that exposes its graph as `.graph`.

```python
import networkx as nx

# Three 5-node clusters joined in a ring
network = nx.relabel_nodes(nx.connected_caveman_graph(3, 5), lambda i: f"gene{i}")
```

## Community Detection Algorithms

### detect_communities()

Detect communities in biological networks using various algorithms. The result is a dictionary mapping each node to an integer community ID:

```python
from metainformant.networks.analysis.community import detect_communities

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
    n_iterations=2,
    random_state=42
)

# Girvan-Newman edge removal, stopped at 3 communities
communities = detect_communities(
    network,
    method='girvan_newman',
    n_communities=3
)

print(f"{len(set(communities.values()))} communities")
```

**Supported Methods** (any other name raises `ValueError`):
- `'louvain'`: Louvain algorithm (needs the optional `python-louvain` package)
- `'leiden'`: Leiden algorithm (needs the optional `leidenalg` and `python-igraph` packages)
- `'greedy'`: NetworkX greedy modularity optimization
- `'girvan_newman'`: Girvan-Newman edge-betweenness removal (`n_communities=` sets the target count)
- `'label_propagation'`: Label propagation
- `'asyn_lpa'`: Asynchronous label propagation
- `'fluid'`: Fluid communities (requires `k=`, the number of communities)

Leiden uses RB configuration modularity so `resolution` controls the partition.
`random_state` seeds the Leiden optimiser directly for reproducible runs.

Louvain hands `random_state` to python-louvain's own `random_state` argument, which builds a random generator of its own: the same seed gives the same partition, and the global `random` and `numpy.random` states are neither seeded nor consumed. With a seed, `randomize` is not passed on (python-louvain treats the two as alternatives); `randomize=False` stays deterministic and ignores the seed.

If the optional package behind `'louvain'` or `'leiden'` is not installed, `detect_communities` logs a warning and falls back instead of failing: Leiden falls back to Louvain when `python-louvain` is available, and everything else falls back to greedy modularity. The Louvain fallback keeps only the options Louvain understands (`resolution`, `randomize`, `random_state`, `weight` and `partition`), so a Leiden-only argument such as `n_iterations` is dropped and the call still returns Louvain communities. The greedy-modularity fallback likewise keeps only `weight`, `resolution`, `cutoff` and `best_n`. A method that raises for any other reason is also replaced by greedy modularity (called without options), so keyword arguments the chosen method does not accept lead to that fallback rather than an error.

Each algorithm is also available on its own and returns a list of communities (a list of node lists, largest first): `louvain_communities`, `leiden_communities`, `greedy_modularity_communities`, `girvan_newman_communities`, `label_propagation_communities`, `asyn_lpa_communities` and `fluid_communities`.

### hierarchical_communities()

Detect a hierarchy of partitions by repeatedly contracting each community into a single super-node and detecting communities again. The result maps each level (starting at 0) to a `{node: community_id}` dictionary. Level 0 holds the original nodes; at every later level the keys are the super-nodes `community_0`, `community_1`, ... (the previous level's communities, numbered from largest to smallest):

```python
from metainformant.networks.analysis.community import hierarchical_communities

levels = hierarchical_communities(network, method='greedy', max_levels=3)

for level, assignment in levels.items():
    print(f"Level {level}: {len(set(assignment.values()))} communities")
```

### modularity()

Calculate modularity score for community assignments (a node-to-community dictionary or a list of node lists):

```python
from metainformant.networks.analysis.community import modularity

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
from metainformant.networks.analysis.community import community_metrics

# Get detailed community statistics
metrics = community_metrics(network, communities)

print(f"Number of communities: {metrics['n_communities']}")
print(f"Community size range: {metrics['community_sizes']['min']}-{metrics['community_sizes']['max']}")
print(f"Largest community: {metrics['community_sizes']['max']} nodes")
print(f"Modularity: {metrics['modularity']:.3f}")
```

**Metrics include:**
- Number of communities (`n_communities`)
- Community size distribution (`community_sizes`: mean, std, min, max and the list of sizes)
- Modularity score (`modularity`)
- Fraction of nodes assigned to a community (`coverage`)
- Internal vs external connectivity (`internal_edge_ratio`, `conductance`)
- Community density (`internal_density`)
- A composite `quality_score`

## Stability and Resolution

Check how reproducible a partition is across runs, and scan resolutions for the best modularity:

```python
from metainformant.networks.analysis.community import community_stability, optimize_resolution

# Repeat detection with different seeds and compare the partitions
stability = community_stability(network, method='louvain', n_runs=5, seed=42)
print(f"Stability (mean NMI between runs): {stability['stability_score']:.3f}")
print(f"Mean modularity: {stability['avg_modularity']:.3f}")

# Scan resolutions between 0.5 and 2.0 and keep the best modularity
scan = optimize_resolution(network, resolution_range=(0.5, 2.0), n_points=5, method='louvain')
print(f"Best resolution: {scan['optimal_resolution']:.2f} (modularity {scan['optimal_modularity']:.3f})")
```

## Functional Analysis

### Community Enrichment

Analyze functional enrichment within communities by testing each community's nodes against a `PathwayNetwork` (see [Pathway Analysis](./pathway.md)):

```python
from metainformant.networks.analysis.pathway import PathwayNetwork, pathway_enrichment

# Get community assignments
communities = detect_communities(network, method='louvain')

# Pathway annotations: pathway ID -> genes
pathway_network = PathwayNetwork(
    name="toy_pathways",
    pathways={
        "process_one": [f"gene{i}" for i in range(0, 5)],
        "process_two": [f"gene{i}" for i in range(5, 10)],
        "process_three": [f"gene{i}" for i in range(10, 15)],
    },
)

# Functional enrichment analysis
for community_id in sorted(set(communities.values())):
    community_nodes = [node for node, comm in communities.items() if comm == community_id]

    # Analyze enrichment for this community (a dict keyed by pathway ID)
    enrichment = pathway_enrichment(
        community_nodes,
        pathway_network,
        background_genes=list(network.nodes()),
        min_overlap=2,
    )
    ranked = sorted(enrichment.items(), key=lambda item: item[1]["p_value"])
    print(f"Community {community_id}: {len(community_nodes)} nodes")
    for pathway_id, result in ranked[:5]:  # Top 5 pathways
        print(f"  {pathway_id}: overlap={result['overlap_size']}, p={result['p_value']:.2e}")
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

## Integration Examples

### With Protein Interaction Networks

```python-snippet
import pandas as pd
from metainformant.networks.interaction.ppi import load_string_interactions
from metainformant.networks.analysis.community import detect_communities
from metainformant.ontology.core import go

# Load PPI network from a STRING-style table (columns protein1, protein2, combined_score)
string_table = pd.read_csv("string_links.tsv", sep="\t")
protein_network = load_string_interactions(string_table, confidence_threshold=700)

# Detect functional modules (protein -> module id)
membership = detect_communities(protein_network, method='leiden')
modules = {}
for protein, module_id in membership.items():
    modules.setdefault(module_id, []).append(protein)

# GO term annotations: GO term -> set of annotated proteins (placeholder values; build yours from an annotation table)
go_annotations = {"GO:0016301": {"P1", "P2", "P3"}, "GO:0005215": {"P4", "P5", "P6"}}

# Functional enrichment of modules
for module_id, proteins in modules.items():
    enriched_terms = go.enrich_genes(proteins, protein_network.get_proteins(), go_annotations)
    print(f"Module {module_id}: {enriched_terms['go_term'].head(3).tolist()}")
```

### With Gene Expression Data

```python-snippet
import pandas as pd
from metainformant.networks.interaction.regulatory_analysis import infer_grn
from metainformant.networks.analysis.community import detect_communities

# expression_data: samples as rows, genes as columns
expression_data = pd.read_csv("expression_matrix.csv", index_col=0)

# Build regulatory network
reg_network = infer_grn(expression_data, list(expression_data.columns), method="correlation", threshold=0.7)

# Find regulatory modules (gene -> module id)
membership = detect_communities(reg_network)
reg_modules = {}
for gene, module_id in membership.items():
    reg_modules.setdefault(module_id, []).append(gene)

# Analyze module expression patterns
for module_id, genes in reg_modules.items():
    module_expression = expression_data[genes]
    # Analyze coordinated expression patterns
```

## Testing

Community detection functionality is tested in `tests/networks/test_networks_community.py`:

```bash
# Run community detection tests
uv run pytest tests/networks/test_networks_community.py -v

# Test specific algorithms
uv run pytest tests/networks/test_networks_community.py::TestCommunityDetection::test_detect_communities_different_algorithms -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
- [Regulatory Networks](./regulatory.md): Gene regulatory network analysis
