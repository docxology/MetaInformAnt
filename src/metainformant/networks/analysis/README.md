# Networks Analysis

Community detection, graph topology analysis, and pathway enrichment for biological networks.

## Contents

| File | Purpose |
|------|---------|
| `graph_core.py` | `BiologicalNetwork` class, network construction from edge lists/DataFrames, file I/O, adjacency conversion, validation |
| `graph_algorithms.py` | Network metrics, similarity, union/intersection, filtering, centrality, shortest paths |
| `graph.py` | Re-export facade combining `graph_core` and `graph_algorithms` for backward compatibility |
| `community.py` | Community detection: Louvain, Leiden, Girvan-Newman, label propagation |
| `pathway.py` | Pathway enrichment, topology analysis, and pathway network construction |

## Key Functions

| Function | Description |
|----------|-------------|
| `louvain_communities()` | Community detection via Louvain algorithm |
| `leiden_communities()` | Community detection via Leiden algorithm |
| `detect_communities()` | Unified interface for multiple detection methods |
| `evaluate_communities()` | Modularity, conductance, and coverage metrics |
| `compare_community_methods()` | Benchmark multiple algorithms on the same graph |
| `pathway_enrichment_analysis()` | Over-representation analysis for gene sets |
| `pathway_topology_analysis()` | Topology metrics of a pathway graph (degree, clustering, density, components) |
| `PathwayNetwork` | Class for pathway database operations and enrichment |

## Usage

```python
from metainformant.networks.analysis.community import detect_communities, evaluate_communities
from metainformant.networks.analysis.pathway import pathway_enrichment_analysis

communities = detect_communities(graph, method="louvain")
metrics = evaluate_communities(graph, communities)
enrichment = pathway_enrichment_analysis(
    genes=gene_list,
    background_genes=all_genes,
    pathways=pathway_db,
)
```
