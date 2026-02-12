# Networks Analysis

Community detection, graph topology analysis, and pathway enrichment for biological networks.

## Contents

| File | Purpose |
|------|---------|
| `community.py` | Community detection: Louvain, Leiden, Girvan-Newman, label propagation |
| `graph.py` | Graph topology metrics and network structure analysis |
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
| `pathway_topology_analysis()` | Topology-aware enrichment scoring |
| `PathwayNetwork` | Class for pathway database operations and enrichment |

## Usage

```python
from metainformant.networks.analysis.community import detect_communities, evaluate_communities
from metainformant.networks.analysis.pathway import pathway_enrichment_analysis

communities = detect_communities(graph, method="louvain")
metrics = evaluate_communities(graph, communities)
enrichment = pathway_enrichment_analysis(gene_list, pathway_db)
```
