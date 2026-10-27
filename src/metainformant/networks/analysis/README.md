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
| `louvain_communities()` | Community detection via Louvain algorithm (`random_state` goes to python-louvain's own `random_state`) |
| `leiden_communities()` | Community detection via Leiden algorithm |
| `detect_communities()` | Unified interface for multiple detection methods |
| `evaluate_communities()` | Modularity, conductance, and coverage metrics |
| `compare_community_methods()` | Benchmark multiple algorithms on the same graph |
| `pathway_enrichment_analysis()` | Over-representation analysis for gene sets (needs scipy) |
| `pathway_enrichment()` | Enrichment of a gene list against a `PathwayNetwork` (`'fisher'` or `'hypergeometric'`, exact with or without scipy) |
| `pathway_similarity_analysis()` | Jaccard, overlap, or cosine similarity of every pathway pair |
| `create_pathway_network()` | Graph of pathways joined by similarity above a threshold |
| `pathway_topology_analysis()` | Topology metrics of a pathway graph (degree, clustering, density, components) |
| `PathwayNetwork` | Class for pathway database operations and enrichment |

### Community detection behavior

- **Seeds.** `louvain_communities(graph, random_state=...)` passes the seed to python-louvain's own `random_state`, so the same seed gives the same partition and the global `random` and `numpy.random` states are never seeded. With a seed, `randomize` is not passed on; `randomize=False` stays deterministic and ignores the seed. `leiden_communities` passes its seed to the Leiden optimiser.
- **Fallbacks in `detect_communities`.** When the package behind `'louvain'` or `'leiden'` is missing, the call falls back instead of failing. Leiden falls back to Louvain (when python-louvain is installed) with only the options Louvain accepts: `resolution`, `randomize`, `random_state`, `weight` and `partition`. Leiden-only options such as `n_iterations` are dropped, so the call returns Louvain communities. Every other fallback goes to greedy modularity, which keeps `weight`, `resolution`, `cutoff` and `best_n`.

### Pathway analysis behavior

- `pathway_enrichment(..., method=...)` accepts `'fisher'` and `'hypergeometric'` and raises `ValueError` for any other method. Both report the one-sided p-value `P(X >= overlap)`. Without scipy the exact upper tail is computed in integer arithmetic, and `correction="fdr"` (Benjamini-Hochberg) still runs.
- `pathway_similarity_analysis()` returns keys of the form `"{name1}_{name2}"`, which cannot be split back into names when a pathway ID contains an underscore. `create_pathway_network()` scores the pairs directly, so such IDs are safe.
- A pathway with no genes has similarity 0.0 with every pathway under all three similarity methods (the cosine of an empty set is defined as 0.0).

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
