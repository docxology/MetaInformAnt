# Information Integration

Information-theoretic analysis applied to biological data integration and network analysis, measuring entropy across DNA, RNA, single-cell, multi-omics, and ML pipeline outputs.

## Contents

| File | Purpose |
|------|---------|
| `integration.py` | Domain-specific entropy analysis for DNA, RNA, single-cell, multi-omics, ML |
| `networks.py` | Network entropy, information flow, community detection, centrality measures |

## Key Functions

| Function | Description |
|----------|-------------|
| `dna_integration()` | Compute entropy and information content of DNA sequence data |
| `rna_integration()` | Entropy analysis of RNA expression matrices |
| `singlecell_integration()` | Information metrics for single-cell expression profiles |
| `multiomics_integration()` | Cross-omics mutual information and integration scores |
| `ml_integration()` | Information-theoretic evaluation of ML feature sets |
| `network_entropy()` | Shannon or von Neumann entropy of a network graph |
| `information_flow()` | Random-walk or diffusion-based information flow between nodes |
| `information_community_detection()` | Community detection via infomap or map equation |
| `network_information_centrality()` | Entropy-based node centrality ranking |
| `network_motif_information()` | Information content of network motif patterns |
| `information_graph_distance()` | Entropy-based distance between two graphs |

## Usage

```python
from metainformant.information.integration.integration import dna_integration
from metainformant.information.integration.networks import network_entropy

results = dna_integration(sequences, method="entropy")
h = network_entropy(graph, method="shannon")
```
