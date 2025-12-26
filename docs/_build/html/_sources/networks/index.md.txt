# Networks: Biological Network Analysis

The networks module provides comprehensive tools for analyzing biological networks including protein interaction networks, regulatory networks, pathway analysis, and community detection.

## Overview

This module handles various types of biological networks and provides algorithms for network construction, analysis, and visualization. It supports multiple network types and provides both basic graph operations and advanced biological network analysis.

## Core Components

### [Graph Analysis](./graph.md)
Fundamental graph operations, metrics, and utilities:
- Network creation from various data sources
- Centrality measures (degree, betweenness, closeness)
- Shortest path analysis
- Network filtering and manipulation
- Performance optimization for large networks

### [Community Detection](./community.md)
Algorithms for identifying functional modules and communities:
- Louvain and Leiden clustering algorithms
- Modularity optimization
- Community functional analysis
- Hierarchical community detection
- Community stability assessment

### [Protein-Protein Interaction Networks](./ppi.md)
Analysis of protein interaction networks:
- STRING database integration
- Interaction prediction methods
- Network topology analysis
- Functional enrichment of modules
- Interaction confidence filtering

### [Gene Regulatory Networks](./regulatory.md)
Regulatory network inference and analysis:
- Network inference from expression data
- Transcription factor metadata
- Regulatory motif detection
- Regulation strength analysis
- Feed-forward and feedback loop identification

### [Pathway Analysis](./pathway.md)
Biological pathway enrichment and analysis:
- KEGG, Reactome, and GO pathway integration
- Enrichment analysis algorithms
- Pathway activity inference
- Multi-omics pathway analysis
- Pathway network construction

## Architecture

```mermaid
flowchart TD
    A[Input Data] --> B[Network Construction]
    B --> C[Graph Analysis]
    C --> D[Community Detection]
    C --> E[Pathway Enrichment]
    C --> F[Regulatory Analysis]
    D --> G[Functional Analysis]
    E --> G
    F --> G
    G --> H[Visualization]
    H --> I[Biological Insights]
```

## Key Features

### Comprehensive Network Types
- **Protein-Protein Interactions**: Physical and functional protein associations
- **Regulatory Networks**: Transcriptional regulation and signaling
- **Pathway Networks**: Metabolic and signaling pathways
- **Correlation Networks**: Gene co-expression and correlation patterns
- **Community Networks**: Functional modules and complexes

### Multiple Analysis Methods
- **Community Detection**: Louvain, Leiden, hierarchical, spectral clustering
- **Centrality Analysis**: Degree, betweenness, closeness, eigenvector centrality
- **Enrichment Analysis**: Hypergeometric, Fisher's exact, binomial tests
- **Network Inference**: Correlation, regression, mutual information, Granger causality

### Integration Capabilities
- **Multi-omics Integration**: Combine networks across different data types
- **Functional Annotation**: Gene Ontology and pathway enrichment
- **Expression Correlation**: Link networks to expression patterns
- **Database Integration**: STRING, KEGG, Reactome, GO databases

## Quick Start

### Basic Network Analysis

```python
from metainformant.networks import create_network, add_edges_from_interactions, network_metrics, centrality_measures

# Create network from interactions
nodes = ["gene1", "gene2", "gene3"]
network = create_network(nodes, directed=False)

# Add edges from interaction list
interactions = [
    ("gene1", "gene2", 0.8),
    ("gene2", "gene3", 0.6),
    ("gene3", "gene1", 0.9)
]

add_edges_from_interactions(network, interactions)

# Analyze network properties
metrics = network_metrics(network)
print(f"Network has {metrics['num_nodes']} nodes and {metrics['num_edges']} edges")

# Calculate centrality measures
centrality = centrality_measures(network)
hub_genes = [node for node, score in centrality['degree'].items()
             if score > centrality['degree'].mean() + centrality['degree'].std()]
```

### Community Detection

```python
from metainformant.networks import detect_communities, modularity

# Detect communities
communities = detect_communities(network, method='leiden', resolution=1.0)

# Evaluate community quality
mod_score = modularity(network, communities)
print(f"Modularity: {mod_score:.3f}")

# Analyze community properties
for community_id in set(communities.values()):
    community_nodes = [node for node, comm in communities.items() if comm == community_id]
    print(f"Community {community_id}: {len(community_nodes)} nodes")
```

### Pathway Enrichment

```python
from metainformant.networks import pathway_enrichment

# Gene list of interest
gene_list = ["HK1", "PFKL", "FBA", "PGK1", "CS", "ACO1"]

# Pathway enrichment analysis
enriched_pathways = pathway_enrichment(
    gene_list,
    pathway_network,
    method="hypergeometric",
    correction="fdr"
)

# Display top enriched pathways
for pathway in enriched_pathways[:5]:
    print(f"{pathway['name']}: p={pathway['p_value']:.2e}")
    print(f"  Overlap: {pathway['overlap_genes']}")
```

## Integration with Other Modules

### With Expression Data

```python
from metainformant.rna import workflow
from metainformant.networks import infer_grn

# Load expression data
expression_data = workflow.extract_expression_patterns(rna_data)

# Infer regulatory network
regulatory_network = infer_grn(
    expression_data,
    method="correlation",
    threshold=0.7
)

# Analyze network modules
modules = detect_communities(regulatory_network)
```

### With Protein Data

```python
from metainformant.networks import ppi, community
from metainformant.ontology import go

# Load PPI network
protein_network = ppi.load_string_network()

# Find functional modules
modules = community.detect_communities(protein_network, method='leiden')

# Functional enrichment
for module_id, proteins in modules.items():
    enriched_terms = go.enrich_proteins(proteins)
    print(f"Module {module_id}: {enriched_terms[:3]}")
```

## Performance Features

- **Scalable Algorithms**: Efficient implementations for large networks
- **Parallel Processing**: Multi-core support for computationally intensive operations
- **Memory Optimization**: Streaming processing for very large networks
- **Sparse Representations**: Memory-efficient storage for sparse networks

## Testing

Comprehensive tests cover all network functionality:

```bash
# Run all network tests
uv run pytest tests/test_networks_*.py -v

# Test specific components
uv run pytest tests/test_networks_graph.py::test_centrality_measures -v
uv run pytest tests/test_networks_community.py::test_detect_communities_louvain -v
```

## Related Documentation

- [Core Graph Operations](./graph.md): Basic network operations and metrics
- [Community Detection](./community.md): Network module identification
- [Protein Interaction Networks](./ppi.md): PPI network analysis
- [Regulatory Networks](./regulatory.md): Gene regulation analysis
- [Pathway Analysis](./pathway.md): Biological pathway enrichment
