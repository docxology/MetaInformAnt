# Networks Module

The `networks` module provides tools for biological network analysis, including protein interaction networks, regulatory networks, pathway analysis, and community detection.

## Overview

This module handles various types of biological networks and provides algorithms for network construction, topology analysis, community detection, and functional enrichment. Supports both directed and undirected networks with weighted edges.

## Key Components

### Core Graph Operations (`graph.py`)
Fundamental network data structures and basic analysis.

**Usage:**
```python
from metainformant.networks import (
    create_network,
    network_metrics,
    centrality_measures,
    shortest_paths
)

# Create network
network = create_network(["A", "B", "C", "D"], directed=False)
network.add_edge("A", "B", weight=0.8)
network.add_edge("B", "C", weight=0.6)

# Calculate metrics
metrics = network_metrics(network)
print(f"Nodes: {metrics['num_nodes']}, Edges: {metrics['num_edges']}")
print(f"Density: {metrics['density']:.3f}")

# Centrality analysis
centralities = centrality_measures(network)
print(f"Hub nodes: {sorted(centralities['degree'].items(), key=lambda x: x[1], reverse=True)[:3]}")
```

### Protein-Protein Interaction Networks (`ppi.py`)
Analysis of protein-protein interaction networks.

**Key Features:**
- PPI network construction from databases (STRING format)
- Network topology analysis and hub identification
- Interaction prediction from feature data
- Functional enrichment analysis

**Usage:**
```python
from metainformant.networks import (
    ProteinNetwork,
    load_string_interactions,
    predict_interactions
)
import numpy as np

# Load from STRING database
ppi_network = load_string_interactions(
    "string_interactions.tsv",
    score_threshold=700,  # High confidence
    limit_organisms=["9606"]  # Human only
)

# Create network with confidence filtering
network = ppi_network.create_network(min_confidence=0.7)

# Network statistics
stats = ppi_network.network_statistics()
print(f"Hub proteins: {stats['hub_proteins']}")

# Predict interactions from features
features = np.random.randn(100, 50)  # 100 proteins, 50 features
protein_ids = [f"P{i:05d}" for i in range(100)]
predicted_ppi = predict_interactions(
    features, protein_ids,
    method="correlation",
    threshold=0.8
)
```

### Regulatory Networks (`regulatory.py`)
Gene regulatory network inference and analysis.

**Key Features:**
- GRN inference from expression data (correlation, mutual information, Granger causality)
- Regulatory motif detection (feed-forward loops, feedback loops, bifans)
- Master regulator identification
- Pathway regulation analysis

**Usage:**
```python
from metainformant.networks import (
    GeneRegulatoryNetwork,
    infer_grn,
    regulatory_motifs
)
import numpy as np

# Infer regulatory network from expression data
expression = np.random.randn(100, 200)  # 100 samples, 200 genes
gene_names = [f"GENE_{i}" for i in range(200)]
tf_list = [f"GENE_{i}" for i in range(50)]  # First 50 are TFs

grn = infer_grn(
    expression, gene_names,
    method="correlation",
    tf_list=tf_list,
    threshold=0.75
)

# Regulatory statistics
stats = grn.regulatory_statistics()
print(f"Master regulators: {stats['master_regulators']}")

# Detect regulatory motifs
motifs = regulatory_motifs(grn, motif_types=["feed_forward_loop", "feedback_loop"])
print(f"Feed-forward loops: {len(motifs['feed_forward_loop'])}")
```

### Pathway Analysis (`pathway.py`)
Biological pathway enrichment and analysis.

**Key Features:**
- Pathway database loading (GMT, CSV, TSV formats)
- Pathway overlap analysis
- Gene set enrichment analysis
- Network-based pathway enrichment

**Usage:**
```python
from metainformant.networks import (
    PathwayNetwork,
    pathway_enrichment,
    load_pathway_database
)

# Create pathway network
pn = PathwayNetwork(name="KEGG_pathways")
pn.add_pathway(
    "path:00010",
    ["GENE1", "GENE2", "GENE3", "GENE4"],
    metadata={"name": "Glycolysis"}
)
pn.add_pathway("path:00020", ["GENE2", "GENE5", "GENE6"])

# Pathway enrichment
query_genes = ["GENE1", "GENE2", "GENE3"]
enrichment = pathway_enrichment(query_genes, pn)

for result in enrichment[:5]:  # Top 5 enriched pathways
    print(f"{result['pathway_id']}: fold_enrichment={result['fold_enrichment']:.2f}")

# Load from file
pathway_db = load_pathway_database("kegg_pathways.gmt", format="gmt")
```

### Community Detection (`community.py`)
Network community/module detection algorithms.

**Key Features:**
- Multiple algorithms (Louvain, Leiden, Greedy)
- Modularity calculation
- Community metrics and analysis

**Usage:**
```python
from metainformant.networks import (
    detect_communities,
    modularity,
    community_metrics
)

# Detect communities using Louvain algorithm
communities = detect_communities(
    network,
    method="louvain",
    resolution=1.0,
    seed=42
)

# Calculate modularity
mod = modularity(network, communities)
print(f"Modularity: {mod:.3f}")

# Community metrics
metrics = community_metrics(network, communities)
print(f"Number of communities: {metrics['num_communities']}")
print(f"Average community size: {metrics['avg_community_size']:.1f}")
```

## Integration with Other Modules

### With Multi-Omics Module
```python
from metainformant.multiomics import MultiOmicsData
from metainformant.networks import predict_interactions

# Use expression data to predict protein interactions
omics_data = MultiOmicsData(proteomics=protein_expression)
features = omics_data.get_layer("proteomics").values.T
protein_ids = omics_data.get_layer("proteomics").columns.tolist()

ppi = predict_interactions(features, protein_ids, method="coexpression")
```

### With Ontology Module
```python
from metainformant.networks import pathway_enrichment
from metainformant.ontology import go

# Enrich network modules with GO terms
network_modules = detect_communities(protein_network)
for module_id, module_genes in network_modules.items():
    pathways = pathway_enrichment(list(module_genes), pathway_network)
    # Further GO enrichment...
```
enriched = pathway_enrichment(gene_list, pathway_db, background_size=20000)
```

### Community Detection (`community.py`)
Network community detection and analysis.

**Key Features:**
- Louvain and Leiden clustering
- Hierarchical clustering
- Community functional analysis
- Community comparison across conditions

**Usage:**
```python
from metainformant.networks import detect_communities, modularity, community_metrics

# Detect communities
communities = detect_communities(network, method="louvain", resolution=1.0)

# Calculate modularity
mod = modularity(network, communities)

# Community metrics
metrics = community_metrics(network, communities)
```

## Integration with Other Modules

### With Ontology Module
```python
from metainformant.networks import detect_communities
from metainformant.ontology import go

# Functional analysis of network modules
communities = detect_communities(protein_network.network, method="louvain")
# Use communities for enrichment analysis
```

### With Visualization Module
```python
from metainformant.networks import ProteinNetwork, load_string_interactions
from metainformant.visualization import plots

# Visualize interaction networks
ppi_network = load_string_interactions("string_file.tsv", score_threshold=400)
network = ppi_network.create_network()
# Use visualization.plots for plotting network graphs
```

## Performance Features

- Efficient graph algorithms for large networks
- Parallel community detection
- Memory-optimized network storage
- Streaming processing for very large networks

## Testing

Comprehensive tests cover:
- Network construction accuracy
- Algorithm correctness validation
- Performance with large networks
- Integration with external databases

## Dependencies

- NetworkX for graph algorithms
- Optional: igraph for high-performance graph operations

This module provides comprehensive tools for biological network analysis and interpretation.
