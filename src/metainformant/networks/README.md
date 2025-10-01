# Networks Module

The `networks` module provides tools for biological network analysis, including protein interaction networks, regulatory networks, and pathway analysis.

## Overview

This module handles various types of biological networks and provides algorithms for network construction, analysis, and visualization.

## Submodules

### Protein-Protein Interaction (`ppi.py`)
Analysis of protein-protein interaction networks.

**Key Features:**
- PPI network construction from databases
- Network topology analysis
- Module detection and functional enrichment
- Interaction prediction and validation

**Usage:**
```python
from metainformant.networks import ppi

# Load PPI network
network = ppi.load_string_network()
modules = ppi.find_modules(network)
enriched = ppi.functional_enrichment(modules)
```

### Regulatory Networks (`regulatory.py`)
Gene regulatory network inference and analysis.

**Key Features:**
- Transcription factor binding site analysis
- Regulatory motif detection
- Network inference from expression data
- Regulatory cascade identification

**Usage:**
```python
from metainformant.networks import regulatory

# Infer regulatory network
expression_data = load_expression_matrix()
regulators = regulatory.identify_regulators(expression_data)
network = regulatory.infer_network(expression_data, regulators)
```

### Pathway Analysis (`pathway.py`)
Biological pathway enrichment and analysis.

**Key Features:**
- Pathway database integration (KEGG, Reactome)
- Enrichment analysis algorithms
- Pathway visualization
- Cross-species pathway mapping

**Usage:**
```python
from metainformant.networks import pathway

# Pathway enrichment
gene_list = ["GENE1", "GENE2", "GENE3"]
enriched_pathways = pathway.enrich_pathways(gene_list, database="kegg")
visualization = pathway.visualize_pathway(enriched_pathways[0])
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
from metainformant.networks import community

# Detect communities
communities = community.louvain_clustering(network)
stability = community.community_stability(communities)
```

## Integration with Other Modules

### With Ontology Module
```python
from metainformant.networks import ppi
from metainformant.ontology import go

# Functional analysis of network modules
modules = ppi.find_modules(protein_network)
enriched_terms = go.enrich_modules(modules)
```

### With Visualization Module
```python
from metainformant.networks import ppi
from metainformant.visualization import network_plots

# Visualize interaction networks
network = ppi.load_network("ppi_data.txt")
network_plots.plot_interaction_network(network, layout="force_directed")
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
