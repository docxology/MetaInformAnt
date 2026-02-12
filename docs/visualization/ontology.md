# Ontology Visualization

This document provides comprehensive documentation for gene ontology and functional annotation visualization capabilities in METAINFORMANT.

## Overview

Ontology visualization includes specialized plots for gene ontology (GO) analysis, semantic similarity, functional enrichment, and biological knowledge representation. These tools create publication-quality figures for functional genomics research.

## Module Functions

### Gene Ontology Structure

#### GO DAG Visualization
```python
from metainformant.ontology.core import go
from metainformant.ontology.visualization import visualization as ont_viz
import networkx as nx

# Load GO ontology
go_graph = go.load_go_obo("go.obo")

# Plot GO directed acyclic graph
ax = ont_viz.plot_go_dag(go_graph, terms=['GO:0008150', 'GO:0003674', 'GO:0005575'])
```

#### GO Term Hierarchy
```python
# Visualize hierarchy from root term
ax = ont_viz.plot_go_term_hierarchy(go_graph, root_term='GO:0008150', max_depth=3)
```

### Semantic Similarity

#### Similarity Matrices
```python
# Plot semantic similarity heatmap
term_ic = {'GO:0008150': 8.5, 'GO:0003674': 9.2, 'GO:0005575': 7.8}
similarity_matrix = np.random.rand(10, 10)  # Pre-computed similarity matrix
term_labels = [f'GO:{i}' for i in range(10)]

ax = ont_viz.plot_semantic_similarity_matrix(similarity_matrix, term_labels)
```

#### Similarity Clustermap
```python
# Create hierarchical clustering of similarities
g = ont_viz.plot_semantic_similarity_clustermap(similarity_matrix, term_labels)
```

### Functional Enrichment

#### Enrichment Bar Plots
```python
# Visualize GO enrichment results
enrichment_results = [
    {'term': 'GO:0008150', 'pvalue': 1e-15, 'fold_change': 3.2},
    {'term': 'GO:0003674', 'pvalue': 1e-12, 'fold_change': 2.8},
    {'term': 'GO:0005575', 'pvalue': 1e-8, 'fold_change': 2.1}
]

ax = ont_viz.plot_go_enrichment_barplot(enrichment_results)
```

#### Enrichment Dot Plots
```python
# Bubble plot with gene counts
ax = ont_viz.plot_go_enrichment_dotplot(enrichment_results)
```

### Information Content

#### IC Profile Plots
```python
# Plot information content across terms
term_ic = {'GO:0008150': 8.5, 'GO:0003674': 9.2, 'GO:0005575': 7.8}
ax = ont_viz.plot_information_content_profile(term_ic)
```

### Ontology Networks

#### Relationship Networks
```python
# Visualize ontology relationship network
ax = ont_viz.plot_ontology_network(go_graph, figsize=(12, 8))
```

#### Functional Annotation Heatmaps
```python
# Plot gene-term annotation matrix
annotation_matrix = np.random.randint(0, 2, (50, 20))  # Genes x GO terms
gene_labels = [f'Gene_{i}' for i in range(50)]
term_labels = [f'GO:{i}' for i in range(20)]

ax = ont_viz.plot_functional_annotation_heatmap(annotation_matrix, term_labels, gene_labels)
```

### Advanced Features

#### Interactive GO Networks
```python
# Create 3D interactive GO network
fig = ont_viz.create_interactive_go_network(go_graph, output_path="go_network.html")
```

## Integration with Ontology Module

### With GO Enrichment Analysis
```python
from metainformant.ontology.core import go
from metainformant.ontology.visualization import visualization as ont_viz

# Perform GO enrichment
gene_list = ['gene1', 'gene2', 'gene3']  # Your genes of interest
background = ['gene1', 'gene2', ..., 'gene1000']  # Background gene set
annotations = {'gene1': {'GO:0008150'}, 'gene2': {'GO:0003674'}}  # Gene annotations

results = go.enrich_genes(gene_list, background, annotations)

# Visualize enrichment
ax = ont_viz.plot_go_enrichment_barplot(results)
```

### With Semantic Similarity
```python
from metainformant.ontology.core import go
from metainformant.ontology.query import query
from metainformant.ontology.visualization import visualization as ont_viz

# Load ontology and compute IC
go_graph = go.load_go_obo("go.obo")
ic_map = query.calculate_ic_map(go_graph)

# Compute similarities
terms = ['GO:0008150', 'GO:0003674', 'GO:0005575']
similarity_matrix = np.zeros((len(terms), len(terms)))

for i, term1 in enumerate(terms):
    for j, term2 in enumerate(terms):
        similarity_matrix[i, j] = go.semantic_similarity(term1, term2, ic_map, go_graph)

# Visualize
ax = ont_viz.plot_semantic_similarity_matrix(similarity_matrix, terms)
```

### With Functional Annotations
```python
from metainformant.ontology.core import go
from metainformant.ontology.visualization import visualization as ont_viz

# Create annotation matrix
genes = ['gene1', 'gene2', 'gene3']
go_terms = ['GO:0008150', 'GO:0003674', 'GO:0005575']

# Binary annotation matrix (gene x term)
annotations = np.random.randint(0, 2, (len(genes), len(go_terms)))

ax = ont_viz.plot_functional_annotation_heatmap(annotations, go_terms, genes)
```

## Output Options

All visualization functions support:
```python
# Save static plots
ax = ont_viz.plot_go_dag(go_graph, output_path="go_dag.png")

# Save interactive networks
fig = ont_viz.create_interactive_go_network(go_graph, output_path="go_network.html")
```

## File Format Support

- **OBO Files**: Open Biological Ontology format for GO
- **GAF Files**: Gene Association File format for annotations
- **JSON/CSV**: Custom annotation and enrichment result formats
- **Network Formats**: GraphML, GML for ontology networks

## Performance Considerations

- **Large Ontologies**: GO has >40,000 terms; consider subgraph visualization
- **Dense Networks**: Use hierarchical layout for complex relationship networks
- **Similarity Matrices**: Work best with <100 terms for readability
- **Interactive Plots**: Require Plotly for web-based 3D networks

## Dependencies

- **Required**: matplotlib, numpy, networkx
- **Optional**: seaborn (enhanced styling), plotly (interactive plots)
- **Integration**: obonet (OBO parsing), goatools (enrichment analysis)

## Color Schemes

Recommended color schemes for ontology data:
```python
# GO aspect colors
go_colors = {
    'biological_process': '#FF6B6B',    # Red
    'molecular_function': '#4ECDC4',    # Teal
    'cellular_component': '#45B7D1'     # Blue
}

# Semantic similarity
similarity_cmap = 'YlOrRd'  # Yellow to red gradient

# Information content
ic_cmap = 'viridis'  # Blue to yellow colormap
```

## Examples

### Complete Ontology Analysis Workflow
```python
from metainformant.ontology.core import go
from metainformant.ontology.query import query
from metainformant.ontology.visualization import visualization as ont_viz
import numpy as np

# Load ontology
go_graph = go.load_go_obo("go.obo")
ic_map = query.calculate_ic_map(go_graph)

# Perform enrichment analysis
gene_list = ['gene1', 'gene2', 'gene3']
background = [f'gene{i}' for i in range(1000)]
annotations = {f'gene{i}': {'GO:0008150', 'GO:0003674'} for i in range(1000)}

enrichment_results = go.enrich_genes(gene_list, background, annotations)

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# GO DAG
ont_viz.plot_go_dag(go_graph, terms=['GO:0008150'], ax=axes[0,0])

# Enrichment results
ont_viz.plot_go_enrichment_barplot(enrichment_results, ax=axes[0,1])

# Semantic similarity (if multiple enriched terms)
if len(enrichment_results) > 1:
    terms = [r['term'] for r in enrichment_results]
    sim_matrix = np.random.rand(len(terms), len(terms))  # Placeholder
    ont_viz.plot_semantic_similarity_matrix(sim_matrix, terms, ax=axes[1,0])

# Information content profile
term_ic = {r['term']: ic_map.get(r['term'], 5.0) for r in enrichment_results}
ont_viz.plot_information_content_profile(term_ic, ax=axes[1,1])

plt.tight_layout()
plt.savefig("ontology_analysis.png", dpi=300, bbox_inches='tight')
```

## Troubleshooting

### Common Issues

1. **Missing GO Terms**: Ensure OBO file is current and terms exist
2. **Empty Enrichment Results**: Check gene annotation quality and background set
3. **Network Layout Issues**: Try different layout algorithms for complex graphs
4. **Memory Errors**: Subsample large ontologies for visualization

### Data Requirements

- **GO Terms**: Valid GO identifiers (GO:####### format)
- **P-values**: -log10 transformed for bar plot scaling
- **Similarity Matrices**: Square matrices with values 0-1
- **IC Values**: Information content in bits (typically 0-15)

## Related Documentation

- **[Ontology Analysis](../ontology/)**: Core ontology analysis functions
- **[GO Enrichment](../ontology/go.md)**: Gene ontology enrichment workflows
- **[Semantic Similarity](../ontology/go.md)**: Similarity calculation methods
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive ontology visualization capabilities integrated with METAINFORMANT's functional annotation workflows.

