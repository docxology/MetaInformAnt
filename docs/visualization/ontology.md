# Ontology Visualization

This document provides comprehensive documentation for gene ontology and functional annotation visualization capabilities in METAINFORMANT.

## Overview

Ontology visualization includes specialized plots for gene ontology (GO) analysis, semantic similarity, functional enrichment, and biological knowledge representation. These tools create publication-quality figures for functional genomics research.

## Module Functions

### Gene Ontology Structure

#### GO DAG Visualization
```python
from metainformant.ontology.core import go
from metainformant.ontology.query.serialize import ontology_to_graph
from metainformant.ontology.visualization import visualization as ont_viz
import networkx as nx

# Load the GO ontology; the DAG and network plots take a networkx DiGraph
# (every edge points from a term to its parent)
ontology = go.load_go_obo("go.obo")
go_graph = ontology_to_graph(ontology)

# Plot GO directed acyclic graph (layout_seed makes the layout reproducible)
ax = ont_viz.plot_go_dag(go_graph, terms=['GO:0008150', 'GO:0003674', 'GO:0005575'], layout_seed=0)
```

#### GO Term Hierarchy
```python
# Visualize the hierarchy below a root term. The plot follows each node's successors, and
# ontology_to_graph points every edge from a term to its parent, so reverse the graph first
ax = ont_viz.plot_go_term_hierarchy(go_graph.reverse(), root_term='GO:0008150', max_depth=3)
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
    {'term': 'GO:0008150', 'pvalue': 1e-15, 'gene_count': 24},
    {'term': 'GO:0003674', 'pvalue': 1e-12, 'gene_count': 18},
    {'term': 'GO:0005575', 'pvalue': 1e-8, 'gene_count': 9}
]

ax = ont_viz.plot_go_enrichment_barplot(enrichment_results)
```

The plot functions take a list of dictionaries with a `term` and a `pvalue` entry (the dot plot also reads `gene_count`, default 1).

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
ax = ont_viz.plot_ontology_network(go_graph, figsize=(12, 8), layout_seed=0)
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
fig = ont_viz.create_interactive_go_network(go_graph, output_path="go_network.html", layout_seed=0)
```

#### Reproducible Layouts
The graph layouts start from random positions, so an unseeded layout draws the same graph differently on every call.
`plot_go_dag`, `plot_ontology_network` and `create_interactive_go_network` take `layout_seed`:

- The default is `None`, which leaves the layout unseeded: it then draws from NumPy's global generator, as networkx
  does.
- With an integer seed the same graph gives the same node positions (`plot_go_dag` also seeds its fallback layout),
  and the global random state is neither seeded nor consumed.
- The pathway similarity network plot, `pathway_network_plot` in `metainformant.ontology.visualization`, takes
  `layout_seed` as well, but its default is `42`, so that plot is reproducible without an argument. The defaults differ
  by design: `visualization.py` stays unseeded unless a seed is given, `plots.py` is seeded unless another seed is given.

```python
from metainformant.ontology.visualization import pathway_network_plot

# network: the dictionary returned by metainformant.ontology.pathway_enrichment.enrichment.pathway_network
fig = pathway_network_plot(network, "pathway_network.png", layout_seed=7)  # layout_seed defaults to 42
```

## Integration with Ontology Module

### With GO Enrichment Analysis
```python
from metainformant.ontology.core import go
from metainformant.ontology.visualization import visualization as ont_viz

# Perform GO enrichment
gene_list = ['gene1', 'gene2', 'gene3']  # Your genes of interest
background = [f'gene{i}' for i in range(1, 1001)]  # Background gene set
annotations = {  # GO term -> genes annotated with it
    'GO:0008150': {'gene1', 'gene2', 'gene10', 'gene11'},
    'GO:0003674': {'gene3', 'gene20', 'gene21'},
}

results = go.enrich_genes(gene_list, background, annotations)  # a DataFrame, one row per tested term

# The plot functions take a list of dictionaries with 'term', 'pvalue' and (for the dot plot) 'gene_count'
enrichment_results = [
    {'term': row.go_term, 'pvalue': row.p_value, 'gene_count': int(row.observed)}
    for row in results.itertuples()
]

# Visualize enrichment
ax = ont_viz.plot_go_enrichment_barplot(enrichment_results)
```

### With Semantic Similarity
```python
import numpy as np
from metainformant.ontology.core import go
from metainformant.ontology.visualization import visualization as ont_viz

# Load the ontology and compute information content (bits) from annotation counts
ontology = go.load_go_obo("go.obo")
term_counts = {'GO:0008150': 40, 'GO:0003674': 25, 'GO:0005575': 10}  # annotations per GO term
term_ic = go.calculate_term_ic(ontology, term_counts, total_annotations=sum(term_counts.values()))
hierarchy = go.build_hierarchy_dict(ontology)  # term -> set of ancestors

# Pairwise Lin similarity (between 0 and 1; Resnik scores are in bits and not bounded by 1)
terms = ['GO:0008150', 'GO:0003674', 'GO:0005575']
similarity_matrix = np.zeros((len(terms), len(terms)))

for i, term1 in enumerate(terms):
    for j, term2 in enumerate(terms):
        similarity_matrix[i, j] = go.semantic_similarity(term1, term2, term_ic, hierarchy, method="lin")

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
import matplotlib.pyplot as plt
import numpy as np
from metainformant.ontology.core import go
from metainformant.ontology.query.serialize import ontology_to_graph
from metainformant.ontology.visualization import visualization as ont_viz

# Load the ontology: the DAG plot needs a networkx DiGraph, the similarity functions the ontology itself
ontology = go.load_go_obo("go.obo")
go_graph = ontology_to_graph(ontology)
hierarchy = go.build_hierarchy_dict(ontology)

# Perform enrichment analysis (annotations: GO term -> genes annotated with it)
gene_list = ['gene1', 'gene2', 'gene3', 'gene350']
background = [f'gene{i}' for i in range(1000)]
annotations = {'GO:0008150': set(background[:400]), 'GO:0003674': set(background[300:700])}

results = go.enrich_genes(gene_list, background, annotations)
enrichment_results = [
    {'term': row.go_term, 'pvalue': row.p_value, 'gene_count': int(row.observed)}
    for row in results.itertuples()
]

# Information content (bits) from the annotation counts
term_counts = {term: len(genes) for term, genes in annotations.items()}
term_ic = go.calculate_term_ic(ontology, term_counts, total_annotations=sum(term_counts.values()))

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# GO DAG
ont_viz.plot_go_dag(go_graph, terms=['GO:0008150'], ax=axes[0,0], layout_seed=0)

# Enrichment results
ont_viz.plot_go_enrichment_barplot(enrichment_results, ax=axes[0,1])

# Semantic similarity of the enriched terms (Lin: values between 0 and 1)
terms = [r['term'] for r in enrichment_results]
sim_matrix = np.array([
    [go.semantic_similarity(t1, t2, term_ic, hierarchy, method="lin") for t2 in terms]
    for t1 in terms
])
ont_viz.plot_semantic_similarity_matrix(sim_matrix, terms, ax=axes[1,0])

# Information content profile
ont_viz.plot_information_content_profile({t: term_ic[t] for t in terms}, ax=axes[1,1])

plt.tight_layout()
plt.savefig("output/ontology_analysis.png", dpi=300, bbox_inches='tight')
```

## Troubleshooting

### Common Issues

1. **Missing GO Terms**: Ensure OBO file is current and terms exist
2. **Empty Enrichment Results**: Check gene annotation quality and background set
3. **Network Layout Issues**: Try different layout algorithms for complex graphs
4. **Memory Errors**: Subsample large ontologies for visualization

### Data Requirements

- **GO Terms**: Valid GO identifiers (GO:####### format)
- **P-values**: raw p-values in the `pvalue` entry (the bar and dot plots apply -log10 themselves)
- **Similarity Matrices**: Square matrices. Lin and Jiang-Conrath scores lie between 0 and 1; Resnik scores are
  information contents in bits and are not bounded by 1
- **IC Values**: Information content in bits (typically 0-15)

## Related Documentation

- **[Ontology Analysis](../ontology/)**: Core ontology analysis functions
- **[GO Enrichment](../ontology/go.md)**: Gene ontology enrichment workflows
- **[Semantic Similarity](../ontology/go.md)**: Similarity calculation methods
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive ontology visualization capabilities integrated with METAINFORMANT's functional annotation workflows.

