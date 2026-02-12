# Ontology Visualization

Comprehensive visualization for ontologies including GO DAG rendering, semantic similarity matrices, enrichment plots, and functional annotation networks. Supports matplotlib, seaborn, and optional Plotly for interactive output.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `visualization` module |
| `visualization.py` | All ontology visualization functions |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_go_dag()` | Render Gene Ontology directed acyclic graph |
| `plot_semantic_similarity_matrix()` | Heatmap of pairwise semantic similarity |
| `plot_go_enrichment_barplot()` | Bar plot of GO enrichment results |
| `plot_go_enrichment_dotplot()` | Dot plot of GO enrichment with size and color encoding |
| `plot_ontology_network()` | Network visualization of ontology terms |
| `plot_information_content_profile()` | Profile of information content across terms |
| `plot_go_term_hierarchy()` | Hierarchical tree visualization of GO terms |
| `plot_functional_annotation_heatmap()` | Heatmap of functional annotations across genes |
| `plot_semantic_similarity_clustermap()` | Clustered heatmap of semantic similarity |
| `create_interactive_go_network()` | Interactive Plotly network of GO terms |

## Usage

```python
from metainformant.ontology.visualization import visualization

visualization.plot_go_enrichment_barplot(enrichment_results, output_path="output/enrichment.png")
```
