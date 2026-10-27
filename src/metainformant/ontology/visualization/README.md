# Ontology Visualization

Comprehensive visualization for ontologies including GO DAG rendering, semantic similarity matrices, enrichment plots, and functional annotation networks. Supports matplotlib, seaborn, and optional Plotly for interactive output.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `enrichment_dotplot`, `enrichment_score_plot` and `pathway_network_plot` from `plots.py`; the `visualization` module is importable as a submodule |
| `visualization.py` | Ontology visualization functions (DAG, similarity, enrichment, network plots) |
| `plots.py` | Enrichment dot plot, enrichment score plot and pathway similarity network plot |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_go_dag()` | Render a Gene Ontology directed acyclic graph (a `networkx.DiGraph`); `layout_seed` fixes the layout |
| `plot_semantic_similarity_matrix()` | Heatmap of pairwise semantic similarity |
| `plot_go_enrichment_barplot()` | Bar plot of GO enrichment results |
| `plot_go_enrichment_dotplot()` | Dot plot of GO enrichment with size and color encoding |
| `plot_ontology_network()` | Network visualization of ontology terms; `layout_seed` fixes the layout |
| `plot_information_content_profile()` | Profile of information content across terms |
| `plot_go_term_hierarchy()` | Hierarchical tree visualization of GO terms |
| `plot_functional_annotation_heatmap()` | Heatmap of functional annotations across genes |
| `plot_semantic_similarity_clustermap()` | Clustered heatmap of semantic similarity |
| `create_interactive_go_network()` | Interactive Plotly network of GO terms; `layout_seed` fixes the layout |
| `pathway_network_plot()` | Pathway similarity network from `pathway_network()` output (`plots.py`); `layout_seed` fixes the layout (default 42) |

## Usage

```python
from metainformant.ontology.visualization import visualization

visualization.plot_go_enrichment_barplot(enrichment_results, output_path="output/enrichment.png")
```

## Layout seeds

The graph layouts start from random positions. `plot_go_dag`, `plot_ontology_network` and
`create_interactive_go_network` (in `visualization.py`) take `layout_seed`, default `None` (unseeded: the layout draws
from NumPy's global generator, as networkx does). `pathway_network_plot` (in `plots.py`) takes `layout_seed` with the
default `42`. The same seed gives the same drawing, and a seed never seeds or consumes the global random state. The
defaults differ by design: the `visualization.py` plots stay unseeded unless a seed is given, `pathway_network_plot` is
seeded unless another seed is given.
