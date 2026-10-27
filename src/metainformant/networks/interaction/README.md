# Networks Interaction

Protein-protein interaction (PPI) networks and gene regulatory network construction and analysis.

## Contents

| File | Purpose |
|------|---------|
| `ppi.py` | Compatibility facade re-exporting the PPI implementation below |
| `_ppi_impl.py` | PPI network loading/saving, hub detection, comparison, enrichment, interaction prediction, `ProteinNetwork` |
| `regulatory.py` | Backward-compatible facade re-exporting `regulatory_core` and `regulatory_analysis` |
| `regulatory_core.py` | `GeneRegulatoryNetwork` class and network-level regulatory operations |
| `regulatory_analysis.py` | GRN inference from expression, motif/cascade detection, validation |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_ppi_network()` | Load a PPI network from a TSV, CSV, JSON (node-link), BioPlex or IntAct file |
| `construct_ppi_network_from_interactions()` | Build network from interaction tuples |
| `ppi_network_analysis()` | Degree distribution, connectivity, clustering, topology, and centrality (the 10 most central proteins, most central first) |
| `find_ppi_hubs()` | Identify highly connected hub proteins |
| `ppi_network_clustering()` | Cluster PPI network into functional modules |
| `ProteinNetwork` | PPI network class: partners, confidence/evidence filtering, complexes, STRING I/O; `ProteinNetwork.from_file()` accepts the same five formats as `load_ppi_network()` |
| `predict_interactions()` | Predict novel protein-protein interactions |
| `save_ppi_network()` | Export network to TSV, CSV, or JSON (node-link) |
| `infer_grn()` | Infer a `GeneRegulatoryNetwork` from a samples x genes matrix (`'correlation'` or `'mutual_info'`; takes `random_state`) |
| `infer_regulatory_network_from_expression()` | Infer a regulatory `networkx.DiGraph` from a genes x samples DataFrame (`'correlation'` or `'mutual_info'`; `'grnboost2'` is not implemented) |
| `regulatory_motifs()` | Feed-forward loops (TF1 -> TF2 -> target plus the direct TF1 -> target edge), auto-regulation, mutual regulation |
| `analyze_regulatory_dynamics()` | Noisy simulation of regulatory dynamics; takes `random_state` |

## Usage

```python
from metainformant.networks.interaction.ppi import find_ppi_hubs, load_ppi_network, ppi_network_analysis

network = load_ppi_network("data/ppi_interactions.tsv")  # format="tsv" (default), "csv", "json", "bioplex" or "intact"
stats = ppi_network_analysis(network)
hubs = find_ppi_hubs(network, percentile=95.0)
```

Gene regulatory network inference with a seed (only the `'mutual_info'` method is random; the global NumPy random state is neither seeded nor consumed):

```python
from metainformant.networks.interaction.regulatory_analysis import infer_grn, regulatory_motifs

grn = infer_grn(expression_matrix, gene_names, method="mutual_info", threshold=0.1, random_state=0)
motifs = regulatory_motifs(grn)
```

`infer_grn` supports only `'correlation'` and `'mutual_info'`; for Granger causality on time series use `metainformant.information.network_info.information_flow.information_flow_network(series, method="granger_causality")`.
