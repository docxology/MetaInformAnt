# Networks Interaction

Protein-protein interaction (PPI) networks and gene regulatory network construction and analysis.

## Contents

| File | Purpose |
|------|---------|
| `ppi.py` | PPI network loading, hub detection, clustering, and interaction prediction |
| `regulatory.py` | Gene regulatory network inference and analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_ppi_network()` | Load PPI network from TSV or other formats |
| `construct_ppi_network_from_interactions()` | Build network from interaction tuples |
| `ppi_network_analysis()` | Degree distribution, centrality, and topology metrics |
| `find_ppi_hubs()` | Identify highly connected hub proteins |
| `ppi_network_clustering()` | Cluster PPI network into functional modules |
| `ProteinNetwork` | Class with interaction prediction (similarity, correlation, ML) |
| `predict_interactions()` | Predict novel protein-protein interactions |
| `save_ppi_network()` | Export network to TSV or other formats |

## Usage

```python
from metainformant.networks.interaction.ppi import load_ppi_network, ppi_network_analysis

network = load_ppi_network("data/ppi_interactions.tsv")
stats = ppi_network_analysis(network)
hubs = find_ppi_hubs(network, top_n=20)
```
