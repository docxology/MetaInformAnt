# SPEC: Network Scripts

Scripts for biological network construction and community analysis.

## Workflows

- `build_ppi_network.py`: Constructs protein-protein interaction networks from database exports.
- `detect_communities.py`: Runs Louvain or Leiden algorithms on network graphs.

## Standards

- **Graph Handling**: Uses `networkx` objects internally as facilitated by `metainformant.networks.graph`.
- **Visualization**: Integrates with `metainformant.visualization.genomics.networks` for high-quality output.
