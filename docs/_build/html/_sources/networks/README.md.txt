# Networks Documentation

This directory contains comprehensive documentation for METAINFORMANT's network analysis capabilities.

## Overview

The networks module provides tools for graph theory, network analysis, and biological network modeling across different domains including protein-protein interactions, regulatory networks, and community detection.

## Documentation Files

### Network Analysis
- **`index.md`**: Networks domain overview and module index
- **`graph.md`**: General graph theory and network analysis
- **`community.md`**: Community detection algorithms (Louvain, Leiden)
- **`ppi.md`**: Protein-protein interaction networks
- **`pathway.md`**: Biological pathway analysis
- **`regulatory.md`**: Gene regulatory network modeling

## Related Source Code

- See `src/metainformant/networks/` for implementation details
- See `tests/test_networks_*.py` for comprehensive test coverage

## Usage Examples

Network analysis across biological domains:

```python
from metainformant.networks import graph, community

# Analyze protein-protein interaction network
ppi_network = graph.load_network("ppi_data.txt")
communities = community.detect_communities(ppi_network, method="louvain")
```

## Integration

Networks analysis integrates with:
- **Protein**: Protein interaction networks
- **Ontology**: Functional annotation networks
- **Visualization**: Network visualization and layout
- **Single-cell**: Cell-cell communication networks

## Testing

Comprehensive tests ensure algorithm correctness:
- Graph algorithms validation
- Community detection verification
- Network metrics calculation
- Integration testing with biological data

This documentation provides complete coverage of METAINFORMANT's network analysis capabilities.
