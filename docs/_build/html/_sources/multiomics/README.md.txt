# Multi-Omics Documentation

This directory contains documentation for METAINFORMANT's multi-omics integration capabilities.

## Overview

The multi-omics module provides tools for integrating and analyzing data from multiple omic layers simultaneously, enabling comprehensive systems-level biological analysis.

## Documentation Files

### Core Multi-Omics Integration
- **`index.md`**: Multi-omics integration overview and module index

## Related Source Code

- See `src/metainformant/multiomics/` for implementation details
- See `tests/test_multiomics_*.py` for comprehensive test coverage

## Usage Examples

Multi-omics integration combines data from different biological domains:

```python
from metainformant.multiomics import integration

# Integrate RNA-seq and proteomics data
integrated_data = integration.combine_omics(rna_data, protein_data)
```

## Integration

Multi-omics analysis integrates with all domain modules:
- **DNA/RNA/Protein**: Molecular data integration
- **Epigenome**: Chromatin state analysis
- **Single-cell**: Multi-omic single-cell analysis
- **Visualization**: Multi-dimensional data visualization

This documentation provides complete coverage of METAINFORMANT's multi-omics integration capabilities.
