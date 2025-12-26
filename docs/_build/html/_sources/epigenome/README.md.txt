# Epigenome Documentation

This directory contains comprehensive documentation for METAINFORMANT's epigenetic modification analysis capabilities.

## Overview

The epigenome domain provides tools for DNA methylation, chromatin structure, and histone modification analysis.

## Documentation Files

### Core Epigenome Tools
- **`index.md`**: Epigenome domain overview and module index

## Related Source Code

- See `src/metainformant/epigenome/` for implementation details
- See `tests/test_epigenome_*.py` for comprehensive test coverage
- See `src/metainformant/epigenome/README.md` for module-specific documentation

## Usage Examples

The epigenome domain supports epigenetic analysis:

```python
from metainformant.epigenome import methylation

# Methylation pattern analysis
methylation_data = methylation.load_bismark_output("methylation_calls.txt")
differential = methylation.find_differential_methylation(methylation_data, groups)
```

## Integration

Epigenome analysis integrates with:
- **DNA analysis** for genomic context
- **RNA analysis** for expression correlation
- **Statistical methods** for methylation analysis
- **Visualization** for epigenetic pattern plotting

## Testing

Comprehensive tests ensure epigenetic reliability:
- Methylation calling algorithm validation
- Chromatin accessibility analysis
- Genomic annotation integration
- Statistical method verification

## Contributing

When adding new epigenome functionality:
1. Update epigenetic analysis documentation
2. Add comprehensive methylation tests
3. Ensure compatibility with epigenomic formats
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's epigenome capabilities.
