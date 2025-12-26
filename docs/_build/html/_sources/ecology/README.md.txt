# Ecology Documentation

This directory contains comprehensive documentation for METAINFORMANT's ecological metadata and community analysis capabilities.

## Overview

The ecology domain provides tools for ecological metadata management, community composition analysis, and environmental data integration.

## Documentation Files

### Core Ecology Tools
- **`index.md`**: Ecology domain overview and module index

## Related Source Code

- See `src/metainformant/ecology/` for implementation details
- See `tests/test_ecology_*.py` for comprehensive test coverage
- See `src/metainformant/ecology/README.md` for module-specific documentation

## Usage Examples

The ecology domain supports ecological analysis:

```python
from metainformant.ecology import community

# Community diversity analysis
abundance_data = load_species_abundance()
diversity = community.shannon_diversity(abundance_data)
composition = community.analyze_composition(abundance_data)
```

## Integration

Ecology analysis integrates with:
- **DNA analysis** for genetic diversity context
- **Statistical methods** for community analysis
- **Visualization** for ecological data plotting
- **Environmental data** for context integration

## Testing

Comprehensive tests ensure ecological reliability:
- Diversity metric calculation validation
- Community analysis algorithm testing
- Environmental data integration
- Statistical method verification

## Contributing

When adding new ecology functionality:
1. Update ecological analysis documentation
2. Add comprehensive diversity metric tests
3. Ensure compatibility with ecological datasets
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's ecology capabilities.
