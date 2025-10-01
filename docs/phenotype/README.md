# Phenotype Documentation

This directory contains comprehensive documentation for METAINFORMANT's phenotypic trait analysis and curation capabilities.

## Overview

The phenotype domain provides tools for morphological and behavioral phenotype data analysis, including automated curation from web sources like AntWiki.

## Documentation Files

### Core Phenotype Tools
- **`index.md`**: Phenotype domain overview and module index
- **`antwiki.md`**: AntWiki database integration and phenotype curation

## Related Source Code

- See `src/metainformant/phenotype/` for implementation details
- See `tests/test_phenotype_*.py` for comprehensive test coverage
- See `src/metainformant/phenotype/README.md` for module-specific documentation

## Usage Examples

The phenotype domain supports trait analysis:

```python
from metainformant.phenotype import antwiki

# Retrieve phenotype data from AntWiki
species_traits = antwiki.get_species_traits("Camponotus pennsylvanicus")
morphological_data = antwiki.extract_morphological_measurements(species_traits)
```

## Integration

Phenotype analysis integrates with:
- **DNA analysis** for genotype-phenotype associations
- **Ontology** for trait functional annotation
- **Statistical methods** for trait correlation analysis
- **Visualization** for trait distribution plotting

## Testing

Comprehensive tests ensure phenotype reliability:
- Web scraping functionality validation
- Data parsing and standardization
- Trait measurement accuracy
- Integration with other modules

## Contributing

When adding new phenotype functionality:
1. Update trait analysis documentation
2. Add comprehensive data validation tests
3. Ensure compatibility with phenotype databases
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's phenotype analysis capabilities.
