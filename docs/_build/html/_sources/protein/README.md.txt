# Protein Documentation

This directory contains comprehensive documentation for METAINFORMANT's protein sequence and structure analysis capabilities.

## Overview

The protein domain provides tools for proteomic analysis, sequence manipulation, structure integration, and functional annotation.

## Documentation Files

### Core Protein Analysis
- **`index.md`**: Protein domain overview and module index
- **`proteomes.md`**: Proteome retrieval and analysis

## Related Source Code

- See `src/metainformant/protein/` for implementation details
- See `tests/test_protein_*.py` for comprehensive test coverage
- See `src/metainformant/protein/README.md` for module-specific documentation

## Usage Examples

The protein domain supports proteomic workflows:

```python
from metainformant.protein import proteomes, uniprot

# Retrieve and analyze proteomes
proteome = proteomes.get_proteome("UP000005640")
sequences = proteomes.extract_sequences(proteome)

# Functional annotation
protein_data = uniprot.get_protein_data("P12345")
domains = interpro.get_domains(protein_data["sequence"])
```

## Integration

Protein analysis integrates with:
- **DNA analysis** for translation and ORF finding
- **RNA workflows** for expression context
- **Structure analysis** for 3D modeling
- **Ontology** for functional annotation

## Testing

Comprehensive tests ensure proteomic reliability:
- Sequence format validation
- Database integration testing
- Structure parsing correctness
- Annotation accuracy verification

## Contributing

When adding new protein functionality:
1. Update sequence and structure documentation
2. Add comprehensive integration tests
3. Ensure compatibility with standard formats
4. Update database integration examples

This documentation provides complete coverage of METAINFORMANT's protein analysis capabilities.
