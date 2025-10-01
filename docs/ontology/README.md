# Ontology Documentation

This directory contains comprehensive documentation for METAINFORMANT's functional annotation and semantic analysis capabilities.

## Overview

The ontology domain provides tools for functional annotation using biological ontologies like Gene Ontology (GO) and semantic similarity analysis.

## Documentation Files

### Core Ontology Tools
- **`index.md`**: Ontology domain overview and module index
- **`go.md`**: Gene Ontology annotation and enrichment analysis

## Related Source Code

- See `src/metainformant/ontology/` for implementation details
- See `tests/test_ontology_*.py` for comprehensive test coverage
- See `src/metainformant/ontology/README.md` for module-specific documentation

## Usage Examples

The ontology domain supports functional annotation:

```python
from metainformant.ontology import go

# Gene Ontology enrichment analysis
gene_list = ["GENE1", "GENE2", "GENE3"]
enriched = go.enrich_genes(gene_list, background_genes)

# Semantic similarity calculation
similarity = go.semantic_similarity("GO:0008150", "GO:0009987")
```

## Integration

Ontology tools integrate with:
- **DNA/RNA analysis** for functional annotation
- **Protein analysis** for domain and family classification
- **Network analysis** for functional enrichment
- **Statistical methods** for enrichment analysis

## Testing

Comprehensive tests ensure annotation reliability:
- Ontology parsing and validation
- Enrichment analysis correctness
- Semantic similarity validation
- Integration with external databases

## Contributing

When adding new ontology functionality:
1. Update annotation and enrichment documentation
2. Add comprehensive validation tests
3. Ensure compatibility with OBO format standards
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's ontology capabilities.
