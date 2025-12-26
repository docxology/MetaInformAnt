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

The ontology domain supports ontology parsing, traversal, and analysis:

```python
from metainformant.ontology import (
    load_go_obo, ancestors, descendants,
    common_ancestors, distance, filter_by_namespace
)

# Load Gene Ontology
onto = load_go_obo("go-basic.obo")

# Traverse hierarchy
ancestors_set = ancestors(onto, "GO:0008150")
descendants_set = descendants(onto, "GO:0008150")

# Find relationships
common = common_ancestors(onto, "GO:0009987", "GO:0008150")
dist = distance(onto, "GO:0009987", "GO:0008150")

# Filter by namespace
bp_onto = filter_by_namespace(onto, "biological_process")
```

**Note**: Enrichment analysis and semantic similarity functions are planned but not yet implemented.

## Integration

Ontology tools integrate with:
- **DNA/RNA analysis** for functional annotation
- **Protein analysis** for domain and family classification
- **Network analysis** for functional enrichment
- **Statistical methods** for enrichment analysis

## Testing

Comprehensive tests ensure reliability:
- Ontology parsing and validation
- Hierarchy traversal correctness
- Error handling and edge cases
- Serialization integrity

## Contributing

When adding new ontology functionality:
1. Update annotation and enrichment documentation
2. Add comprehensive validation tests
3. Ensure compatibility with OBO format standards
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's ontology capabilities.
