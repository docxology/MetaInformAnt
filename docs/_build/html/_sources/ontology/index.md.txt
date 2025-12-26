# Ontology Module Overview

Functional annotation and ontology analysis utilities for biological ontologies like Gene Ontology (GO).

## Features

- **OBO Format Parsing**: Parse OBO format ontology files with support for:
  - Basic fields (id, name, namespace, definition)
  - Alternative IDs
  - Multiple relationship types (is_a, part_of, regulates, etc.)
  - Synonyms and cross-references
  - GO subsets

- **Hierarchy Traversal**: 
  - Ancestor and descendant queries
  - Common ancestor finding
  - Path finding and distance calculations
  - Subgraph extraction

- **Query Functions**:
  - Term lookup by name
  - Namespace filtering
  - Root and leaf term identification

- **Validation**: 
  - Ontology integrity checks
  - Cycle detection
  - Orphaned term detection

- **Serialization**: 
  - Save/load ontologies to/from JSON

- **Performance**: 
  - In-memory caching for expensive operations
  - Efficient BFS traversal algorithms

## Documentation

- [Gene Ontology (GO) Analysis](./go.md) - GO-specific functions and usage

## Module Reference

See `src/metainformant/ontology/README.md` for comprehensive API documentation and examples.
