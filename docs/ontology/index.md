# Ontology Module Overview

Functional annotation and ontology analysis utilities for biological ontologies like Gene Ontology (GO).

## Features

- **OBO Format Parsing**: Parse OBO format ontology files with support for:
  - Basic fields (id, name, namespace, definition)
  - Alternative IDs
  - Multiple relationship types (is_a, part_of, regulates, etc.)
  - Synonyms (quoted text, scope and dbxrefs; escapes round-trip) and cross-references
  - GO subsets

- **Hierarchy Traversal**: 
  - Ancestor and descendant queries
  - Common ancestor finding
  - Path finding and distance calculations
  - Subgraph, sub-ontology and namespace extraction (the extracted ontology keeps the parent/child links among its terms)

- **Query Functions**:
  - Term lookup by name
  - Namespace filtering
  - Root and leaf term identification

- **Validation**: 
  - Ontology integrity checks
  - Cycle detection
  - Orphaned term detection

- **Serialization**: 
  - Save/load ontologies to/from OBO and JSON (ontologies parsed from OBO round-trip through JSON)

- **Enrichment**:
  - GO over-representation analysis (`enrich_genes`) and information-content semantic similarity
  - Pathway ORA, preranked GSEA with a seeded gene-set permutation null, and pathway similarity networks

- **Performance**: 
  - Per-ontology caching of ancestor and descendant queries (refreshed when terms or relationships are added)
  - Efficient BFS traversal algorithms

## Documentation

- [Gene Ontology (GO) Analysis](./go.md) - GO-specific functions and usage
- [OBO Format Parsing](./obo.md) - OBO parser, definitions and synonyms
- [Ontology Querying and Traversal](./query.md) - traversal, sub-ontologies, caching
- [Pathway Enrichment Analysis](./pathway_enrichment.md) - ORA, GSEA, pathway networks

## Module Reference

See `src/metainformant/ontology/README.md` for comprehensive API documentation and examples.
