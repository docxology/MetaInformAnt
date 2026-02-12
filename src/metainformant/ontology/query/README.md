# Ontology Query

Querying, traversal, and serialization of ontology data structures. Provides ancestor/descendant discovery, path finding, subontology extraction, and file I/O in OBO and JSON formats.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `query` and `serialize` modules |
| `query.py` | Ontology graph traversal, search, statistics, and validation |
| `serialize.py` | Save/load ontologies in OBO, JSON, and NetworkX graph formats |

## Key Functions

| Function | Description |
|----------|-------------|
| `query.ancestors()` | Find all ancestors of a term via a relationship type |
| `query.descendants()` | Find all descendants of a term |
| `query.common_ancestors()` | Compute shared ancestors between two terms |
| `query.shortest_path()` | Find shortest path between two terms |
| `query.get_roots()` | Identify root terms with no parents |
| `query.get_leaves()` | Identify leaf terms with no children |
| `query.information_content()` | Compute information content of a term |
| `query.find_terms_by_name()` | Search terms by name pattern |
| `query.validate_ontology_integrity()` | Check ontology for structural issues |
| `serialize.save_ontology()` | Save ontology to OBO or JSON file |
| `serialize.load_ontology()` | Load ontology from OBO or JSON file |
| `serialize.ontology_to_graph()` | Convert ontology to NetworkX graph |
| `serialize.merge_ontologies()` | Merge multiple ontologies with conflict resolution |

## Usage

```python
from metainformant.ontology.query import query, serialize

anc = query.ancestors(ontology, "GO:0008150")
serialize.save_ontology(ontology, "output/go.json", format="json")
```
