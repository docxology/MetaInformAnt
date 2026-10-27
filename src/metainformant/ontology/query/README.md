# Ontology Query

Querying, traversal, and serialization of ontology data structures. Provides ancestor/descendant discovery, path finding, subontology extraction, and file I/O in OBO and JSON formats.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports the `query` and `serialize` modules and every public function of both, so `from metainformant.ontology.query import ancestors` works |
| `query.py` | Ontology graph traversal, search, statistics, and validation |
| `serialize.py` | Save/load ontologies in OBO, JSON, and NetworkX graph formats |

## Key Functions

| Function | Description |
|----------|-------------|
| `query.ancestors()` | Find all ancestors of a term via a relationship type (the term itself is excluded) |
| `query.descendants()` | Find all descendants of a term (the term itself is excluded) |
| `query.common_ancestors()` | Compute shared ancestors between two terms |
| `query.shortest_path()` | Find shortest path between two terms |
| `query.get_roots()` | Identify root terms with no parents |
| `query.get_leaves()` | Identify leaf terms with no children |
| `query.information_content()` | Compute information content of a term |
| `query.find_terms_by_name()` | Search terms by name pattern |
| `query.filter_by_namespace()` | New ontology with one namespace's terms and the parent/child links among them (`get_roots`, `get_parents`, `ancestors` work on it) |
| `query.subgraph()` | New ontology with the given terms and, for `is_a`, the parent/child links among them |
| `query.get_subontology()` | New ontology with the given terms and all their descendants, with the induced parent/child links (a term whose parents were left out is a root) |
| `query.clear_cache()` | Empty the `ancestors`/`descendants` cache (needed after editing `parents_of`/`children_of` in place) |
| `query.validate_ontology_integrity()` | Check ontology for structural issues |
| `serialize.save_ontology()` | Save ontology to OBO or JSON file (OBO definitions and synonyms are written escaped, definitions followed by their dbxref list) |
| `serialize.load_ontology()` | Load ontology from OBO or JSON file; a JSON file saved from an OBO-parsed ontology loads with every `is_a` link once |
| `serialize.ontology_to_graph()` | Convert ontology to NetworkX graph |
| `serialize.merge_ontologies()` | Merge multiple ontologies with conflict resolution |

## Usage

```python
from metainformant.ontology.query import query, serialize

anc = query.ancestors(ontology, "GO:0008150")
serialize.save_ontology(ontology, "output/go.json", format="json")
```

The `ancestors` and `descendants` results are cached per ontology object. An entry holds a weak reference to its
ontology, so it is dropped with the ontology, and it is refreshed when `add_term` or `add_relationship` changes the
ontology's structure; only an in-place edit of the parent or child sets needs `clear_cache()`.

Every public function is also importable from the package itself:

```python
from metainformant.ontology.query import ancestors, descendants, load_ontology, save_ontology

anc = ancestors(ontology, "GO:0008150")
```
