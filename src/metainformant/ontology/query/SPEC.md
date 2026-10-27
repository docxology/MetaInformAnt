# Specification: query

## 🎯 Scope
Ontology querying, traversal, and serialization.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `query.py`
- `serialize.py`

### Behavior Contracts

#### `query.py`
- Query cache — `ancestors` and `descendants` cache their results per ontology object. An entry holds a weak reference
  to its ontology (it disappears with it and is never served to another object that reuses the address) and the
  structure stamp `(len(terms), len(relationships), len(parents_of), len(children_of))`; a different stamp discards the
  entry, so `Ontology.add_term` and `Ontology.add_relationship` refresh the cache. An in-place edit of a parent or child
  set that changes no container size is not detected and needs `clear_cache()`.
- `filter_by_namespace`, `subgraph` (for `relation_type="is_a"`) and `get_subontology` (for `relation_type="is_a"`)
  return an ontology whose `parents_of` / `children_of` hold the `is_a` links with both ends among the kept terms (plus
  kept `is_a` relationships), so `get_roots`, `get_parents`, `get_children` and `get_ancestors` work on the result. A
  term whose parents were all left out has no entry and is a root. `get_subontology` reads `root_terms` once, so any
  iterable (a generator included) is accepted.

#### `serialize.py`
- `load_ontology(path, format="obo" | "json")` — the JSON loader rebuilds each term field by field (`Term.metadata` is
  restored as saved, including the `relationships` list of OBO-parsed terms) and merges the saved relationship list
  with the `is_a` relationships that `add_term` creates, so every `is_a` link comes back once, with its saved metadata.
- `save_ontology(onto, path, format="obo")` writes synonyms as `synonym: "<escaped text>" EXACT []`.
