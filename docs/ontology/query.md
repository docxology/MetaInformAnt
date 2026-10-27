# Ontology Querying and Traversal

Functions for querying, traversing, and analysing ontology graphs: ancestor
and descendant retrieval, common ancestors, shortest paths, sub-ontology
extraction, term search, information content, integrity validation, and
caching.

## Key Concepts

**Graph traversal.** `ancestors` and `descendants` collect the transitive
closure of the parent or child links, while `path_to_root`, `shortest_path`
and `distance` use breadth-first search. The traversal functions accept a
`relation_type` parameter (default `"is_a"`) to restrict which edge types are
followed. `get_roots` and `get_leaves` accept it too but always use the
`is_a` parent and child links.

**Information content** (IC) measures term specificity. Computed as
`-log(|descendants(t)| / |all_terms|)`, where `descendants(t)` excludes `t`
itself and `corpus_size` replaces `|all_terms|` when given. Among terms that
have descendants, fewer descendants give a higher IC; a term with no
descendants (a leaf) has IC `0.0`.

**Caching** is enabled by default. Query results for `ancestors` and
`descendants` are cached per ontology object. A cache entry holds a weak
reference to its ontology, so it is dropped together with the ontology and can
never be served to a different ontology that reuses its memory address. An
entry is also discarded when the ontology's structure changes: `add_term` and
`add_relationship` refresh it automatically, because the check compares the
number of terms, relationships, and parent and child entries. Editing a parent
or child set in place without adding or removing an entry is not detected, so
call `clear_cache()` after such an edit. Use `clear_cache()` and
`set_cache_enabled()` to manage the cache.

## Function Reference

### ancestors

```python-snippet
def ancestors(
    onto: Ontology, term_id: str, relation_type: str = "is_a",
) -> Set[str]
```

Returns all ancestors of a term reachable via the specified relationship type.
The term itself is not included: a term is not its own ancestor, so a root
term has no ancestors (in a cyclic ontology, a term on the cycle is reachable
from itself). Use `ancestors(onto, term_id) | {term_id}` when the term itself
is wanted. Raises `TermNotFoundError` if term_id is absent.

### descendants

```python-snippet
def descendants(
    onto: Ontology, term_id: str, relation_type: str = "is_a",
) -> Set[str]
```

Returns all descendants of a term reachable via the specified relationship
type. The term itself is not included, so a leaf term has no descendants; use
`descendants(onto, term_id) | {term_id}` when the term itself is wanted. Raises
`TermNotFoundError` if term_id is absent.

### common_ancestors

```python-snippet
def common_ancestors(
    onto: Ontology, term1: str, term2: str, relation_type: str = "is_a",
) -> Set[str]
```

Intersection of the ancestor sets of two terms. Because a term is not its own
ancestor, neither term appears in the result: for a term and its parent, the
result holds only the parent's ancestors.

### most_informative_common_ancestor

```python-snippet
def most_informative_common_ancestor(
    onto: Ontology, term1: str, term2: str,
    ic_map: Dict[str, float], relation_type: str = "is_a",
) -> str
```

Returns the common ancestor with the highest information content. Raises
`ValueError` when no common ancestors exist, for example when one of the terms
is a root (a root has no ancestors, and neither term counts as its own
ancestor).

### path_to_root

```python-snippet
def path_to_root(
    onto: Ontology, term_id: str, relation_type: str = "is_a",
) -> List[str]
```

Returns a shortest path (fewest edges) from a term to a root (a term with no
parents), found by breadth-first search, with the term first and the root last.
A root returns `[term_id]`. When several roots are equally close, which path is
returned is not specified. Returns an empty list when no root is reachable.

### shortest_path

```python-snippet
def shortest_path(
    onto: Ontology, term1: str, term2: str, relation_type: str = "is_a",
) -> List[str]
```

BFS-based undirected shortest path between two terms. Returns an empty list
if no path exists.

### distance

```python-snippet
def distance(
    onto: Ontology, term1: str, term2: str, relation_type: str = "is_a",
) -> int
```

Minimum edge count between two terms, or -1 if not connected.

### get_subontology

```python-snippet
def get_subontology(
    onto: Ontology, root_terms: Iterable[str], relation_type: str = "is_a",
) -> Ontology
```

Extracts a sub-ontology containing the specified root terms, all of their
descendants, and the relationships of type `relation_type` between them. For
`relation_type="is_a"` the result also carries the parent and child links among
its terms, so `get_roots`, `get_parents`, `get_children` and `get_ancestors`
work on it. Parents outside the sub-ontology are not carried over: a root term
whose own parents were left out is a root of the sub-ontology, and
`get_ancestors` stays inside it. `root_terms` can be any iterable, a generator
included. Raises `TermNotFoundError` if a root term is absent.

### subgraph

```python-snippet
def subgraph(
    onto: Ontology, term_ids: List[str], relation_type: str = "is_a",
) -> Ontology
```

Creates a new ontology containing only the specified terms (IDs missing from
`onto` are skipped) and the `relation_type` relationships that connect them.
For `relation_type="is_a"` the result also carries the parent and child links
among those terms, so `get_roots`, `get_parents`, `get_children`, `ancestors`,
`descendants` and `path_to_root` work on it. Only links with both ends in the
subgraph are kept: a term whose parents were all left out is a root of the
subgraph. For any other `relation_type` the result has no parent and child
indexes, so the `is_a` queries find no links on it, while `shortest_path` and
`distance` read the relationship list and still work.

### find_terms_by_name / find_terms_by_namespace / filter_by_namespace

```python-snippet
def find_terms_by_name(onto: Ontology, name_pattern: str, case_sensitive: bool = False) -> List[str]
def find_terms_by_namespace(onto: Ontology, namespace: str) -> List[str]
def filter_by_namespace(onto: Ontology, namespace: str) -> Ontology
```

`find_terms_by_name` does substring search against term names.
`find_terms_by_namespace` returns term IDs in a namespace.
`filter_by_namespace` returns a new `Ontology` with only that namespace's terms,
the relationships among them, and the parent and child links among them, so
`get_roots`, `get_parents`, `get_children` and `ancestors` work on it. A term
whose parents all belong to other namespaces becomes a root of the filtered
ontology.

### get_roots / get_leaves

```python-snippet
def get_roots(onto: Ontology, relation_type: str = "is_a") -> Set[str]
def get_leaves(onto: Ontology, relation_type: str = "is_a") -> Set[str]
```

Root terms have no parents; leaf terms have no children. `relation_type` is
accepted, but the `is_a` parent and child indexes are always used.

### information_content / calculate_ic_map

```python-snippet
def information_content(onto: Ontology, term_id: str, corpus_size: int | None = None) -> float
def calculate_ic_map(onto: Ontology, corpus_size: int | None = None) -> Dict[str, float]
```

IC for a single term or for all terms in the ontology.

### get_subontology_stats

```python-snippet
def get_subontology_stats(onto: Ontology) -> Dict[str, Any]
```

Comprehensive statistics: term/relationship/obsolete counts, namespace
distribution, relationship type distribution, depth statistics, root/leaf
counts, and definition coverage.

### validate_ontology_integrity

```python-snippet
def validate_ontology_integrity(onto: Ontology) -> Tuple[bool, List[str]]
```

Checks for dangling references (relationships pointing to missing terms),
self-loops, and duplicate relationships.

### Cache Management

```python-snippet
def clear_cache() -> None
def set_cache_enabled(enabled: bool) -> None
```

`clear_cache` empties the cache of every ontology; it is needed after editing
`parents_of` or `children_of` in place (see Caching above). Adding terms or
relationships with `add_term` and `add_relationship`, and garbage-collecting an
ontology, need no call. `set_cache_enabled(False)` turns caching off.

## Usage Example

```python
from metainformant.ontology.core.obo import parse_obo
from metainformant.ontology.query.query import (
    ancestors, descendants, shortest_path, calculate_ic_map,
    most_informative_common_ancestor, get_subontology_stats,
)

ontology = parse_obo("data/go-basic.obo")

# Ancestor traversal
anc = ancestors(ontology, "GO:0006915")  # apoptotic process
print(f"Ancestors: {len(anc)}")

# Semantic similarity via MICA
ic_map = calculate_ic_map(ontology)
mica = most_informative_common_ancestor(ontology, "GO:0006915", "GO:0006917", ic_map)

# Shortest path
path = shortest_path(ontology, "GO:0006915", "GO:0008150")
print(f"Path length: {len(path) - 1}")

# Statistics
stats = get_subontology_stats(ontology)
print(f"Terms: {stats['num_terms']}, Max depth: {stats['depth_stats']['max']}")
```

## Related Modules

- `metainformant.ontology.core.obo` -- OBO file parsing
- `metainformant.ontology.core.types` -- data structures
- `metainformant.ontology.pathway_enrichment` -- enrichment analysis
- `metainformant.ontology.visualization.visualization` -- ontology plotting
