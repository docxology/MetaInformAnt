# Ontology Querying and Traversal

Functions for querying, traversing, and analysing ontology graphs: ancestor
and descendant retrieval, common ancestors, shortest paths, sub-ontology
extraction, term search, information content, integrity validation, and
caching.

## Key Concepts

**Graph traversal** uses breadth-first search over `Relationship` edges.
All traversal functions accept a `relation_type` parameter (default `"is_a"`)
to restrict which edge types are followed.

**Information content** (IC) measures term specificity. Computed as
`-log(|descendants(t)| / |all_terms|)`. Higher IC indicates a more specific
(less frequently annotated) term.

**Caching** is enabled by default. Query results for `ancestors` and
`descendants` are stored in a module-level dictionary. Use `clear_cache()` and
`set_cache_enabled()` to manage.

## Function Reference

### ancestors

```python
def ancestors(
    onto: Ontology, term_id: str, relation_type: str = "is_a",
) -> Set[str]
```

Returns all ancestors of a term (including the term itself) reachable via the
specified relationship type. Raises `TermNotFoundError` if term_id is absent.

### descendants

```python
def descendants(
    onto: Ontology, term_id: str, relation_type: str = "is_a",
) -> Set[str]
```

Returns all descendants of a term (including the term itself).

### common_ancestors

```python
def common_ancestors(
    onto: Ontology, term1: str, term2: str, relation_type: str = "is_a",
) -> Set[str]
```

Intersection of ancestor sets for two terms.

### most_informative_common_ancestor

```python
def most_informative_common_ancestor(
    onto: Ontology, term1: str, term2: str,
    ic_map: Dict[str, float], relation_type: str = "is_a",
) -> str
```

Returns the common ancestor with the highest information content. Raises
`ValueError` when no common ancestors exist.

### path_to_root

```python
def path_to_root(
    onto: Ontology, term_id: str, relation_type: str = "is_a",
) -> List[str]
```

Returns one path from a term to a root (term with no parents). Follows the
first parent at each step.

### shortest_path

```python
def shortest_path(
    onto: Ontology, term1: str, term2: str, relation_type: str = "is_a",
) -> List[str]
```

BFS-based undirected shortest path between two terms. Returns an empty list
if no path exists.

### distance

```python
def distance(
    onto: Ontology, term1: str, term2: str, relation_type: str = "is_a",
) -> int
```

Minimum edge count between two terms, or -1 if not connected.

### get_subontology

```python
def get_subontology(
    onto: Ontology, root_terms: Iterable[str], relation_type: str = "is_a",
) -> Ontology
```

Extracts a sub-ontology containing all descendants of the specified root
terms and the relationships between them.

### subgraph

```python
def subgraph(
    onto: Ontology, term_ids: List[str], relation_type: str = "is_a",
) -> Ontology
```

Creates a new ontology containing only the specified terms and the
relationships that connect them.

### find_terms_by_name / find_terms_by_namespace / filter_by_namespace

```python
def find_terms_by_name(onto: Ontology, name_pattern: str, case_sensitive: bool = False) -> List[str]
def find_terms_by_namespace(onto: Ontology, namespace: str) -> List[str]
def filter_by_namespace(onto: Ontology, namespace: str) -> Ontology
```

`find_terms_by_name` does substring search against term names.
`find_terms_by_namespace` returns term IDs in a namespace.
`filter_by_namespace` returns a new `Ontology` with only that namespace.

### get_roots / get_leaves

```python
def get_roots(onto: Ontology, relation_type: str = "is_a") -> Set[str]
def get_leaves(onto: Ontology, relation_type: str = "is_a") -> Set[str]
```

Root terms have no parents; leaf terms have no children.

### information_content / calculate_ic_map

```python
def information_content(onto: Ontology, term_id: str, corpus_size: int | None = None) -> float
def calculate_ic_map(onto: Ontology, corpus_size: int | None = None) -> Dict[str, float]
```

IC for a single term or for all terms in the ontology.

### get_subontology_stats

```python
def get_subontology_stats(onto: Ontology) -> Dict[str, Any]
```

Comprehensive statistics: term/relationship/obsolete counts, namespace
distribution, relationship type distribution, depth statistics, root/leaf
counts, and definition coverage.

### validate_ontology_integrity

```python
def validate_ontology_integrity(onto: Ontology) -> Tuple[bool, List[str]]
```

Checks for dangling references (relationships pointing to missing terms),
self-loops, and duplicate relationships.

### Cache Management

```python
def clear_cache() -> None
def set_cache_enabled(enabled: bool) -> None
```

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
