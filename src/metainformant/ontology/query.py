from __future__ import annotations

import time
from collections import deque
from typing import Iterable, List, Set

from metainformant.core.logging import get_logger

from .types import Ontology

logger = get_logger(__name__)

# In-memory cache for traversal operations
_cache: dict[str, tuple[Set[str], float]] = {}
_cache_ttl: float = 3600.0  # 1 hour default TTL
_cache_enabled: bool = True


def _cache_key(term_id: str, direction: str, ontology_id: str) -> str:
    """Generate cache key for traversal operation."""
    return f"{ontology_id}:{direction}:{term_id}"


def _get_cached(term_id: str, direction: str, ontology_id: str) -> Set[str] | None:
    """Get cached result if available and not expired."""
    if not _cache_enabled:
        return None
    
    key = _cache_key(term_id, direction, ontology_id)
    if key not in _cache:
        return None
    
    result, timestamp = _cache[key]
    if time.time() - timestamp > _cache_ttl:
        del _cache[key]
        return None
    
    return result


def _set_cached(term_id: str, direction: str, ontology_id: str, result: Set[str]) -> None:
    """Cache traversal result."""
    if not _cache_enabled:
        return
    
    key = _cache_key(term_id, direction, ontology_id)
    _cache[key] = (result, time.time())


def clear_cache() -> None:
    """Clear all cached traversal results."""
    global _cache
    _cache.clear()
    logger.debug("Cleared ontology traversal cache")


def set_cache_enabled(enabled: bool) -> None:
    """Enable or disable caching for traversal operations.
    
    Args:
        enabled: If True, enable caching; if False, disable and clear cache
    """
    global _cache_enabled
    _cache_enabled = enabled
    if not enabled:
        clear_cache()


def set_cache_ttl(seconds: float) -> None:
    """Set cache time-to-live in seconds.
    
    Args:
        seconds: TTL in seconds (0 or negative to disable TTL)
    """
    global _cache_ttl
    _cache_ttl = seconds


def ancestors(onto: Ontology, term_id: str, use_cache: bool = True) -> Set[str]:
    """Get all ancestor terms (transitive parents) of a given term.
    
    Returns the complete set of parent terms reachable via is_a relationships,
    excluding the term itself. Useful for finding all broader terms.
    
    Args:
        onto: Ontology object containing terms
        term_id: Identifier of the term (e.g., "GO:0008150")
        use_cache: If True, use cached result if available
        
    Returns:
        Set of ancestor term IDs. Returns empty set if term not found or has no parents.
        
    Raises:
        ValueError: If term_id is empty or not found in ontology
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> ancestors(onto, "GO:0008150")  # biological_process
        {'GO:0008150', 'GO:0003674', ...}
    """
    if not term_id:
        raise ValueError("term_id cannot be empty")
    
    if not onto.has_term(term_id):
        raise ValueError(f"Term {term_id} not found in ontology")
    
    # Try cache first
    ontology_id = str(id(onto))  # Use object ID as ontology identifier
    if use_cache:
        cached = _get_cached(term_id, "ancestors", ontology_id)
        if cached is not None:
            return cached
    
    visited: Set[str] = set()
    queue: deque[str] = deque(onto.parents_of.get(term_id, set()))
    
    while queue:
        node = queue.popleft()
        if node in visited:
            continue
        visited.add(node)
        queue.extend(onto.parents_of.get(node, set()))
    
    if len(visited) > 1000:
        logger.debug(f"Large ancestor set for {term_id}: {len(visited)} terms")
    
    # Cache result
    if use_cache:
        _set_cached(term_id, "ancestors", ontology_id, visited)
    
    return visited


def descendants(onto: Ontology, term_id: str, use_cache: bool = True) -> Set[str]:
    """Get all descendant terms (transitive children) of a given term.
    
    Returns the complete set of child terms reachable via is_a relationships,
    excluding the term itself. Useful for finding all more specific terms.
    
    Args:
        onto: Ontology object containing terms
        term_id: Identifier of the term (e.g., "GO:0008150")
        use_cache: If True, use cached result if available
        
    Returns:
        Set of descendant term IDs. Returns empty set if term not found or has no children.
        
    Raises:
        ValueError: If term_id is empty or not found in ontology
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> descs = descendants(onto, "GO:0008150")  # biological_process
        >>> len(descs) > 0
        True
    """
    if not term_id:
        raise ValueError("term_id cannot be empty")
    
    if not onto.has_term(term_id):
        raise ValueError(f"Term {term_id} not found in ontology")
    
    # Try cache first
    ontology_id = str(id(onto))
    if use_cache:
        cached = _get_cached(term_id, "descendants", ontology_id)
        if cached is not None:
            return cached
    
    visited: Set[str] = set()
    queue: deque[str] = deque(onto.children_of.get(term_id, set()))
    
    while queue:
        node = queue.popleft()
        if node in visited:
            continue
        visited.add(node)
        queue.extend(onto.children_of.get(node, set()))
    
    if len(visited) > 1000:
        logger.debug(f"Large descendant set for {term_id}: {len(visited)} terms")
    
    # Cache result
    if use_cache:
        _set_cached(term_id, "descendants", ontology_id, visited)
    
    return visited


def subgraph(onto: Ontology, roots: Iterable[str]) -> Ontology:
    """Extract subgraph of ontology rooted at specified terms.
    
    Creates a new ontology containing only the specified root terms and
    all their descendants (transitively). Useful for focusing analysis
    on specific ontology branches.
    
    Args:
        onto: Source ontology object
        roots: Iterable of term IDs to use as root nodes. All descendants
            of these terms will be included in the subgraph.
            
    Returns:
        New Ontology object containing only the root terms and their
        descendants. Term relationships are preserved within the subgraph.
        
    Raises:
        ValueError: If any root term_id is empty or not found in ontology
        
    Examples:
        >>> onto = load_go_obo("go-basic.obo")
        >>> sub_onto = subgraph(onto, ["GO:0008150"])  # biological_process
        >>> sub_onto.num_terms() < onto.num_terms()
        True
        >>> "GO:0008150" in sub_onto.terms
        True
    """
    keep: Set[str] = set()
    for r in roots:
        if not r:
            raise ValueError("Root term_id cannot be empty")
        if not onto.has_term(r):
            raise ValueError(f"Root term {r} not found in ontology")
        keep.add(r)
        keep.update(descendants(onto, r))

    new = Ontology()
    for tid, term in onto.terms.items():
        if tid in keep:
            parents = [p for p in term.is_a_parents if p in keep]
            # Filter relationships to only include terms in keep set
            filtered_rels = {
                rel_type: [t for t in rel_terms if t in keep]
                for rel_type, rel_terms in term.relationships.items()
            }
            new.add_term(
                type(term)(
                    term_id=term.term_id,
                    name=term.name,
                    namespace=term.namespace,
                    definition=term.definition,
                    alt_ids=list(term.alt_ids),
                    is_a_parents=parents,
                    relationships=filtered_rels,
                    synonyms=list(term.synonyms),
                    xrefs=list(term.xrefs),
                    subsets=list(term.subsets),
                )
            )
    return new


def common_ancestors(onto: Ontology, term1: str, term2: str) -> Set[str]:
    """Get common ancestor terms of two terms.
    
    Finds all terms that are ancestors of both term1 and term2.
    
    Args:
        onto: Ontology object containing terms
        term1: First term identifier
        term2: Second term identifier
        
    Returns:
        Set of common ancestor term IDs
        
    Raises:
        ValueError: If either term_id is empty or not found
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> common = common_ancestors(onto, "GO:0009987", "GO:0008150")
        >>> len(common) > 0
        True
    """
    anc1 = ancestors(onto, term1)
    anc2 = ancestors(onto, term2)
    # Also include the terms themselves as potential common ancestors
    anc1.add(term1)
    anc2.add(term2)
    return anc1 & anc2


def path_to_root(onto: Ontology, term_id: str) -> List[str]:
    """Get shortest path from term to root (term with no parents).
    
    Returns a list of term IDs representing the path from the given term
    to a root term, following parent relationships. If multiple paths exist,
    returns one of the shortest paths.
    
    Args:
        onto: Ontology object containing terms
        term_id: Starting term identifier
        
    Returns:
        List of term IDs from term_id to a root, including term_id itself.
        Returns [term_id] if term has no parents.
        
    Raises:
        ValueError: If term_id is empty or not found
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> path = path_to_root(onto, "GO:0009987")
        >>> len(path) >= 1
        True
        >>> path[0] == "GO:0009987"
        True
    """
    if not term_id:
        raise ValueError("term_id cannot be empty")
    
    if not onto.has_term(term_id):
        raise ValueError(f"Term {term_id} not found in ontology")
    
    path: List[str] = [term_id]
    current = term_id
    
    # Follow parents until we reach a term with no parents
    while onto.parents_of.get(current):
        # Get first parent (arbitrary choice if multiple parents)
        parent = next(iter(onto.parents_of[current]))
        if parent in path:  # Cycle detection
            logger.warning(f"Cycle detected in path_to_root for {term_id}")
            break
        path.append(parent)
        current = parent
    
    return path


def distance(onto: Ontology, term1: str, term2: str) -> int | None:
    """Calculate shortest path distance between two terms.
    
    Uses common ancestors to find shortest path. Distance is the number of
    edges in the shortest path, or None if no path exists.
    
    Args:
        onto: Ontology object containing terms
        term1: First term identifier
        term2: Second term identifier
        
    Returns:
        Integer distance (number of edges), or None if no path exists
        
    Raises:
        ValueError: If either term_id is empty or not found
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> dist = distance(onto, "GO:0009987", "GO:0008150")
        >>> dist is not None
        True
    """
    if not term1 or not term2:
        raise ValueError("term_ids cannot be empty")
    
    if not onto.has_term(term1):
        raise ValueError(f"Term {term1} not found in ontology")
    if not onto.has_term(term2):
        raise ValueError(f"Term {term2} not found in ontology")
    
    if term1 == term2:
        return 0
    
    # Find common ancestors
    common = common_ancestors(onto, term1, term2)
    if not common:
        return None
    
    # Find lowest common ancestor (LCA) - term with maximum depth
    # We approximate by finding the common ancestor with the longest path from root
    lca = None
    max_depth = -1
    
    for ca in common:
        path = path_to_root(onto, ca)
        depth = len(path) - 1
        if depth > max_depth:
            max_depth = depth
            lca = ca
    
    if lca is None:
        return None
    
    # Calculate distance: path from term1 to LCA + path from term2 to LCA
    path1 = path_to_root(onto, term1)
    path2 = path_to_root(onto, term2)
    
    # Find positions of LCA in each path
    try:
        idx1 = path1.index(lca)
        idx2 = path2.index(lca)
        return idx1 + idx2
    except ValueError:
        return None


def find_term_by_name(onto: Ontology, name: str, namespace: str | None = None) -> List[str]:
    """Find term IDs by name (case-insensitive partial match).
    
    Searches for terms matching the given name. Returns all matching term IDs.
    
    Args:
        onto: Ontology object containing terms
        name: Term name to search for (case-insensitive substring match)
        namespace: Optional namespace filter
        
    Returns:
        List of matching term IDs
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> matches = find_term_by_name(onto, "biological process")
        >>> len(matches) > 0
        True
    """
    if not name:
        return []
    
    name_lower = name.lower()
    matches: List[str] = []
    
    for term_id, term in onto.terms.items():
        if namespace and term.namespace != namespace:
            continue
        if name_lower in term.name.lower():
            matches.append(term_id)
    
    return matches


def filter_by_namespace(onto: Ontology, namespace: str) -> Ontology:
    """Filter ontology to include only terms in specified namespace.
    
    Creates a new ontology containing only terms from the given namespace,
    preserving relationships within that namespace.
    
    Args:
        onto: Source ontology object
        namespace: Namespace to filter by (e.g., "biological_process")
        
    Returns:
        New Ontology object containing only terms from the namespace
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> bp_onto = filter_by_namespace(onto, "biological_process")
        >>> bp_onto.num_terms() < onto.num_terms()
        True
    """
    if not namespace:
        raise ValueError("namespace cannot be empty")
    
    new = Ontology()
    for term_id, term in onto.terms.items():
        if term.namespace == namespace:
            # Only include parents that are also in this namespace
            parents = [p for p in term.is_a_parents if onto.get_namespace(p) == namespace]
            # Filter relationships to only include terms in same namespace
            filtered_rels = {
                rel_type: [t for t in rel_terms if onto.get_namespace(t) == namespace]
                for rel_type, rel_terms in term.relationships.items()
            }
            new.add_term(
                type(term)(
                    term_id=term.term_id,
                    name=term.name,
                    namespace=term.namespace,
                    definition=term.definition,
                    alt_ids=list(term.alt_ids),
                    is_a_parents=parents,
                    relationships=filtered_rels,
                    synonyms=list(term.synonyms),
                    xrefs=list(term.xrefs),
                    subsets=list(term.subsets),
                )
            )
    return new


def get_roots(onto: Ontology, namespace: str | None = None) -> Set[str]:
    """Get root terms (terms with no parents).
    
    Args:
        onto: Ontology object containing terms
        namespace: Optional namespace filter
        
    Returns:
        Set of root term IDs (terms with no parents)
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> roots = get_roots(onto)
        >>> len(roots) > 0
        True
    """
    roots: Set[str] = set()
    
    for term_id, term in onto.terms.items():
        if namespace and term.namespace != namespace:
            continue
        if not onto.parents_of.get(term_id):
            roots.add(term_id)
    
    return roots


def get_leaves(onto: Ontology, namespace: str | None = None) -> Set[str]:
    """Get leaf terms (terms with no children).
    
    Args:
        onto: Ontology object containing terms
        namespace: Optional namespace filter
        
    Returns:
        Set of leaf term IDs (terms with no children)
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> leaves = get_leaves(onto)
        >>> len(leaves) > 0
        True
    """
    leaves: Set[str] = set()
    
    for term_id, term in onto.terms.items():
        if namespace and term.namespace != namespace:
            continue
        if not onto.children_of.get(term_id):
            leaves.add(term_id)
    
    return leaves
