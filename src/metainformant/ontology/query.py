from __future__ import annotations

from collections import deque
from typing import Iterable, List, Set

from .types import Ontology


def ancestors(onto: Ontology, term_id: str) -> Set[str]:
    """Get all ancestor terms (transitive parents) of a given term.
    
    Returns the complete set of parent terms reachable via is_a relationships,
    excluding the term itself. Useful for finding all broader terms.
    
    Args:
        onto: Ontology object containing terms
        term_id: Identifier of the term (e.g., "GO:0008150")
        
    Returns:
        Set of ancestor term IDs. Returns empty set if term not found or has no parents.
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> ancestors(onto, "GO:0008150")  # biological_process
        {'GO:0008150', 'GO:0003674', ...}
    """
    visited: Set[str] = set()
    queue: deque[str] = deque(onto.parents_of.get(term_id, set()))
    while queue:
        node = queue.popleft()
        if node in visited:
            continue
        visited.add(node)
        queue.extend(onto.parents_of.get(node, set()))
    return visited


def descendants(onto: Ontology, term_id: str) -> Set[str]:
    """Get all descendant terms (transitive children) of a given term.
    
    Returns the complete set of child terms reachable via is_a relationships,
    excluding the term itself. Useful for finding all more specific terms.
    
    Args:
        onto: Ontology object containing terms
        term_id: Identifier of the term (e.g., "GO:0008150")
        
    Returns:
        Set of descendant term IDs. Returns empty set if term not found or has no children.
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> descs = descendants(onto, "GO:0008150")  # biological_process
        >>> len(descs) > 0
        True
    """
    visited: Set[str] = set()
    queue: deque[str] = deque(onto.children_of.get(term_id, set()))
    while queue:
        node = queue.popleft()
        if node in visited:
            continue
        visited.add(node)
        queue.extend(onto.children_of.get(node, set()))
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
        keep.add(r)
        keep.update(descendants(onto, r))

    new = Ontology()
    for tid, term in onto.terms.items():
        if tid in keep:
            parents = [p for p in term.is_a_parents if p in keep]
            new.add_term(
                type(term)(
                    term_id=term.term_id,
                    name=term.name,
                    namespace=term.namespace,
                    definition=term.definition,
                    alt_ids=list(term.alt_ids),
                    is_a_parents=parents,
                )
            )
    return new
