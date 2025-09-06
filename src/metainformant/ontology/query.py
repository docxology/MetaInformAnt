from __future__ import annotations

from collections import deque
from typing import Iterable, List, Set

from .types import Ontology


def ancestors(onto: Ontology, term_id: str) -> Set[str]:
    """Return transitive parent set (excluding self)."""
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
    """Return transitive child set (excluding self)."""
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
    """Extract subgraph induced by all descendants of roots (including roots)."""
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
