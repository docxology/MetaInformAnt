from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set


@dataclass(slots=True)
class Term:
    """A minimal ontology term.

    Captures essential GO-like fields sufficient for many analyses and tests.
    """

    term_id: str
    name: str
    namespace: str | None = None
    definition: str | None = None
    alt_ids: List[str] = field(default_factory=list)
    is_a_parents: List[str] = field(default_factory=list)


@dataclass(slots=True)
class Ontology:
    """In-memory directed acyclic graph of ontology terms.

    The graph is represented as adjacency maps for fast traversal.
    """

    terms: Dict[str, Term] = field(default_factory=dict)
    parents_of: Dict[str, Set[str]] = field(default_factory=dict)
    children_of: Dict[str, Set[str]] = field(default_factory=dict)

    def add_term(self, term: Term) -> None:
        self.terms[term.term_id] = term
        self.parents_of.setdefault(term.term_id, set()).update(term.is_a_parents)
        for parent_id in term.is_a_parents:
            self.children_of.setdefault(parent_id, set()).add(term.term_id)

    def has_term(self, term_id: str) -> bool:
        return term_id in self.terms

    def num_terms(self) -> int:
        return len(self.terms)
