from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set


@dataclass(slots=True)
class Term:
    """Represents an ontology term (e.g., Gene Ontology term).
    
    Minimal but sufficient representation for ontology terms with essential
    fields for hierarchy traversal and basic analysis.
    
    Attributes:
        term_id: Unique term identifier (e.g., "GO:0008150")
        name: Human-readable term name
        namespace: Ontology namespace (e.g., "biological_process", "molecular_function")
        definition: Text definition of the term
        alt_ids: List of alternative identifiers for this term
        is_a_parents: List of parent term IDs (is_a relationships)
        
    Examples:
        >>> term = Term(
        ...     term_id="GO:0008150",
        ...     name="biological_process",
        ...     namespace="biological_process"
        ... )
        >>> term.term_id
        'GO:0008150'
    """

    term_id: str
    name: str
    namespace: str | None = None
    definition: str | None = None
    alt_ids: List[str] = field(default_factory=list)
    is_a_parents: List[str] = field(default_factory=list)


@dataclass(slots=True)
class Ontology:
    """In-memory representation of an ontology as a directed acyclic graph.
    
    Efficient structure for ontology terms and their hierarchical relationships.
    Uses adjacency maps for fast parent/child lookups and traversal.
    
    Attributes:
        terms: Dictionary mapping term_id -> Term objects
        parents_of: Dictionary mapping term_id -> set of parent term IDs
        children_of: Dictionary mapping term_id -> set of child term IDs
        
    Examples:
        >>> onto = Ontology()
        >>> term = Term(term_id="GO:001", name="test", is_a_parents=["GO:002"])
        >>> onto.add_term(term)
        >>> onto.has_term("GO:001")
        True
        >>> onto.num_terms()
        1
    """

    terms: Dict[str, Term] = field(default_factory=dict)
    parents_of: Dict[str, Set[str]] = field(default_factory=dict)
    children_of: Dict[str, Set[str]] = field(default_factory=dict)

    def add_term(self, term: Term) -> None:
        """Add a term to the ontology and update parent/child relationships.
        
        Args:
            term: Term object to add
        """
        self.terms[term.term_id] = term
        self.parents_of.setdefault(term.term_id, set()).update(term.is_a_parents)
        for parent_id in term.is_a_parents:
            self.children_of.setdefault(parent_id, set()).add(term.term_id)

    def has_term(self, term_id: str) -> bool:
        """Check if term exists in ontology.
        
        Args:
            term_id: Term identifier to check
            
        Returns:
            True if term exists, False otherwise
        """
        return term_id in self.terms

    def num_terms(self) -> int:
        """Get total number of terms in ontology.
        
        Returns:
            Count of terms
        """
        return len(self.terms)
