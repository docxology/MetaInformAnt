from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set

from metainformant.core.errors import ValidationError


@dataclass(slots=True)
class Term:
    """Represents an ontology term (e.g., Gene Ontology term).
    
    Comprehensive representation for ontology terms with essential fields
    for hierarchy traversal and analysis. Supports multiple relationship types,
    synonyms, cross-references, and subsets.
    
    Attributes:
        term_id: Unique term identifier (e.g., "GO:0008150")
        name: Human-readable term name
        namespace: Ontology namespace (e.g., "biological_process", "molecular_function")
        definition: Text definition of the term
        alt_ids: List of alternative identifiers for this term
        is_a_parents: List of parent term IDs (is_a relationships)
        relationships: Dictionary mapping relationship type -> list of term IDs
            (e.g., {"part_of": ["GO:001"], "regulates": ["GO:002"]})
        synonyms: List of synonym strings (alternative names)
        xrefs: List of cross-reference strings (database references)
        subsets: List of subset names (e.g., GO slims)
        
    Examples:
        >>> term = Term(
        ...     term_id="GO:0008150",
        ...     name="biological_process",
        ...     namespace="biological_process"
        ... )
        >>> term.term_id
        'GO:0008150'
        >>> repr(term)
        "Term(term_id='GO:0008150', name='biological_process', ...)"
    """

    term_id: str
    name: str
    namespace: str | None = None
    definition: str | None = None
    alt_ids: List[str] = field(default_factory=list)
    is_a_parents: List[str] = field(default_factory=list)
    relationships: Dict[str, List[str]] = field(default_factory=dict)
    synonyms: List[str] = field(default_factory=list)
    xrefs: List[str] = field(default_factory=list)
    subsets: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        """Return string representation of Term."""
        parts = [f"term_id='{self.term_id}'", f"name='{self.name}'"]
        if self.namespace:
            parts.append(f"namespace='{self.namespace}'")
        if self.definition:
            parts.append(f"definition='{self.definition[:50]}...'" if len(self.definition) > 50 else f"definition='{self.definition}'")
        if self.alt_ids:
            parts.append(f"alt_ids={len(self.alt_ids)} alt_ids")
        if self.is_a_parents:
            parts.append(f"is_a_parents={len(self.is_a_parents)} parents")
        if self.relationships:
            parts.append(f"relationships={len(self.relationships)} types")
        if self.synonyms:
            parts.append(f"synonyms={len(self.synonyms)}")
        if self.xrefs:
            parts.append(f"xrefs={len(self.xrefs)}")
        if self.subsets:
            parts.append(f"subsets={len(self.subsets)}")
        return f"Term({', '.join(parts)})"


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
        
        Updates both is_a relationships and other relationship types stored in
        term.relationships. Maintains parent/child mappings for is_a relationships.
        
        Args:
            term: Term object to add
            
        Raises:
            ValidationError: If term_id is empty, name is empty, or term_id already exists
        """
        if not term.term_id:
            raise ValidationError("Term ID cannot be empty")
        if not term.name:
            raise ValidationError(f"Term name cannot be empty for term {term.term_id}")
        if term.term_id in self.terms:
            raise ValidationError(f"Term {term.term_id} already exists in ontology")
        
        # Validate parent references exist (warn but don't fail - parents may be added later)
        for parent_id in term.is_a_parents:
            if parent_id not in self.terms and parent_id not in self.parents_of:
                # This is acceptable - parent may be added later or may be external reference
                pass
        
        self.terms[term.term_id] = term
        self.parents_of.setdefault(term.term_id, set()).update(term.is_a_parents)
        for parent_id in term.is_a_parents:
            self.children_of.setdefault(parent_id, set()).add(term.term_id)
    
    def get_relationships(self, term_id: str, rel_type: str | None = None) -> Dict[str, Set[str]] | Set[str]:
        """Get relationships for a term.
        
        Args:
            term_id: Term identifier
            rel_type: Optional relationship type filter (e.g., "part_of", "regulates")
                If None, returns all relationships as a dict.
                If specified, returns only that relationship type as a set.
        
        Returns:
            If rel_type is None: Dictionary mapping relationship type -> set of term IDs
            If rel_type is specified: Set of term IDs for that relationship type
            
        Raises:
            KeyError: If term_id not found
        """
        if term_id not in self.terms:
            raise KeyError(f"Term {term_id} not found in ontology")
        
        term = self.terms[term_id]
        all_rels: Dict[str, Set[str]] = {}
        
        # Include is_a relationships
        if term.is_a_parents:
            all_rels["is_a"] = set(term.is_a_parents)
        
        # Include other relationship types
        for rel_type_name, rel_terms in term.relationships.items():
            all_rels[rel_type_name] = set(rel_terms)
        
        if rel_type is None:
            return all_rels
        else:
            return all_rels.get(rel_type, set())

    def get_term(self, term_id: str) -> Term | None:
        """Get term by ID.
        
        Args:
            term_id: Term identifier to retrieve
            
        Returns:
            Term object if found, None otherwise
        """
        return self.terms.get(term_id)

    def get_namespace(self, term_id: str) -> str | None:
        """Get namespace for a term.
        
        Args:
            term_id: Term identifier
            
        Returns:
            Namespace string if term exists and has namespace, None otherwise
        """
        term = self.terms.get(term_id)
        return term.namespace if term else None

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

    def validate(self) -> tuple[bool, List[str]]:
        """Validate ontology integrity (check for cycles, orphaned terms).
        
        Checks for:
        - Cycles in parent-child relationships
        - Orphaned terms (terms referenced as parents but not in ontology)
        - Terms with invalid parent references
        
        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        errors: List[str] = []
        
        # Check for cycles using DFS
        visited: Set[str] = set()
        rec_stack: Set[str] = set()
        
        def has_cycle(node: str) -> bool:
            """Check for cycles starting from node."""
            visited.add(node)
            rec_stack.add(node)
            
            for child in self.children_of.get(node, set()):
                if child not in self.terms:
                    errors.append(f"Child term {child} referenced by {node} but not in ontology")
                elif child not in visited:
                    if has_cycle(child):
                        return True
                elif child in rec_stack:
                    errors.append(f"Cycle detected: {node} -> {child} (and ancestors)")
                    return True
            
            rec_stack.remove(node)
            return False
        
        # Check all terms for cycles
        for term_id in self.terms:
            if term_id not in visited:
                has_cycle(term_id)
        
        # Check for orphaned parent references
        for term_id, term in self.terms.items():
            for parent_id in term.is_a_parents:
                if parent_id not in self.terms:
                    errors.append(f"Term {term_id} references parent {parent_id} which is not in ontology")
        
        return len(errors) == 0, errors
