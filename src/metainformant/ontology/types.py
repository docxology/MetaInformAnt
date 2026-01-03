"""Ontology data types and structures.

This module defines the core data structures for representing ontologies,
terms, and relationships in METAINFORMANT.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Any
from metainformant.core import logging, validation

logger = logging.get_logger(__name__)


@dataclass
class Term:
    """Represents an ontology term.

    An ontology term contains metadata about a biological concept,
    such as a gene, process, or molecular function.
    """
    id: str
    name: Optional[str] = None
    definition: Optional[str] = None
    namespace: Optional[str] = None
    synonyms: List[str] = field(default_factory=list)
    xrefs: List[str] = field(default_factory=list)
    is_obsolete: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate term after initialization."""
        validation.validate_not_none(self.id, "id")
        validation.validate_not_empty(self.id, "id")

    def add_synonym(self, synonym: str) -> None:
        """Add a synonym to this term."""
        if synonym not in self.synonyms:
            self.synonyms.append(synonym)

    def add_xref(self, xref: str) -> None:
        """Add a cross-reference to this term."""
        if xref not in self.xrefs:
            self.xrefs.append(xref)

    def __str__(self) -> str:
        """String representation of the term."""
        return f"{self.id}: {self.name or 'Unnamed'}"


@dataclass
class Relationship:
    """Represents a relationship between ontology terms.

    Relationships define how terms are connected in an ontology hierarchy,
    such as "is_a", "part_of", or "regulates".
    """
    source: str
    target: str
    relation_type: str
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate relationship after initialization."""
        validation.validate_not_none(self.source, "source")
        validation.validate_not_none(self.target, "target")
        validation.validate_not_none(self.relation_type, "relation_type")
        validation.validate_not_empty(self.source, "source")
        validation.validate_not_empty(self.target, "target")
        validation.validate_not_empty(self.relation_type, "relation_type")

    def __str__(self) -> str:
        """String representation of the relationship."""
        return f"{self.source} {self.relation_type} {self.target}"


@dataclass
class Ontology:
    """Represents a complete ontology.

    An ontology contains terms and their relationships, along with
    metadata about the ontology itself.
    """
    terms: Dict[str, Term] = field(default_factory=dict)
    relationships: List[Relationship] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def add_term(self, term: Term) -> None:
        """Add a term to the ontology."""
        self.terms[term.id] = term
        logger.debug(f"Added term {term.id} to ontology")

    def add_relationship(self, relationship: Relationship) -> None:
        """Add a relationship to the ontology."""
        self.relationships.append(relationship)
        logger.debug(f"Added relationship {relationship}")

    def get_term(self, term_id: str) -> Optional[Term]:
        """Get a term by ID."""
        return self.terms.get(term_id)

    def get_children(self, term_id: str, relation_type: str = "is_a") -> Set[str]:
        """Get all child terms of a given term."""
        children = set()
        for rel in self.relationships:
            if rel.target == term_id and rel.relation_type == relation_type:
                children.add(rel.source)
        return children

    def get_parents(self, term_id: str, relation_type: str = "is_a") -> Set[str]:
        """Get all parent terms of a given term."""
        parents = set()
        for rel in self.relationships:
            if rel.source == term_id and rel.relation_type == relation_type:
                parents.add(rel.target)
        return parents

    def get_ancestors(self, term_id: str, relation_type: str = "is_a") -> Set[str]:
        """Get all ancestor terms (recursive parents)."""
        ancestors = set()
        to_visit = list(self.get_parents(term_id, relation_type))

        while to_visit:
            current = to_visit.pop()
            if current not in ancestors:
                ancestors.add(current)
                to_visit.extend(self.get_parents(current, relation_type))

        return ancestors

    def get_descendants(self, term_id: str, relation_type: str = "is_a") -> Set[str]:
        """Get all descendant terms (recursive children)."""
        descendants = set()
        to_visit = list(self.get_children(term_id, relation_type))

        while to_visit:
            current = to_visit.pop()
            if current not in descendants:
                descendants.add(current)
                to_visit.extend(self.get_children(current, relation_type))

        return descendants

    def get_roots(self, relation_type: str = "is_a") -> Set[str]:
        """Get all root terms (terms with no parents)."""
        all_terms = set(self.terms.keys())
        terms_with_parents = {rel.source for rel in self.relationships
                            if rel.relation_type == relation_type}
        return all_terms - terms_with_parents

    def get_leaves(self, relation_type: str = "is_a") -> Set[str]:
        """Get all leaf terms (terms with no children)."""
        all_terms = set(self.terms.keys())
        terms_with_children = {rel.target for rel in self.relationships
                             if rel.relation_type == relation_type}
        return all_terms - terms_with_children

    def __len__(self) -> int:
        """Return the number of terms in the ontology."""
        return len(self.terms)

    def __str__(self) -> str:
        """String representation of the ontology."""
        return f"Ontology with {len(self.terms)} terms and {len(self.relationships)} relationships"


def create_term(
    id: str,
    name: Optional[str] = None,
    definition: Optional[str] = None,
    namespace: Optional[str] = None,
    synonyms: Optional[List[str]] = None,
    xrefs: Optional[List[str]] = None,
    is_obsolete: bool = False,
    **metadata: Any
) -> Term:
    """Create a new ontology term.

    Args:
        id: Unique identifier for the term.
        name: Human-readable name for the term.
        definition: Definition or description of the term.
        namespace: Ontology namespace (e.g., "biological_process").
        synonyms: Alternative names for the term.
        xrefs: Cross-references to other databases.
        is_obsolete: Whether this term is obsolete.
        **metadata: Additional metadata for the term.

    Returns:
        A new Term instance.

    Examples:
        >>> term = create_term(
        ...     id="GO:0008150",
        ...     name="biological_process",
        ...     definition="Any process that is carried out by living organisms.",
        ...     namespace="biological_process"
        ... )
        >>> str(term)
        'GO:0008150: biological_process'
    """
    return Term(
        id=id,
        name=name,
        definition=definition,
        namespace=namespace,
        synonyms=synonyms or [],
        xrefs=xrefs or [],
        is_obsolete=is_obsolete,
        metadata=metadata
    )


def create_relationship(
    source: str,
    target: str,
    relation_type: str,
    **metadata: Any
) -> Relationship:
    """Create a new relationship between ontology terms.

    Args:
        source: ID of the source term.
        target: ID of the target term.
        relation_type: Type of relationship (e.g., "is_a", "part_of").
        **metadata: Additional metadata for the relationship.

    Returns:
        A new Relationship instance.

    Examples:
        >>> rel = create_relationship(
        ...     source="GO:0009987",
        ...     target="GO:0008150",
        ...     relation_type="is_a"
        ... )
        >>> str(rel)
        'GO:0009987 is_a GO:0008150'
    """
    return Relationship(
        source=source,
        target=target,
        relation_type=relation_type,
        metadata=metadata
    )


def create_ontology(
    terms: Optional[Dict[str, Term]] = None,
    relationships: Optional[List[Relationship]] = None,
    **metadata: Any
) -> Ontology:
    """Create a new ontology.

    Args:
        terms: Dictionary of term ID to Term objects.
        relationships: List of Relationship objects.
        **metadata: Additional metadata for the ontology.

    Returns:
        A new Ontology instance.

    Examples:
        >>> ontology = create_ontology(
        ...     terms={"GO:0008150": create_term("GO:0008150", name="biological_process")},
        ...     metadata={"format-version": "1.2"}
        ... )
        >>> len(ontology)
        1
    """
    return Ontology(
        terms=terms or {},
        relationships=relationships or [],
        metadata=metadata
    )





