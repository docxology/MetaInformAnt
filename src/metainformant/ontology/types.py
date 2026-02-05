"""Ontology data types and structures.

This module defines the core data structures for representing ontologies,
terms, and relationships in METAINFORMANT.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Any, Tuple
from metainformant.core import logging
from metainformant.core.utils import errors

logger = logging.get_logger(__name__)


@dataclass
class Term:
    """Represents an ontology term.

    An ontology term contains metadata about a biological concept,
    such as a gene, process, or molecular function.
    """

    term_id: str = ""
    name: Optional[str] = None
    namespace: Optional[str] = None
    definition: Optional[str] = None
    synonyms: List[str] = field(default_factory=list)
    xrefs: List[str] = field(default_factory=list)
    alt_ids: List[str] = field(default_factory=list)
    is_a_parents: List[str] = field(default_factory=list)
    relationships: Dict[str, List[str]] = field(default_factory=dict)
    subsets: List[str] = field(default_factory=list)
    is_obsolete: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)
    # Backward compatibility alias
    id: str = ""  # Alias for term_id

    def __post_init__(self) -> None:
        """Validate term after initialization."""
        # Handle id/term_id aliases - prefer term_id
        if self.term_id and not self.id:
            self.id = self.term_id
        elif self.id and not self.term_id:
            self.term_id = self.id

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
        return f"{self.term_id}: {self.name or 'Unnamed'}"

    def __repr__(self) -> str:
        """Detailed representation of the term."""
        return f"Term(term_id={self.term_id!r}, name={self.name!r}, namespace={self.namespace!r})"


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
        if not self.source:
            raise errors.ValidationError("source cannot be empty")
        if not self.target:
            raise errors.ValidationError("target cannot be empty")
        if not self.relation_type:
            raise errors.ValidationError("relation_type cannot be empty")

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
    parents_of: Dict[str, Set[str]] = field(default_factory=dict)
    children_of: Dict[str, Set[str]] = field(default_factory=dict)
    relationships: List[Relationship] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def add_term(self, term: Term) -> None:
        """Add a term to the ontology with validation."""
        # Validate term
        if not term.term_id:
            raise errors.ValidationError("Term ID cannot be empty")
        if not term.name:
            raise errors.ValidationError("Term name cannot be empty")
        if term.term_id in self.terms:
            raise errors.ValidationError(f"Term '{term.term_id}' already exists in ontology")

        self.terms[term.term_id] = term

        # Update parent/child mappings from is_a_parents
        for parent_id in term.is_a_parents:
            if term.term_id not in self.parents_of:
                self.parents_of[term.term_id] = set()
            self.parents_of[term.term_id].add(parent_id)

            if parent_id not in self.children_of:
                self.children_of[parent_id] = set()
            self.children_of[parent_id].add(term.term_id)

        logger.debug(f"Added term {term.term_id} to ontology")

    def has_term(self, term_id: str) -> bool:
        """Check if a term exists in the ontology."""
        return term_id in self.terms

    def get_term(self, term_id: str) -> Optional[Term]:
        """Get a term by ID."""
        return self.terms.get(term_id)

    def get_namespace(self, term_id: str) -> Optional[str]:
        """Get the namespace of a term."""
        term = self.terms.get(term_id)
        return term.namespace if term else None

    def get_relationships(self, term_id: str, rel_type: Optional[str] = None) -> Dict[str, Set[str]] | Set[str]:
        """Get relationships for a term.

        Args:
            term_id: The term ID to get relationships for
            rel_type: Optional relationship type to filter by

        Returns:
            If rel_type is None, returns dict of all relationship types to targets
            If rel_type is specified, returns set of targets for that relationship type

        Raises:
            KeyError: If term not found
        """
        if term_id not in self.terms:
            raise KeyError(f"Term '{term_id}' not found in ontology")

        term = self.terms[term_id]
        all_rels: Dict[str, Set[str]] = {}

        # Add is_a relationships from is_a_parents
        if term.is_a_parents:
            all_rels["is_a"] = set(term.is_a_parents)

        # Add other relationships
        for rel_name, targets in term.relationships.items():
            all_rels[rel_name] = set(targets)

        if rel_type is not None:
            return all_rels.get(rel_type, set())

        return all_rels

    def validate(self) -> Tuple[bool, List[str]]:
        """Validate the ontology structure.

        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        validation_errors: List[str] = []

        for term_id, term in self.terms.items():
            # Check that all parents exist
            for parent_id in term.is_a_parents:
                if parent_id not in self.terms:
                    validation_errors.append(f"Term '{term_id}' references non-existent parent '{parent_id}'")

            # Check relationships
            for rel_type, targets in term.relationships.items():
                for target_id in targets:
                    if target_id not in self.terms:
                        validation_errors.append(
                            f"Term '{term_id}' has {rel_type} relationship to non-existent term '{target_id}'"
                        )

        is_valid = len(validation_errors) == 0
        return is_valid, validation_errors

    def add_relationship(self, relationship: Relationship) -> None:
        """Add a relationship to the ontology."""
        self.relationships.append(relationship)

        # Update parent/child mappings if it's an is_a relationship
        if relationship.relation_type == "is_a":
            if relationship.source not in self.parents_of:
                self.parents_of[relationship.source] = set()
            self.parents_of[relationship.source].add(relationship.target)

            if relationship.target not in self.children_of:
                self.children_of[relationship.target] = set()
            self.children_of[relationship.target].add(relationship.source)

        logger.debug(f"Added relationship {relationship}")

    def get_children(self, term_id: str, relation_type: str = "is_a") -> Set[str]:
        """Get all child terms of a given term."""
        if relation_type == "is_a":
            return self.children_of.get(term_id, set())
        children = set()
        for rel in self.relationships:
            if rel.target == term_id and rel.relation_type == relation_type:
                children.add(rel.source)
        return children

    def get_parents(self, term_id: str, relation_type: str = "is_a") -> Set[str]:
        """Get all parent terms of a given term."""
        if relation_type == "is_a":
            return self.parents_of.get(term_id, set())
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
        terms_with_parents = set(self.parents_of.keys())
        return all_terms - terms_with_parents

    def get_leaves(self, relation_type: str = "is_a") -> Set[str]:
        """Get all leaf terms (terms with no children)."""
        all_terms = set(self.terms.keys())
        terms_with_children = set(self.children_of.keys())
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
    **metadata: Any,
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
        term_id=id,
        name=name,
        definition=definition,
        namespace=namespace,
        synonyms=synonyms or [],
        xrefs=xrefs or [],
        is_obsolete=is_obsolete,
        metadata=metadata,
    )


def create_relationship(source: str, target: str, relation_type: str, **metadata: Any) -> Relationship:
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
    return Relationship(source=source, target=target, relation_type=relation_type, metadata=metadata)


def create_ontology(
    terms: Optional[Dict[str, Term]] = None, relationships: Optional[List[Relationship]] = None, **metadata: Any
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
    return Ontology(terms=terms or {}, relationships=relationships or [], metadata=metadata)
