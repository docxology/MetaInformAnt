"""Tests for ontology types module.

All tests follow real-implementation policy and use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.ontology.core.types import Ontology, Term


class TestTerm:
    """Test Term dataclass."""

    def test_term_creation_basic(self):
        """Test basic term creation."""
        term = Term(term_id="GO:0008150", name="biological_process", namespace="biological_process")
        assert term.term_id == "GO:0008150"
        assert term.name == "biological_process"
        assert term.namespace == "biological_process"

    def test_term_creation_with_optional_fields(self):
        """Test term creation with optional fields."""
        term = Term(
            term_id="GO:0003674",
            name="molecular_function",
            namespace="molecular_function",
            definition="A molecular function",
            alt_ids=["GO:0003675"],
            is_a_parents=["GO:0008150"],
        )
        assert term.definition == "A molecular function"
        assert term.alt_ids == ["GO:0003675"]
        assert term.is_a_parents == ["GO:0008150"]

    def test_term_creation_defaults(self):
        """Test term creation with default values."""
        term = Term(term_id="GO:0008150", name="biological_process")
        assert term.namespace is None
        assert term.definition is None
        assert term.alt_ids == []
        assert term.is_a_parents == []


class TestOntology:
    """Test Ontology dataclass."""

    def test_ontology_creation_empty(self):
        """Test creating empty ontology."""
        onto = Ontology()
        assert onto.terms == {}
        assert onto.parents_of == {}
        assert onto.children_of == {}

    def test_ontology_creation_with_terms(self):
        """Test creating ontology with terms."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])

        onto = Ontology(terms={"A": term_a, "B": term_b}, parents_of={"B": {"A"}}, children_of={"A": {"B"}})
        assert len(onto.terms) == 2
        assert "A" in onto.terms
        assert "B" in onto.terms

    def test_ontology_parents_of(self):
        """Test parents_of mapping."""
        onto = Ontology(terms={}, parents_of={"B": {"A"}, "C": {"B"}}, children_of={})
        assert "A" in onto.parents_of["B"]
        assert "B" in onto.parents_of["C"]

    def test_ontology_children_of(self):
        """Test children_of mapping."""
        onto = Ontology(terms={}, parents_of={}, children_of={"A": {"B"}, "B": {"C"}})
        assert "B" in onto.children_of["A"]
        assert "C" in onto.children_of["B"]

    def test_ontology_get_term(self):
        """Test getting term by ID via dictionary access."""
        term = Term("GO:0008150", "biological_process", namespace="biological_process")
        onto = Ontology(terms={"GO:0008150": term})
        retrieved = onto.terms.get("GO:0008150")
        assert retrieved == term

    def test_ontology_get_term_missing(self):
        """Test getting non-existent term."""
        onto = Ontology(terms={})
        retrieved = onto.terms.get("NONEXISTENT")
        assert retrieved is None

    def test_ontology_has_term(self):
        """Test checking if term exists."""
        term = Term("GO:0008150", "biological_process")
        onto = Ontology(terms={"GO:0008150": term})
        assert onto.has_term("GO:0008150") is True
        assert onto.has_term("NONEXISTENT") is False

    def test_ontology_get_term_method(self):
        """Test getting term by ID."""
        term = Term("GO:0008150", "biological_process")
        onto = Ontology(terms={"GO:0008150": term})
        retrieved = onto.get_term("GO:0008150")
        assert retrieved == term
        assert onto.get_term("NONEXISTENT") is None

    def test_ontology_get_namespace(self):
        """Test getting namespace for term."""
        term = Term("GO:0008150", "biological_process", namespace="biological_process")
        onto = Ontology(terms={"GO:0008150": term})
        assert onto.get_namespace("GO:0008150") == "biological_process"
        assert onto.get_namespace("NONEXISTENT") is None

    def test_ontology_add_term_validation(self):
        """Test add_term validation."""
        onto = Ontology()
        from metainformant.core.utils.errors import ValidationError

        # Empty term_id
        with pytest.raises(ValidationError, match="Term ID cannot be empty"):
            onto.add_term(Term("", "name"))

        # Empty name
        with pytest.raises(ValidationError, match="Term name cannot be empty"):
            onto.add_term(Term("GO:001", ""))

        # Duplicate term_id
        term = Term("GO:001", "name")
        onto.add_term(term)
        with pytest.raises(ValidationError, match="already exists"):
            onto.add_term(term)

    def test_ontology_validate(self):
        """Test ontology validation."""
        onto = Ontology()

        # Valid ontology
        term = Term("GO:001", "name")
        onto.add_term(term)
        is_valid, errors = onto.validate()
        assert is_valid
        assert len(errors) == 0

        # Orphaned parent reference
        term2 = Term("GO:002", "name2", is_a_parents=["GO:999"])
        onto.add_term(term2)
        is_valid, errors = onto.validate()
        assert not is_valid
        assert any("GO:999" in error for error in errors)

    def test_term_repr(self):
        """Test Term __repr__."""
        term = Term("GO:001", "name", namespace="test")
        repr_str = repr(term)
        assert "GO:001" in repr_str
        assert "name" in repr_str
        assert "test" in repr_str

    def test_term_with_relationships(self):
        """Test Term with relationships."""
        term = Term(
            term_id="GO:001",
            name="name",
            relationships={"part_of": ["GO:002"], "regulates": ["GO:003"]},
            synonyms=["syn1"],
            xrefs=["xref:001"],
            subsets=["subset1"],
        )
        assert "part_of" in term.relationships
        assert term.synonyms == ["syn1"]
        assert term.xrefs == ["xref:001"]
        assert term.subsets == ["subset1"]

    def test_ontology_get_relationships(self):
        """Test get_relationships method."""
        term = Term(
            term_id="GO:001",
            name="name",
            is_a_parents=["GO:002"],
            relationships={"part_of": ["GO:003"]},
        )
        onto = Ontology()
        onto.add_term(term)

        all_rels = onto.get_relationships("GO:001")
        assert "is_a" in all_rels
        assert "part_of" in all_rels
        assert "GO:002" in all_rels["is_a"]
        assert "GO:003" in all_rels["part_of"]

        part_of_rels = onto.get_relationships("GO:001", rel_type="part_of")
        assert isinstance(part_of_rels, set)
        assert "GO:003" in part_of_rels

    def test_ontology_get_relationships_missing_term(self):
        """Test get_relationships with missing term."""
        onto = Ontology()
        with pytest.raises(KeyError, match="not found"):
            onto.get_relationships("NONEXISTENT")


class TestFactories:
    """Test create_term / create_relationship / create_ontology factories."""

    def test_create_term_factory(self):
        from metainformant.ontology.core.types import create_term

        term = create_term(
            id="GO:0008150",
            name="biological_process",
            namespace="biological_process",
            synonyms=["bp"],
            extra_key="extra",
        )
        assert term.id == "GO:0008150"
        assert term.term_id == "GO:0008150"
        assert term.synonyms == ["bp"]
        assert term.metadata == {"extra_key": "extra"}
        assert str(term) == "GO:0008150: biological_process"

    def test_create_relationship_factory(self):
        from metainformant.ontology.core.types import create_relationship

        rel = create_relationship("GO:0009987", "GO:0008150", "is_a", evidence="IEA")
        assert rel.source == "GO:0009987"
        assert rel.metadata == {"evidence": "IEA"}
        assert str(rel) == "GO:0009987 is_a GO:0008150"

    def test_create_relationship_validation(self):
        from metainformant.core.utils.errors import ValidationError
        from metainformant.ontology.core.types import create_relationship

        with pytest.raises(ValidationError, match="source cannot be empty"):
            create_relationship("", "GO:1", "is_a")

    def test_create_ontology_factory(self):
        from metainformant.ontology.core.types import create_ontology, create_relationship, create_term

        t = create_term(id="GO:1", name="child", is_a_parents=["GO:0"])
        root = create_term(id="GO:0", name="root")
        onto = create_ontology(
            terms={"GO:0": root, "GO:1": t},
            relationships=[create_relationship("GO:1", "GO:0", "is_a")],
            format_version="1.2",
        )
        assert len(onto) == 2
        assert onto.metadata["format_version"] == "1.2"
        assert onto.parents_of["GO:1"] == {"GO:0"}


class TestOntologyTraversal:
    """Test ancestor/descendant and non-is_a branches of Ontology methods."""

    def _ontology(self) -> Ontology:
        onto = Ontology()
        onto.add_term(Term("A", "root"))
        onto.add_term(Term("B", "child", is_a_parents=["A"]))
        onto.add_term(Term("C", "grandchild", is_a_parents=["B"]))
        return onto

    def test_get_ancestors_and_descendants(self):
        onto = self._ontology()
        assert onto.get_ancestors("C") == {"A", "B"}
        assert onto.get_descendants("A") == {"B", "C"}

    def test_get_children_and_parents_non_is_a(self):
        from metainformant.ontology.core.types import create_relationship

        onto = self._ontology()
        onto.add_relationship(create_relationship("B", "A", "part_of"))
        assert onto.get_children("A", relation_type="part_of") == {"B"}
        assert onto.get_parents("B", relation_type="part_of") == {"A"}

    def test_add_relationship_updates_is_a_mappings(self):
        from metainformant.ontology.core.types import create_relationship

        onto = self._ontology()
        onto.add_term(Term("D", "other", is_a_parents=[]))
        onto.add_relationship(create_relationship("D", "A", "is_a"))
        assert onto.parents_of["D"] == {"A"}
        assert "D" in onto.children_of["A"]

    def test_str_and_len(self):
        onto = self._ontology()
        assert len(onto) == 3
        assert "3 terms" in str(onto)

    def test_get_roots_and_leaves(self):
        onto = self._ontology()
        assert onto.get_roots() == {"A"}
        assert onto.get_leaves() == {"C"}
