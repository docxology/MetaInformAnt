"""Tests for ontology types module.

All tests follow NO_MOCKING_POLICY and use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.ontology.types import Ontology, Term


class TestTerm:
    """Test Term dataclass."""

    def test_term_creation_basic(self):
        """Test basic term creation."""
        term = Term(
            term_id="GO:0008150",
            name="biological_process",
            namespace="biological_process"
        )
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
            is_a_parents=["GO:0008150"]
        )
        assert term.definition == "A molecular function"
        assert term.alt_ids == ["GO:0003675"]
        assert term.is_a_parents == ["GO:0008150"]

    def test_term_creation_defaults(self):
        """Test term creation with default values."""
        term = Term(
            term_id="GO:0008150",
            name="biological_process"
        )
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
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b},
            parents_of={"B": {"A"}},
            children_of={"A": {"B"}}
        )
        assert len(onto.terms) == 2
        assert "A" in onto.terms
        assert "B" in onto.terms

    def test_ontology_parents_of(self):
        """Test parents_of mapping."""
        onto = Ontology(
            terms={},
            parents_of={"B": {"A"}, "C": {"B"}},
            children_of={}
        )
        assert "A" in onto.parents_of["B"]
        assert "B" in onto.parents_of["C"]

    def test_ontology_children_of(self):
        """Test children_of mapping."""
        onto = Ontology(
            terms={},
            parents_of={},
            children_of={"A": {"B"}, "B": {"C"}}
        )
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

