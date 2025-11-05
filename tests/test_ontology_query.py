"""Tests for ontology query module.

All tests follow NO_MOCKING_POLICY and use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.ontology.query import ancestors, descendants, subgraph
from metainformant.ontology.types import Ontology, Term


class TestAncestors:
    """Test ancestors query function."""

    def test_ancestors_simple_hierarchy(self):
        """Test finding ancestors in simple hierarchy."""
        # Create simple ontology: A -> B -> C
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        term_c = Term("C", "Term C", namespace="test", is_a_parents=["B"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b, "C": term_c},
            parents_of={"B": {"A"}, "C": {"B"}},
            children_of={"A": {"B"}, "B": {"C"}}
        )
        
        # Ancestors of C should include B and A
        anc = ancestors(onto, "C")
        assert "B" in anc
        assert "A" in anc
        assert "C" not in anc  # Excludes self

    def test_ancestors_no_parents(self):
        """Test ancestors of root term."""
        term_a = Term("A", "Term A", namespace="test")
        onto = Ontology(
            terms={"A": term_a},
            parents_of={},
            children_of={}
        )
        
        anc = ancestors(onto, "A")
        assert anc == set()  # No ancestors

    def test_ancestors_missing_term(self):
        """Test ancestors of non-existent term."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        anc = ancestors(onto, "NONEXISTENT")
        assert anc == set()


class TestDescendants:
    """Test descendants query function."""

    def test_descendants_simple_hierarchy(self):
        """Test finding descendants in simple hierarchy."""
        # Create simple ontology: A -> B -> C
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        term_c = Term("C", "Term C", namespace="test", is_a_parents=["B"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b, "C": term_c},
            parents_of={"B": {"A"}, "C": {"B"}},
            children_of={"A": {"B"}, "B": {"C"}}
        )
        
        # Descendants of A should include B and C
        desc = descendants(onto, "A")
        assert "B" in desc
        assert "C" in desc
        assert "A" not in desc  # Excludes self

    def test_descendants_no_children(self):
        """Test descendants of leaf term."""
        term_c = Term("C", "Term C", namespace="test")
        onto = Ontology(
            terms={"C": term_c},
            parents_of={},
            children_of={}
        )
        
        desc = descendants(onto, "C")
        assert desc == set()  # No descendants

    def test_descendants_missing_term(self):
        """Test descendants of non-existent term."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        desc = descendants(onto, "NONEXISTENT")
        assert desc == set()


class TestSubgraph:
    """Test subgraph extraction function."""

    def test_subgraph_single_term(self):
        """Test extracting subgraph with single term."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b},
            parents_of={"B": {"A"}},
            children_of={"A": {"B"}}
        )
        
        sub = subgraph(onto, ["A"])
        assert "A" in sub.terms
        assert len(sub.terms) == 1

    def test_subgraph_multiple_terms(self):
        """Test extracting subgraph with multiple terms."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        term_c = Term("C", "Term C", namespace="test", is_a_parents=["B"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b, "C": term_c},
            parents_of={"B": {"A"}, "C": {"B"}},
            children_of={"A": {"B"}, "B": {"C"}}
        )
        
        sub = subgraph(onto, ["A", "C"])
        assert "A" in sub.terms
        assert "C" in sub.terms
        assert len(sub.terms) == 2

    def test_subgraph_empty_seed(self):
        """Test subgraph with empty seed set."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        sub = subgraph(onto, [])
        assert len(sub.terms) == 0


