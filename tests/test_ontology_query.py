"""Tests for ontology query module.

All tests follow NO_MOCKING_POLICY and use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.ontology.query import (
    ancestors, descendants, subgraph,
    common_ancestors, path_to_root, distance,
    find_term_by_name, filter_by_namespace,
    get_roots, get_leaves, clear_cache, set_cache_enabled
)
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
        """Test ancestors of non-existent term raises ValueError."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        with pytest.raises(ValueError, match="not found in ontology"):
            ancestors(onto, "NONEXISTENT")
    
    def test_ancestors_empty_term_id(self):
        """Test ancestors with empty term_id raises ValueError."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        with pytest.raises(ValueError, match="cannot be empty"):
            ancestors(onto, "")


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
        """Test descendants of non-existent term raises ValueError."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        with pytest.raises(ValueError, match="not found in ontology"):
            descendants(onto, "NONEXISTENT")
    
    def test_descendants_empty_term_id(self):
        """Test descendants with empty term_id raises ValueError."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        with pytest.raises(ValueError, match="cannot be empty"):
            descendants(onto, "")


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
    
    def test_subgraph_invalid_root(self):
        """Test subgraph with invalid root raises ValueError."""
        onto = Ontology(terms={}, parents_of={}, children_of={})
        with pytest.raises(ValueError, match="not found in ontology"):
            subgraph(onto, ["INVALID"])


class TestCommonAncestors:
    """Test common_ancestors function."""
    
    def test_common_ancestors_simple(self):
        """Test finding common ancestors."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        term_c = Term("C", "Term C", namespace="test", is_a_parents=["A"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b, "C": term_c},
            parents_of={"B": {"A"}, "C": {"A"}},
            children_of={"A": {"B", "C"}}
        )
        
        common = common_ancestors(onto, "B", "C")
        assert "A" in common
        assert "B" in common  # Terms themselves are included
        assert "C" in common


class TestPathToRoot:
    """Test path_to_root function."""
    
    def test_path_to_root_simple(self):
        """Test finding path to root."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b},
            parents_of={"B": {"A"}},
            children_of={"A": {"B"}}
        )
        
        path = path_to_root(onto, "B")
        assert path == ["B", "A"]
        assert path[0] == "B"
    
    def test_path_to_root_no_parents(self):
        """Test path to root for term with no parents."""
        term_a = Term("A", "Term A", namespace="test")
        onto = Ontology(terms={"A": term_a}, parents_of={}, children_of={})
        
        path = path_to_root(onto, "A")
        assert path == ["A"]


class TestDistance:
    """Test distance function."""
    
    def test_distance_simple(self):
        """Test calculating distance between terms."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        term_c = Term("C", "Term C", namespace="test", is_a_parents=["B"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b, "C": term_c},
            parents_of={"B": {"A"}, "C": {"B"}},
            children_of={"A": {"B"}, "B": {"C"}}
        )
        
        dist = distance(onto, "C", "A")
        assert dist == 2  # C -> B -> A
    
    def test_distance_same_term(self):
        """Test distance to self is zero."""
        term_a = Term("A", "Term A", namespace="test")
        onto = Ontology(terms={"A": term_a}, parents_of={}, children_of={})
        
        dist = distance(onto, "A", "A")
        assert dist == 0


class TestFindTermByName:
    """Test find_term_by_name function."""
    
    def test_find_term_by_name_exact(self):
        """Test finding term by exact name."""
        term = Term("GO:001", "biological process", namespace="test")
        onto = Ontology(terms={"GO:001": term}, parents_of={}, children_of={})
        
        matches = find_term_by_name(onto, "biological process")
        assert "GO:001" in matches
    
    def test_find_term_by_name_partial(self):
        """Test finding term by partial name."""
        term = Term("GO:001", "biological process", namespace="test")
        onto = Ontology(terms={"GO:001": term}, parents_of={}, children_of={})
        
        matches = find_term_by_name(onto, "biological")
        assert "GO:001" in matches
    
    def test_find_term_by_name_case_insensitive(self):
        """Test case-insensitive matching."""
        term = Term("GO:001", "Biological Process", namespace="test")
        onto = Ontology(terms={"GO:001": term}, parents_of={}, children_of={})
        
        matches = find_term_by_name(onto, "biological process")
        assert "GO:001" in matches
    
    def test_find_term_by_name_with_namespace(self):
        """Test finding term with namespace filter."""
        term1 = Term("GO:001", "process", namespace="biological_process")
        term2 = Term("GO:002", "process", namespace="molecular_function")
        onto = Ontology(
            terms={"GO:001": term1, "GO:002": term2},
            parents_of={}, children_of={}
        )
        
        matches = find_term_by_name(onto, "process", namespace="biological_process")
        assert "GO:001" in matches
        assert "GO:002" not in matches


class TestFilterByNamespace:
    """Test filter_by_namespace function."""
    
    def test_filter_by_namespace(self):
        """Test filtering ontology by namespace."""
        term1 = Term("GO:001", "term1", namespace="biological_process")
        term2 = Term("GO:002", "term2", namespace="molecular_function")
        onto = Ontology(
            terms={"GO:001": term1, "GO:002": term2},
            parents_of={}, children_of={}
        )
        
        filtered = filter_by_namespace(onto, "biological_process")
        assert "GO:001" in filtered.terms
        assert "GO:002" not in filtered.terms
        assert filtered.num_terms() == 1


class TestGetRootsAndLeaves:
    """Test get_roots and get_leaves functions."""
    
    def test_get_roots(self):
        """Test getting root terms."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b},
            parents_of={"B": {"A"}},
            children_of={"A": {"B"}}
        )
        
        roots = get_roots(onto)
        assert "A" in roots
        assert "B" not in roots
    
    def test_get_leaves(self):
        """Test getting leaf terms."""
        term_a = Term("A", "Term A", namespace="test")
        term_b = Term("B", "Term B", namespace="test", is_a_parents=["A"])
        
        onto = Ontology(
            terms={"A": term_a, "B": term_b},
            parents_of={"B": {"A"}},
            children_of={"A": {"B"}}
        )
        
        leaves = get_leaves(onto)
        assert "B" in leaves
        assert "A" not in leaves


class TestCaching:
    """Test caching functionality."""
    
    def test_cache_enabled_disabled(self):
        """Test enabling/disabling cache."""
        set_cache_enabled(True)
        set_cache_enabled(False)
        clear_cache()
    
    def test_clear_cache(self):
        """Test clearing cache."""
        clear_cache()  # Should not raise


