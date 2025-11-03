"""Comprehensive tests for ontology module.

Tests cover OBO parsing, ontology traversal, term queries, and hierarchy operations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.ontology.go import load_go_obo, write_go_summary
from metainformant.ontology.obo import parse_obo
from metainformant.ontology.query import ancestors, descendants, subgraph
from metainformant.ontology.types import Ontology, Term


class TestOntologyTypes:
    """Tests for Ontology and Term classes."""

    def test_term_creation(self):
        """Test creating ontology term."""
        term = Term(
            term_id="GO:0008150",
            name="biological_process",
            namespace="biological_process",
            definition="Any process accomplished by biological systems",
            is_a_parents=["GO:0003674"],
        )

        assert term.term_id == "GO:0008150"
        assert term.name == "biological_process"
        assert len(term.is_a_parents) == 1

    def test_ontology_creation(self):
        """Test creating ontology."""
        onto = Ontology()
        assert onto.num_terms() == 0

        term = Term(term_id="GO:001", name="test_term", is_a_parents=[])
        onto.add_term(term)

        assert onto.num_terms() == 1
        assert onto.has_term("GO:001")

    def test_ontology_hierarchy(self):
        """Test ontology parent-child relationships."""
        onto = Ontology()

        # Parent term
        parent = Term(term_id="GO:001", name="parent", is_a_parents=[])
        onto.add_term(parent)

        # Child term
        child = Term(term_id="GO:002", name="child", is_a_parents=["GO:001"])
        onto.add_term(child)

        assert "GO:001" in onto.children_of
        assert "GO:002" in onto.children_of["GO:001"]
        assert "GO:001" in onto.parents_of["GO:002"]


class TestOntologyQueries:
    """Tests for ontology query functions."""

    def test_ancestors_query(self):
        """Test ancestor retrieval."""
        onto = Ontology()

        # Create hierarchy: root -> intermediate -> leaf
        root = Term(term_id="GO:001", name="root", is_a_parents=[])
        intermediate = Term(term_id="GO:002", name="intermediate", is_a_parents=["GO:001"])
        leaf = Term(term_id="GO:003", name="leaf", is_a_parents=["GO:002"])

        onto.add_term(root)
        onto.add_term(intermediate)
        onto.add_term(leaf)

        # Get ancestors of leaf
        leaf_ancestors = ancestors(onto, "GO:003")

        assert "GO:001" in leaf_ancestors
        assert "GO:002" in leaf_ancestors
        assert "GO:003" not in leaf_ancestors  # Excludes self

    def test_descendants_query(self):
        """Test descendant retrieval."""
        onto = Ontology()

        root = Term(term_id="GO:001", name="root", is_a_parents=[])
        child1 = Term(term_id="GO:002", name="child1", is_a_parents=["GO:001"])
        child2 = Term(term_id="GO:003", name="child2", is_a_parents=["GO:001"])

        onto.add_term(root)
        onto.add_term(child1)
        onto.add_term(child2)

        root_descendants = descendants(onto, "GO:001")

        assert "GO:002" in root_descendants
        assert "GO:003" in root_descendants
        assert "GO:001" not in root_descendants  # Excludes self

    def test_subgraph_extraction(self):
        """Test subgraph extraction."""
        onto = Ontology()

        # Create small hierarchy
        root = Term(term_id="GO:001", name="root", is_a_parents=[])
        child1 = Term(term_id="GO:002", name="child1", is_a_parents=["GO:001"])
        child2 = Term(term_id="GO:003", name="child2", is_a_parents=["GO:001"])
        unrelated = Term(term_id="GO:004", name="unrelated", is_a_parents=[])

        onto.add_term(root)
        onto.add_term(child1)
        onto.add_term(child2)
        onto.add_term(unrelated)

        # Extract subgraph rooted at GO:001
        sub_onto = subgraph(onto, ["GO:001"])

        assert sub_onto.has_term("GO:001")
        assert sub_onto.has_term("GO:002")
        assert sub_onto.has_term("GO:003")
        assert not sub_onto.has_term("GO:004")  # Unrelated term excluded


class TestOBOParsing:
    """Tests for OBO file parsing."""

    def test_parse_simple_obo(self, tmp_path: Path):
        """Test parsing simple OBO file."""
        obo_content = """[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
def: "Any process accomplished by biological systems." [GOC:go_curators]
is_a: GO:0003674

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process
"""

        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        onto = parse_obo(obo_file)

        assert onto.num_terms() == 2
        assert onto.has_term("GO:0008150")
        assert onto.has_term("GO:0009987")

        # Check parent relationship
        assert "GO:0003674" in onto.parents_of["GO:0008150"]
        assert "GO:0008150" in onto.parents_of["GO:0009987"]

    def test_parse_obo_with_alt_ids(self, tmp_path: Path):
        """Test parsing OBO with alternative IDs."""
        obo_content = """[Term]
id: GO:0008150
name: biological_process
alt_id: GO:0000004
is_a: GO:0003674
"""

        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        onto = parse_obo(obo_file)
        term = onto.terms["GO:0008150"]

        assert "GO:0000004" in term.alt_ids


class TestGOFunctions:
    """Tests for GO-specific functions."""

    def test_load_go_obo(self, tmp_path: Path):
        """Test loading GO OBO file."""
        obo_content = """[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
is_a: GO:0003674
"""

        obo_file = tmp_path / "go.obo"
        obo_file.write_text(obo_content)

        onto = load_go_obo(obo_file)

        assert onto.num_terms() == 1
        assert onto.has_term("GO:0008150")

    def test_write_go_summary(self, tmp_path: Path):
        """Test writing GO summary."""
        onto = Ontology()
        term = Term(term_id="GO:001", name="test", is_a_parents=[])
        onto.add_term(term)

        summary_path = write_go_summary(onto, dest=tmp_path / "summary.json")

        assert summary_path.exists()
        import json

        summary = json.loads(summary_path.read_text())
        assert summary["num_terms"] == 1


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_ancestors_missing_term(self):
        """Test ancestors for non-existent term."""
        onto = Ontology()
        ancestors_set = ancestors(onto, "GO:9999999")

        assert ancestors_set == set()

    def test_descendants_missing_term(self):
        """Test descendants for non-existent term."""
        onto = Ontology()
        descendants_set = descendants(onto, "GO:9999999")

        assert descendants_set == set()

    def test_empty_ontology(self):
        """Test operations on empty ontology."""
        onto = Ontology()

        assert onto.num_terms() == 0
        assert not onto.has_term("GO:0008150")

        ancestors_set = ancestors(onto, "GO:0008150")
        assert ancestors_set == set()



