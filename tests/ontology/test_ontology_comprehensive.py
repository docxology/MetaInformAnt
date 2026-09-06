"""Comprehensive tests for ontology module.

Tests cover OBO parsing, ontology traversal, term queries, and hierarchy operations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.ontology.core.go import load_go_obo, write_go_summary
from metainformant.ontology.core.obo import parse_obo
from metainformant.ontology.core.types import Ontology, Term
from metainformant.ontology.query.query import ancestors, descendants, get_subontology


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
        assert len(onto) == 0

        # Create with terms
        term = Term(id="GO:001", name="test")
        onto = Ontology(terms={"GO:001": term})
        assert len(onto) == 1
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
        sub_onto = get_subontology(onto, ["GO:001"])

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

        assert len(onto) == 2
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

        assert len(onto) == 1
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
        assert summary["term_count"] == 1


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_ancestors_missing_term(self):
        """Test ancestors for non-existent term raises ValueError."""
        from metainformant.core.utils.errors import TermNotFoundError

        onto = Ontology()
        with pytest.raises(TermNotFoundError, match="not found in ontology"):
            ancestors(onto, "GO:9999999")

    def test_descendants_missing_term(self):
        """Test descendants for non-existent term raises ValueError."""
        from metainformant.core.utils.errors import TermNotFoundError

        onto = Ontology()
        with pytest.raises(TermNotFoundError, match="not found in ontology"):
            descendants(onto, "GO:9999999")

    def test_empty_ontology(self):
        """Test operations on empty ontology."""
        onto = Ontology()

        from metainformant.core.utils.errors import TermNotFoundError

        assert len(onto) == 0
        assert not onto.has_term("GO:0008150")

        with pytest.raises(TermNotFoundError, match="not found in ontology"):
            ancestors(onto, "GO:0008150")


class TestGoEnrichment:
    """Tests for enrich_genes in core.go."""

    def _annotations(self):
        return {
            "GO:0001": {f"GENE{i}" for i in range(1, 11)},
            "GO:0002": {f"GENE{i}" for i in range(11, 31)},
        }

    def test_fisher_enrichment_finds_enriched_term(self):
        from metainformant.ontology.core.go import enrich_genes

        genes = [f"GENE{i}" for i in range(1, 6)]
        background = [f"GENE{i}" for i in range(1, 51)]
        df = enrich_genes(genes, background, self._annotations(), method="fisher")
        assert not df.empty
        top = df.iloc[0]
        assert top["go_term"] == "GO:0001"
        assert top["observed"] == 5
        assert top["p_value"] < 0.05
        assert 0.0 <= top["bonferroni_p"] <= 1.0

    def test_hypergeometric_method_matches_fisher_direction(self):
        from metainformant.ontology.core.go import enrich_genes

        genes = [f"GENE{i}" for i in range(1, 6)]
        background = [f"GENE{i}" for i in range(1, 51)]
        df = enrich_genes(genes, background, self._annotations(), method="hypergeometric")
        assert not df.empty
        assert df.iloc[0]["go_term"] == "GO:0001"

    def test_invalid_method_raises(self):
        from metainformant.ontology.core.go import enrich_genes

        with pytest.raises(ValueError, match="Unsupported method"):
            enrich_genes(["GENE1"], None, self._annotations(), method="chi2")

    def test_genes_outside_background_filtered(self):
        from metainformant.ontology.core.go import enrich_genes

        df = enrich_genes(["NOGENE1", "NOGENE2"], ["GENE1"], self._annotations())
        assert df.empty


class TestSemanticSimilarity:
    """Tests for semantic_similarity in core.go."""

    def _inputs(self):
        term_ic = {"A": 0.1, "B": 1.0, "C": 2.0}
        hierarchy = {"B": {"A"}, "C": {"A", "B"}}
        return term_ic, hierarchy

    def test_resnik_is_mica_ic(self):
        from metainformant.ontology.core.go import semantic_similarity

        term_ic, hierarchy = self._inputs()
        assert semantic_similarity("B", "C", term_ic, hierarchy, method="resnik") == pytest.approx(1.0)

    def test_lin(self):
        from metainformant.ontology.core.go import semantic_similarity

        term_ic, hierarchy = self._inputs()
        assert semantic_similarity("B", "C", term_ic, hierarchy, method="lin") == pytest.approx(2 / 3)

    def test_jiang_conrath_bounded(self):
        from metainformant.ontology.core.go import semantic_similarity

        term_ic, hierarchy = self._inputs()
        sim = semantic_similarity("B", "C", term_ic, hierarchy, method="jiang-conrath")
        assert 0.0 <= sim <= 1.0
        assert sim == pytest.approx(0.75)

    def test_missing_term_returns_zero(self):
        from metainformant.ontology.core.go import semantic_similarity

        term_ic, hierarchy = self._inputs()
        assert semantic_similarity("B", "MISSING", term_ic, hierarchy) == 0.0

    def test_invalid_method_raises(self):
        from metainformant.ontology.core.go import semantic_similarity

        term_ic, hierarchy = self._inputs()
        with pytest.raises(ValueError, match="Unsupported method"):
            semantic_similarity("B", "C", term_ic, hierarchy, method="cosine")


class TestTermIcAndHierarchy:
    """Tests for calculate_term_ic and build_hierarchy_dict."""

    def test_calculate_term_ic(self):
        import math

        from metainformant.ontology.core.go import calculate_term_ic

        onto = Ontology()
        onto.add_term(Term(term_id="A", name="root"))
        onto.add_term(Term(term_id="B", name="child", is_a_parents=["A"]))
        ic = calculate_term_ic(onto, {"B": 1}, total_annotations=2)
        assert ic["B"] == pytest.approx(-math.log2(0.5))
        assert ic["A"] == 0.0

    def test_build_hierarchy_dict(self):
        from metainformant.ontology.core.go import build_hierarchy_dict

        onto = Ontology()
        onto.add_term(Term(term_id="A", name="root"))
        onto.add_term(Term(term_id="B", name="child", is_a_parents=["A"]))
        hierarchy = build_hierarchy_dict(onto)
        assert hierarchy["B"] == {"A"}
        assert hierarchy["A"] == set()


class TestCountGoScripts:
    """Tests for count_go_scripts."""

    def test_counts_go_named_scripts(self, tmp_path: Path):
        from metainformant.ontology.core.go import count_go_scripts

        (tmp_path / "go_analysis.py").write_text("print('go')")
        (tmp_path / "helper.r").write_text("x <- 1")
        (tmp_path / "notes.txt").write_text("not a script")
        assert count_go_scripts(tmp_path) == 1
