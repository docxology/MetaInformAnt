"""Tests for metainformant.core.utils.newick (shared Newick parser)."""

from __future__ import annotations

import pytest

from metainformant.core.utils.newick import from_newick


class TestFromNewick:
    def test_simple_two_leaf_tree_with_root_length(self) -> None:
        tree = from_newick("(A:0.1,B:0.2):0.3;")
        assert tree["A"] is None
        assert tree["B"] is None
        internal = [k for k, v in tree.items() if isinstance(v, dict)]
        assert len(internal) == 1
        assert tree[internal[0]] == {"A": 0.1, "B": 0.2}

    def test_nested_clades(self) -> None:
        tree = from_newick("((A:0.1,B:0.2):0.3,C:0.4);")
        leaves = {k for k, v in tree.items() if v is None}
        assert leaves == {"A", "B", "C"}
        internals = {k: v for k, v in tree.items() if isinstance(v, dict)}
        # Inner clade children and root children
        assert any(v == {"A": 0.1, "B": 0.2} for v in internals.values())
        root_children = next(v for v in internals.values() if "C" in v)
        assert root_children == {"Internal_1": 0.3, "C": 0.4}

    def test_unlabeled_internal_nodes_get_generated_names(self) -> None:
        tree = from_newick("((A,B),(C,D));")
        assert "Internal_0" in tree
        assert "Internal_1" in tree

    def test_missing_branch_length_defaults_to_zero(self) -> None:
        tree = from_newick("(A,B):0.0;")
        internal = next(v for v in tree.values() if isinstance(v, dict))
        assert internal == {"A": 0.0, "B": 0.0}

    def test_no_trailing_semicolon_accepted(self) -> None:
        tree = from_newick("(A:1.0,B:2.0)")
        assert tree["A"] is None
        assert tree["B"] is None

    def test_empty_string_raises(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            from_newick("")

    def test_whitespace_only_raises(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            from_newick("   ")

    def test_semicolon_only_raises(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            from_newick(";")

    def test_unmatched_parenthesis_raises(self) -> None:
        with pytest.raises(ValueError, match="Unmatched"):
            from_newick("(A:0.1,B:0.2;")

    def test_empty_leaf_label_raises(self) -> None:
        with pytest.raises(ValueError, match="Empty leaf label"):
            from_newick("(A:0.1,:0.2);")

    def test_matches_dna_reexport(self) -> None:
        """The dna.phylogeny re-export stays behaviorally identical to core."""
        from metainformant.dna.phylogeny.tree_analysis import from_newick as dna_from_newick

        s = "((Apis:0.05,Bombus:0.06):0.1,Nasonia:0.2);"
        assert dna_from_newick(s) == from_newick(s)
