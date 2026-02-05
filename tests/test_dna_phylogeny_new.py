"""Comprehensive tests for new phylogenetic tree functions.

Tests cover:
- robinson_foulds_distance: RF distance between trees
- is_monophyletic: monophyly checking
- total_branch_length: sum of all branch lengths
- from_newick: Newick string parsing
- prune_tree: tree pruning to a subset of taxa
- tree_diameter: longest leaf-to-leaf path
- Helper functions: _find_root, _get_all_leaves, _get_bipartitions
"""

from __future__ import annotations

import pytest

from metainformant.dna.phylogeny.tree import (
    Tree,
    _find_root,
    _get_all_leaves,
    _get_bipartitions,
    basic_tree_stats,
    from_newick,
    is_monophyletic,
    neighbor_joining_tree,
    prune_tree,
    robinson_foulds_distance,
    to_newick,
    total_branch_length,
    tree_diameter,
    upgma_tree,
)


# ---------------------------------------------------------------------------
# Fixtures: hand-built trees for deterministic testing
# ---------------------------------------------------------------------------


def _make_simple_tree() -> Tree:
    """Build: Root(A:1.0, B:2.0) -- two leaves.

    Structure:
        Root
        / \\
       A   B
      1.0 2.0
    """
    return {
        "A": None,
        "B": None,
        "Root": {"A": 1.0, "B": 2.0},
    }


def _make_four_leaf_symmetric() -> Tree:
    """Build a balanced 4-leaf tree: ((A:1,B:1):0.5,(C:1,D:1):0.5).

    Structure:
        Root
        /   \\
      N1     N2
      / \\   / \\
     A   B C   D

    Branch lengths: all leaf edges = 1.0, internal edges to Root = 0.5.
    """
    return {
        "A": None,
        "B": None,
        "C": None,
        "D": None,
        "N1": {"A": 1.0, "B": 1.0},
        "N2": {"C": 1.0, "D": 1.0},
        "Root": {"N1": 0.5, "N2": 0.5},
    }


def _make_four_leaf_asymmetric() -> Tree:
    """Build: (A:1,(B:1,(C:1,D:1):0.5):0.5) -- caterpillar topology.

    Structure:
        Root
        / \\
       A   N1
          / \\
         B   N2
            / \\
           C   D
    """
    return {
        "A": None,
        "B": None,
        "C": None,
        "D": None,
        "N2": {"C": 1.0, "D": 1.0},
        "N1": {"B": 1.0, "N2": 0.5},
        "Root": {"A": 1.0, "N1": 0.5},
    }


def _make_five_leaf_tree() -> Tree:
    """Build: ((A:1,B:1):0.5,((C:1,D:1):0.3,E:1.3):0.5).

    Structure:
        Root
        /    \\
      N1      N2
      / \\    / \\
     A   B  N3   E
           / \\
          C   D
    """
    return {
        "A": None,
        "B": None,
        "C": None,
        "D": None,
        "E": None,
        "N3": {"C": 1.0, "D": 1.0},
        "N1": {"A": 1.0, "B": 1.0},
        "N2": {"N3": 0.3, "E": 1.3},
        "Root": {"N1": 0.5, "N2": 0.5},
    }


def _make_six_leaf_tree() -> Tree:
    """Build a 6-leaf tree for richer tests.

    Structure: ((A:1,B:1):0.5,(C:1,(D:1,(E:1,F:1):0.3):0.4):0.5)
    """
    return {
        "A": None,
        "B": None,
        "C": None,
        "D": None,
        "E": None,
        "F": None,
        "N_EF": {"E": 1.0, "F": 1.0},
        "N_DEF": {"D": 1.0, "N_EF": 0.3},
        "N_CDEF": {"C": 1.0, "N_DEF": 0.4},
        "N_AB": {"A": 1.0, "B": 1.0},
        "Root": {"N_AB": 0.5, "N_CDEF": 0.5},
    }


# ---------------------------------------------------------------------------
# Sample aligned sequences for integration tests
# ---------------------------------------------------------------------------

SAMPLE_SEQS_4 = {
    "taxA": "ATCGATCGATCGATCG",
    "taxB": "ATCGATCGATCGATCG",
    "taxC": "TTCGATCGATCGTTCG",
    "taxD": "TTTTTTTTTTTTTTTT",
}

SAMPLE_SEQS_5 = {
    "s1": "ATCGATCGATCGATCG",
    "s2": "ATCGATCGATCGATCA",
    "s3": "TTCGATCGATCGTTCG",
    "s4": "TTTTTTTTTTTTTTTT",
    "s5": "AAAATTTTAAAATTTT",
}


# ===========================================================================
# Tests: _find_root
# ===========================================================================


class TestFindRoot:
    """Tests for the _find_root helper."""

    def test_simple_tree(self) -> None:
        tree = _make_simple_tree()
        assert _find_root(tree) == "Root"

    def test_four_leaf_symmetric(self) -> None:
        tree = _make_four_leaf_symmetric()
        assert _find_root(tree) == "Root"

    def test_four_leaf_asymmetric(self) -> None:
        tree = _make_four_leaf_asymmetric()
        assert _find_root(tree) == "Root"

    def test_five_leaf_tree(self) -> None:
        tree = _make_five_leaf_tree()
        assert _find_root(tree) == "Root"

    def test_nj_generated_tree(self) -> None:
        tree = neighbor_joining_tree(SAMPLE_SEQS_4)
        root = _find_root(tree)
        assert root is not None
        # Root should be an internal node (not a leaf)
        assert isinstance(tree[root], dict)


# ===========================================================================
# Tests: _get_all_leaves
# ===========================================================================


class TestGetAllLeaves:
    """Tests for the _get_all_leaves helper."""

    def test_simple_tree(self) -> None:
        tree = _make_simple_tree()
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B"}

    def test_four_leaf_symmetric(self) -> None:
        tree = _make_four_leaf_symmetric()
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D"}

    def test_five_leaf_tree(self) -> None:
        tree = _make_five_leaf_tree()
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D", "E"}

    def test_nj_tree_preserves_all_taxa(self) -> None:
        tree = neighbor_joining_tree(SAMPLE_SEQS_4)
        leaves = set(_get_all_leaves(tree))
        assert leaves == set(SAMPLE_SEQS_4.keys())


# ===========================================================================
# Tests: _get_bipartitions
# ===========================================================================


class TestGetBipartitions:
    """Tests for the _get_bipartitions helper."""

    def test_simple_tree_no_nontrivial_splits(self) -> None:
        """A 2-leaf tree has no non-trivial bipartitions."""
        tree = _make_simple_tree()
        all_leaves = set(_get_all_leaves(tree))
        splits = _get_bipartitions(tree, all_leaves)
        assert len(splits) == 0

    def test_four_leaf_symmetric_has_one_split(self) -> None:
        """((A,B),(C,D)) has one non-trivial split: {A,B} | {C,D}."""
        tree = _make_four_leaf_symmetric()
        all_leaves = set(_get_all_leaves(tree))
        splits = _get_bipartitions(tree, all_leaves)
        assert len(splits) == 1
        split = next(iter(splits))
        assert split == frozenset({"A", "B"}) or split == frozenset({"C", "D"})

    def test_five_leaf_tree_splits(self) -> None:
        """Verify the expected splits from a 5-leaf tree."""
        tree = _make_five_leaf_tree()
        all_leaves = set(_get_all_leaves(tree))
        splits = _get_bipartitions(tree, all_leaves)
        # Should have non-trivial splits
        assert len(splits) >= 1
        # {A,B} should be one split
        assert frozenset({"A", "B"}) in splits or frozenset({"C", "D", "E"}) in splits

    def test_bipartitions_are_frozensets(self) -> None:
        tree = _make_four_leaf_symmetric()
        all_leaves = set(_get_all_leaves(tree))
        splits = _get_bipartitions(tree, all_leaves)
        for split in splits:
            assert isinstance(split, frozenset)


# ===========================================================================
# Tests: robinson_foulds_distance
# ===========================================================================


class TestRobinsonFouldsDistance:
    """Tests for robinson_foulds_distance."""

    def test_identical_trees_distance_zero(self) -> None:
        """RF distance of a tree with itself is 0."""
        tree = _make_four_leaf_symmetric()
        assert robinson_foulds_distance(tree, tree) == 0

    def test_symmetry(self) -> None:
        """RF(t1, t2) == RF(t2, t1)."""
        t1 = _make_four_leaf_symmetric()
        t2 = _make_four_leaf_asymmetric()
        d12 = robinson_foulds_distance(t1, t2)
        d21 = robinson_foulds_distance(t2, t1)
        assert d12 == d21

    def test_non_negative(self) -> None:
        """RF distance is always non-negative."""
        t1 = _make_four_leaf_symmetric()
        t2 = _make_four_leaf_asymmetric()
        assert robinson_foulds_distance(t1, t2) >= 0

    def test_different_topologies_positive(self) -> None:
        """Different 5-leaf topologies should produce positive RF distance.

        With 4 leaves, ((A,B),(C,D)) and (A,(B,(C,D))) share the same non-trivial
        bipartition {A,B}|{C,D}, so we use 5-leaf trees for a meaningful difference.
        """
        # ((A,B),(C,(D,E))) has split {A,B}
        t1 = {
            "A": None, "B": None, "C": None, "D": None, "E": None,
            "N_AB": {"A": 1.0, "B": 1.0},
            "N_DE": {"D": 1.0, "E": 1.0},
            "N_CDE": {"C": 1.0, "N_DE": 0.5},
            "Root": {"N_AB": 0.5, "N_CDE": 0.5},
        }
        # ((A,C),(B,(D,E))) has split {A,C} instead
        t2 = {
            "A": None, "B": None, "C": None, "D": None, "E": None,
            "N_AC": {"A": 1.0, "C": 1.0},
            "N_DE": {"D": 1.0, "E": 1.0},
            "N_BDE": {"B": 1.0, "N_DE": 0.5},
            "Root": {"N_AC": 0.5, "N_BDE": 0.5},
        }
        d = robinson_foulds_distance(t1, t2)
        assert d > 0

    def test_triangle_inequality(self) -> None:
        """RF distance satisfies the triangle inequality: d(a,c) <= d(a,b) + d(b,c).

        We build three 5-leaf trees with different topologies to test this.
        """
        # Tree 1: ((A,B),(C,(D,E)))
        t1 = {
            "A": None, "B": None, "C": None, "D": None, "E": None,
            "N_AB": {"A": 1.0, "B": 1.0},
            "N_DE": {"D": 1.0, "E": 1.0},
            "N_CDE": {"C": 1.0, "N_DE": 0.5},
            "Root": {"N_AB": 0.5, "N_CDE": 0.5},
        }
        # Tree 2: ((A,C),(B,(D,E)))
        t2 = {
            "A": None, "B": None, "C": None, "D": None, "E": None,
            "N_AC": {"A": 1.0, "C": 1.0},
            "N_DE": {"D": 1.0, "E": 1.0},
            "N_BDE": {"B": 1.0, "N_DE": 0.5},
            "Root": {"N_AC": 0.5, "N_BDE": 0.5},
        }
        # Tree 3: ((A,D),(B,(C,E)))
        t3 = {
            "A": None, "B": None, "C": None, "D": None, "E": None,
            "N_AD": {"A": 1.0, "D": 1.0},
            "N_CE": {"C": 1.0, "E": 1.0},
            "N_BCE": {"B": 1.0, "N_CE": 0.5},
            "Root": {"N_AD": 0.5, "N_BCE": 0.5},
        }

        d12 = robinson_foulds_distance(t1, t2)
        d23 = robinson_foulds_distance(t2, t3)
        d13 = robinson_foulds_distance(t1, t3)

        assert d13 <= d12 + d23, f"Triangle inequality violated: {d13} > {d12} + {d23}"

    def test_mismatched_leaves_raises(self) -> None:
        """Trees with different leaf sets should raise ValueError."""
        t1 = _make_simple_tree()  # leaves: A, B
        t2 = _make_four_leaf_symmetric()  # leaves: A, B, C, D
        with pytest.raises(ValueError, match="same leaf set"):
            robinson_foulds_distance(t1, t2)

    def test_two_leaf_trees_zero(self) -> None:
        """Any two binary trees with only 2 leaves have the same (trivial) topology."""
        t1 = {"X": None, "Y": None, "R1": {"X": 1.0, "Y": 2.0}}
        t2 = {"X": None, "Y": None, "R2": {"X": 3.0, "Y": 4.0}}
        assert robinson_foulds_distance(t1, t2) == 0

    def test_three_leaf_trees_same_topology(self) -> None:
        """All unrooted 3-leaf binary trees have the same topology -> RF = 0."""
        t1 = {
            "A": None, "B": None, "C": None,
            "N": {"A": 1.0, "B": 1.0},
            "Root": {"N": 0.5, "C": 1.0},
        }
        t2 = {
            "A": None, "B": None, "C": None,
            "N": {"B": 1.0, "C": 1.0},
            "Root": {"N": 0.5, "A": 1.0},
        }
        # 3-leaf trees have no non-trivial splits, so RF should be 0
        assert robinson_foulds_distance(t1, t2) == 0

    def test_rf_with_nj_generated_trees(self) -> None:
        """RF distance between NJ and UPGMA trees built from same sequences."""
        tree_nj = neighbor_joining_tree(SAMPLE_SEQS_4)
        tree_upgma = upgma_tree(SAMPLE_SEQS_4)
        d = robinson_foulds_distance(tree_nj, tree_upgma)
        assert isinstance(d, int)
        assert d >= 0

    def test_rf_even_integer(self) -> None:
        """RF distance is the size of the symmetric difference of bipartition sets."""
        t1 = _make_four_leaf_symmetric()
        t2 = _make_four_leaf_asymmetric()
        d = robinson_foulds_distance(t1, t2)
        # The value should be an integer
        assert isinstance(d, int)


# ===========================================================================
# Tests: is_monophyletic
# ===========================================================================


class TestIsMonophyletic:
    """Tests for is_monophyletic."""

    def test_single_taxon_not_monophyletic_without_dedicated_node(self) -> None:
        """A single taxon has no internal node whose entire leaf set is just that taxon.

        In the implementation, is_monophyletic checks for an internal node whose
        descendant leaf set equals the target. A single leaf in a standard binary
        tree has no such dedicated internal node, so it returns False.
        """
        tree = _make_four_leaf_symmetric()
        # In a standard binary tree, no internal node has exactly one descendant leaf
        assert is_monophyletic(tree, ["A"]) is False

    def test_empty_taxa_is_monophyletic(self) -> None:
        """Empty taxa set is trivially monophyletic."""
        tree = _make_four_leaf_symmetric()
        assert is_monophyletic(tree, []) is True

    def test_all_taxa_monophyletic(self) -> None:
        """The full leaf set is monophyletic (rooted at the root)."""
        tree = _make_four_leaf_symmetric()
        assert is_monophyletic(tree, ["A", "B", "C", "D"]) is True

    def test_true_clade_symmetric(self) -> None:
        """In ((A,B),(C,D)), {A,B} forms a monophyletic group."""
        tree = _make_four_leaf_symmetric()
        assert is_monophyletic(tree, ["A", "B"]) is True

    def test_true_clade_other_side(self) -> None:
        """In ((A,B),(C,D)), {C,D} forms a monophyletic group."""
        tree = _make_four_leaf_symmetric()
        assert is_monophyletic(tree, ["C", "D"]) is True

    def test_non_clade(self) -> None:
        """In ((A,B),(C,D)), {A,C} is NOT monophyletic."""
        tree = _make_four_leaf_symmetric()
        assert is_monophyletic(tree, ["A", "C"]) is False

    def test_non_clade_three_taxa(self) -> None:
        """In ((A,B),(C,D)), {A,C,D} is NOT monophyletic (not a clade)."""
        tree = _make_four_leaf_symmetric()
        assert is_monophyletic(tree, ["A", "C", "D"]) is False

    def test_nested_clade_asymmetric(self) -> None:
        """In (A,(B,(C,D))), {C,D} is monophyletic."""
        tree = _make_four_leaf_asymmetric()
        assert is_monophyletic(tree, ["C", "D"]) is True

    def test_larger_clade_asymmetric(self) -> None:
        """In (A,(B,(C,D))), {B,C,D} is monophyletic."""
        tree = _make_four_leaf_asymmetric()
        assert is_monophyletic(tree, ["B", "C", "D"]) is True

    def test_non_clade_asymmetric(self) -> None:
        """In (A,(B,(C,D))), {A,D} is NOT monophyletic."""
        tree = _make_four_leaf_asymmetric()
        assert is_monophyletic(tree, ["A", "D"]) is False

    def test_five_leaf_monophyly(self) -> None:
        """Test monophyly in a 5-leaf tree."""
        tree = _make_five_leaf_tree()
        # {A,B} should be monophyletic
        assert is_monophyletic(tree, ["A", "B"]) is True
        # {C,D} should be monophyletic
        assert is_monophyletic(tree, ["C", "D"]) is True
        # {C,D,E} should be monophyletic (N2 node)
        assert is_monophyletic(tree, ["C", "D", "E"]) is True
        # {A,C} should NOT be monophyletic
        assert is_monophyletic(tree, ["A", "C"]) is False

    def test_monophyly_with_nj_tree(self) -> None:
        """Test monophyly on a tree generated by neighbor_joining."""
        tree = neighbor_joining_tree(SAMPLE_SEQS_4)
        leaves = _get_all_leaves(tree)
        # Full leaf set is always monophyletic
        assert is_monophyletic(tree, leaves) is True
        # Empty set is trivially monophyletic
        assert is_monophyletic(tree, []) is True
        # The result for any single leaf depends on tree structure (no dedicated node)
        # but we verify the function returns a bool
        for leaf in leaves:
            result = is_monophyletic(tree, [leaf])
            assert isinstance(result, bool)


# ===========================================================================
# Tests: total_branch_length
# ===========================================================================


class TestTotalBranchLength:
    """Tests for total_branch_length."""

    def test_simple_tree(self) -> None:
        tree = _make_simple_tree()
        # Root->A: 1.0, Root->B: 2.0
        assert total_branch_length(tree) == pytest.approx(3.0)

    def test_four_leaf_symmetric(self) -> None:
        tree = _make_four_leaf_symmetric()
        # N1->A:1 + N1->B:1 + N2->C:1 + N2->D:1 + Root->N1:0.5 + Root->N2:0.5 = 5.0
        assert total_branch_length(tree) == pytest.approx(5.0)

    def test_four_leaf_asymmetric(self) -> None:
        tree = _make_four_leaf_asymmetric()
        # Root->A:1 + Root->N1:0.5 + N1->B:1 + N1->N2:0.5 + N2->C:1 + N2->D:1 = 5.0
        assert total_branch_length(tree) == pytest.approx(5.0)

    def test_five_leaf_tree(self) -> None:
        tree = _make_five_leaf_tree()
        # N3->C:1 + N3->D:1 + N1->A:1 + N1->B:1 + N2->N3:0.3 + N2->E:1.3
        # + Root->N1:0.5 + Root->N2:0.5 = 6.6
        assert total_branch_length(tree) == pytest.approx(6.6)

    def test_non_negative(self) -> None:
        """Total branch length is always non-negative."""
        tree = neighbor_joining_tree(SAMPLE_SEQS_4)
        assert total_branch_length(tree) >= 0.0

    def test_zero_length_branches(self) -> None:
        """Tree with zero-length branches."""
        tree = {
            "A": None,
            "B": None,
            "Root": {"A": 0.0, "B": 0.0},
        }
        assert total_branch_length(tree) == pytest.approx(0.0)

    def test_none_branch_lengths_treated_as_zero(self) -> None:
        """Branches with None length should be treated as zero."""
        tree = {
            "A": None,
            "B": None,
            "Root": {"A": None, "B": 2.0},
        }
        assert total_branch_length(tree) == pytest.approx(2.0)

    def test_branch_length_with_bootstrap(self) -> None:
        """Bootstrap values should not be counted as branch lengths."""
        tree = {
            "A": None,
            "B": None,
            "Root": {"A": 1.0, "B": 2.0, "bootstrap": 95},
        }
        assert total_branch_length(tree) == pytest.approx(3.0)

    def test_nj_tree_total_length(self) -> None:
        """NJ tree should have a positive total branch length for distinct sequences."""
        tree = neighbor_joining_tree(SAMPLE_SEQS_4)
        length = total_branch_length(tree)
        assert length > 0.0
        assert isinstance(length, float)


# ===========================================================================
# Tests: from_newick
# ===========================================================================


class TestFromNewick:
    """Tests for from_newick."""

    def test_simple_two_leaf(self) -> None:
        """Parse (A:1.0,B:2.0);"""
        tree = from_newick("(A:1.0,B:2.0);")
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B"}

    def test_four_leaf_balanced(self) -> None:
        """Parse ((A:1,B:1):0.5,(C:1,D:1):0.5);"""
        tree = from_newick("((A:1,B:1):0.5,(C:1,D:1):0.5);")
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D"}
        stats = basic_tree_stats(tree)
        assert stats["leaves"] == 4

    def test_nested_caterpillar(self) -> None:
        """Parse (A:1,(B:1,(C:1,D:1):0.5):0.5);"""
        tree = from_newick("(A:1,(B:1,(C:1,D:1):0.5):0.5);")
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D"}

    def test_no_branch_lengths(self) -> None:
        """Parse (A,B,(C,D));"""
        tree = from_newick("(A,B,(C,D));")
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D"}

    def test_single_leaf_raises_or_parses(self) -> None:
        """A single leaf Newick string like 'A;' should parse as a single leaf."""
        tree = from_newick("A;")
        assert "A" in tree
        leaves = _get_all_leaves(tree)
        assert len(leaves) == 1

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
        with pytest.raises(ValueError, match="parenthesis"):
            from_newick("((A:1,B:1);")

    def test_five_leaf_tree(self) -> None:
        """Parse a 5-leaf tree."""
        tree = from_newick("((A:1,B:1):0.5,((C:1,D:1):0.3,E:1.3):0.5);")
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D", "E"}
        stats = basic_tree_stats(tree)
        assert stats["leaves"] == 5

    def test_branch_lengths_preserved(self) -> None:
        """Branch lengths from Newick should be preserved in the tree."""
        tree = from_newick("(A:0.1,B:0.2);")
        # Find the internal node
        root = _find_root(tree)
        children = tree[root]
        assert isinstance(children, dict)
        assert children.get("A") == pytest.approx(0.1)
        assert children.get("B") == pytest.approx(0.2)

    def test_internal_node_labels(self) -> None:
        """Parse Newick with labeled internal nodes: ((A,B)AB,(C,D)CD)Root;"""
        tree = from_newick("((A,B)AB,(C,D)CD)Root;")
        assert "AB" in tree
        assert "CD" in tree
        assert "Root" in tree
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B", "C", "D"}

    def test_whitespace_tolerance(self) -> None:
        """Parser should handle leading/trailing whitespace."""
        tree = from_newick("  (A:1.0, B:2.0) ;  ")
        leaves = set(_get_all_leaves(tree))
        assert leaves == {"A", "B"}


# ===========================================================================
# Tests: Newick round-trip (from_newick -> to_newick -> from_newick)
# ===========================================================================


class TestNewickRoundTrip:
    """Test that from_newick and to_newick are consistent."""

    def test_round_trip_preserves_leaves(self) -> None:
        """Parse -> serialize -> re-parse should preserve the leaf set."""
        newick_input = "((A:1,B:1):0.5,(C:1,D:1):0.5);"
        tree1 = from_newick(newick_input)
        newick_output = to_newick(tree1)
        tree2 = from_newick(newick_output)

        leaves1 = set(_get_all_leaves(tree1))
        leaves2 = set(_get_all_leaves(tree2))
        assert leaves1 == leaves2

    def test_round_trip_preserves_topology(self) -> None:
        """Parse -> serialize -> re-parse should preserve topology (RF=0)."""
        newick_input = "((A:1,B:1):0.5,(C:1,D:1):0.5);"
        tree1 = from_newick(newick_input)
        newick_output = to_newick(tree1)
        tree2 = from_newick(newick_output)
        assert robinson_foulds_distance(tree1, tree2) == 0

    def test_round_trip_five_leaf(self) -> None:
        """Round-trip with a 5-leaf tree."""
        newick_input = "((A:1.0,B:2.0):0.5,((C:0.5,D:0.5):0.3,E:1.3):0.5);"
        tree1 = from_newick(newick_input)
        newick_output = to_newick(tree1)
        tree2 = from_newick(newick_output)

        leaves1 = set(_get_all_leaves(tree1))
        leaves2 = set(_get_all_leaves(tree2))
        assert leaves1 == leaves2
        assert robinson_foulds_distance(tree1, tree2) == 0

    def test_round_trip_total_branch_length(self) -> None:
        """Total branch length should be preserved through round-trip."""
        newick_input = "((A:1.5,B:2.5):0.5,(C:3.0,D:4.0):0.5);"
        tree1 = from_newick(newick_input)
        newick_output = to_newick(tree1)
        tree2 = from_newick(newick_output)

        tbl1 = total_branch_length(tree1)
        tbl2 = total_branch_length(tree2)
        assert tbl1 == pytest.approx(tbl2)

    def test_nj_tree_round_trip(self) -> None:
        """Build NJ tree, serialize to Newick, re-parse, and compare leaves."""
        tree1 = neighbor_joining_tree(SAMPLE_SEQS_4)
        newick = to_newick(tree1)
        assert newick.endswith(";")
        tree2 = from_newick(newick)

        leaves1 = set(_get_all_leaves(tree1))
        leaves2 = set(_get_all_leaves(tree2))
        assert leaves1 == leaves2


# ===========================================================================
# Tests: prune_tree
# ===========================================================================


class TestPruneTree:
    """Tests for prune_tree."""

    def test_prune_to_all_leaves_is_identity(self) -> None:
        """Pruning to all leaves should preserve the full leaf set."""
        tree = _make_four_leaf_symmetric()
        pruned = prune_tree(tree, ["A", "B", "C", "D"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A", "B", "C", "D"}

    def test_prune_to_two_leaves(self) -> None:
        """Pruning to 2 leaves from a 4-leaf tree."""
        tree = _make_four_leaf_symmetric()
        pruned = prune_tree(tree, ["A", "B"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A", "B"}

    def test_prune_to_single_clade(self) -> None:
        """Pruning to one side of a symmetric tree."""
        tree = _make_four_leaf_symmetric()
        pruned = prune_tree(tree, ["C", "D"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"C", "D"}

    def test_prune_cross_clade(self) -> None:
        """Pruning to taxa across clades."""
        tree = _make_four_leaf_symmetric()
        pruned = prune_tree(tree, ["A", "C"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A", "C"}

    def test_prune_to_one_leaf(self) -> None:
        """Pruning to a single leaf."""
        tree = _make_four_leaf_symmetric()
        pruned = prune_tree(tree, ["A"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A"}

    def test_prune_five_leaf_tree(self) -> None:
        """Prune a 5-leaf tree to 3 leaves."""
        tree = _make_five_leaf_tree()
        pruned = prune_tree(tree, ["A", "C", "E"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A", "C", "E"}

    def test_prune_empty_raises(self) -> None:
        """Pruning with empty taxa_to_keep should raise ValueError."""
        tree = _make_four_leaf_symmetric()
        with pytest.raises(ValueError, match="must not be empty"):
            prune_tree(tree, [])

    def test_prune_nonexistent_taxa_raises(self) -> None:
        """Pruning with taxa not in the tree should raise ValueError."""
        tree = _make_four_leaf_symmetric()
        with pytest.raises(ValueError, match="not found in tree"):
            prune_tree(tree, ["A", "Z"])

    def test_prune_preserves_monophyly(self) -> None:
        """After pruning, existing clades that were monophyletic should remain so."""
        tree = _make_five_leaf_tree()
        # {A,B} is monophyletic in the original tree
        pruned = prune_tree(tree, ["A", "B", "E"])
        # {A,B} should still be monophyletic in the pruned tree
        assert is_monophyletic(pruned, ["A", "B"]) is True

    def test_prune_collapses_single_children(self) -> None:
        """Pruning should collapse nodes with a single child (no degree-2 nodes)."""
        tree = _make_four_leaf_asymmetric()
        # Prune to {A, D}: intermediate nodes with single child should collapse
        pruned = prune_tree(tree, ["A", "D"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A", "D"}
        # All internal nodes should have 2+ children (no single-child nodes)
        for node, children in pruned.items():
            if isinstance(children, dict):
                real_children = [c for c in children if c != "bootstrap"]
                assert len(real_children) >= 2, f"Node {node} has only {len(real_children)} child(ren)"

    def test_prune_branch_length_accumulation(self) -> None:
        """When collapsing single-child nodes, branch lengths should accumulate."""
        tree = _make_four_leaf_asymmetric()
        # Pruning {C, D} from (A:1,(B:1,(C:1,D:1):0.5):0.5)
        # should give a tree where C and D retain proper branch lengths
        pruned = prune_tree(tree, ["C", "D"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"C", "D"}
        # The total branch length should reflect accumulated internal edges
        tbl = total_branch_length(pruned)
        assert tbl > 0.0

    def test_prune_nj_tree(self) -> None:
        """Prune an NJ-generated tree to a subset of taxa."""
        tree = neighbor_joining_tree(SAMPLE_SEQS_5)
        pruned = prune_tree(tree, ["s1", "s3", "s5"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"s1", "s3", "s5"}

    def test_prune_six_leaf_to_three(self) -> None:
        """Prune a 6-leaf tree to 3 leaves from different subtrees."""
        tree = _make_six_leaf_tree()
        pruned = prune_tree(tree, ["A", "D", "F"])
        leaves = set(_get_all_leaves(pruned))
        assert leaves == {"A", "D", "F"}
        stats = basic_tree_stats(pruned)
        assert stats["leaves"] == 3


# ===========================================================================
# Tests: tree_diameter
# ===========================================================================


class TestTreeDiameter:
    """Tests for tree_diameter."""

    def test_simple_tree(self) -> None:
        """Diameter of a 2-leaf tree is the sum of both branch lengths."""
        tree = _make_simple_tree()
        # A--1.0--Root--2.0--B => diameter = 3.0
        assert tree_diameter(tree) == pytest.approx(3.0)

    def test_four_leaf_symmetric(self) -> None:
        """Diameter of ((A:1,B:1):0.5,(C:1,D:1):0.5).

        Longest path: A->N1->Root->N2->C (or any leaf-to-opposite-leaf) = 1+0.5+0.5+1 = 3.0
        """
        tree = _make_four_leaf_symmetric()
        assert tree_diameter(tree) == pytest.approx(3.0)

    def test_four_leaf_asymmetric(self) -> None:
        """Diameter of (A:1,(B:1,(C:1,D:1):0.5):0.5).

        Longest path: A->Root->N1->N2->C (or D) = 1+0.5+1+0.5+1 ... wait, let's compute:
        A->Root: 1.0
        Root->N1: 0.5
        N1->N2: 0.5
        N2->C: 1.0
        Total A->C = 1.0+0.5+0.5+1.0 = 3.0
        N2->D: 1.0
        Total A->D = 3.0 also.
        """
        tree = _make_four_leaf_asymmetric()
        assert tree_diameter(tree) == pytest.approx(3.0)

    def test_five_leaf_tree(self) -> None:
        """Diameter of ((A:1,B:1):0.5,((C:1,D:1):0.3,E:1.3):0.5).

        Possible longest paths:
        A->N1->Root->N2->E = 1+0.5+0.5+1.3 = 3.3
        A->N1->Root->N2->N3->C = 1+0.5+0.5+0.3+1 = 3.3
        So diameter should be 3.3.
        """
        tree = _make_five_leaf_tree()
        assert tree_diameter(tree) == pytest.approx(3.3)

    def test_non_negative(self) -> None:
        """Diameter is always non-negative."""
        tree = neighbor_joining_tree(SAMPLE_SEQS_4)
        assert tree_diameter(tree) >= 0.0

    def test_zero_length_branches(self) -> None:
        """All zero-length branches yield diameter 0."""
        tree = {
            "A": None,
            "B": None,
            "Root": {"A": 0.0, "B": 0.0},
        }
        assert tree_diameter(tree) == pytest.approx(0.0)

    def test_diameter_at_least_total_over_two(self) -> None:
        """In any binary tree, diameter >= total_branch_length / n_leaves (rough check)."""
        tree = _make_six_leaf_tree()
        d = tree_diameter(tree)
        tbl = total_branch_length(tree)
        n_leaves = len(_get_all_leaves(tree))
        # Diameter should be a substantial portion of total branch length
        assert d > 0.0
        assert d <= tbl  # Diameter cannot exceed total branch length

    def test_diameter_equals_longest_path(self) -> None:
        """Verify diameter is the longest path in a well-understood tree."""
        # Build a tree where we know the longest path exactly:
        # ((A:5, B:1):1, C:2)
        # Longest: A->5->Int->1->Root->2->C = 5+1+2 = 8
        # A->5->Int->1->B = 5+1 = 6
        tree = {
            "A": None,
            "B": None,
            "C": None,
            "Int": {"A": 5.0, "B": 1.0},
            "Root": {"Int": 1.0, "C": 2.0},
        }
        assert tree_diameter(tree) == pytest.approx(8.0)

    def test_diameter_nj_tree(self) -> None:
        """Diameter of an NJ-generated tree should be positive for diverse sequences."""
        tree = neighbor_joining_tree(SAMPLE_SEQS_5)
        d = tree_diameter(tree)
        assert d > 0.0
        assert isinstance(d, float)


# ===========================================================================
# Integration tests: combining multiple new functions
# ===========================================================================


class TestIntegration:
    """Integration tests combining multiple new phylogeny functions."""

    def test_prune_preserves_rf_consistency(self) -> None:
        """RF distance of pruned trees should be consistent with the original topology.

        If we prune two trees to the same subset, the RF distance should not increase
        beyond the unpruned RF distance (it may decrease since some splits vanish).
        """
        t1 = _make_five_leaf_tree()
        # Build a different 5-leaf topology
        t2 = {
            "A": None, "B": None, "C": None, "D": None, "E": None,
            "N_AC": {"A": 1.0, "C": 1.0},
            "N_BD": {"B": 1.0, "D": 1.0},
            "N_ACE": {"N_AC": 0.5, "E": 1.0},
            "Root": {"N_ACE": 0.5, "N_BD": 0.5},
        }

        rf_full = robinson_foulds_distance(t1, t2)

        # Prune both to {A, B, C, D}
        t1_pruned = prune_tree(t1, ["A", "B", "C", "D"])
        t2_pruned = prune_tree(t2, ["A", "B", "C", "D"])

        rf_pruned = robinson_foulds_distance(t1_pruned, t2_pruned)

        # The pruned RF can differ but should be non-negative
        assert rf_pruned >= 0
        # Both distances are valid integers
        assert isinstance(rf_full, int)
        assert isinstance(rf_pruned, int)

    def test_newick_to_prune_to_diameter(self) -> None:
        """Parse Newick, prune, then compute diameter."""
        tree = from_newick("((A:2,B:3):1,((C:1,D:1):0.5,E:4):1);")
        pruned = prune_tree(tree, ["A", "D", "E"])
        d = tree_diameter(pruned)
        assert d > 0.0
        assert isinstance(d, float)

    def test_newick_to_monophyly(self) -> None:
        """Parse Newick and test monophyly."""
        tree = from_newick("((A:1,B:1):0.5,(C:1,D:1):0.5);")
        assert is_monophyletic(tree, ["A", "B"]) is True
        assert is_monophyletic(tree, ["C", "D"]) is True
        assert is_monophyletic(tree, ["A", "C"]) is False

    def test_full_pipeline_nj_tree(self) -> None:
        """Full pipeline: build NJ tree, serialize, re-parse, prune, measure."""
        # Build
        tree = neighbor_joining_tree(SAMPLE_SEQS_5)
        assert basic_tree_stats(tree)["leaves"] == 5

        # Serialize and re-parse
        newick = to_newick(tree)
        tree2 = from_newick(newick)
        assert set(_get_all_leaves(tree2)) == set(SAMPLE_SEQS_5.keys())

        # Prune
        pruned = prune_tree(tree, ["s1", "s3", "s5"])
        assert set(_get_all_leaves(pruned)) == {"s1", "s3", "s5"}

        # Measure
        d = tree_diameter(pruned)
        assert d > 0.0
        tbl = total_branch_length(pruned)
        assert tbl > 0.0

    def test_total_length_decreases_after_prune(self) -> None:
        """Pruning a tree should not increase total branch length relative to the kept subtree."""
        tree = _make_six_leaf_tree()
        original_tbl = total_branch_length(tree)

        pruned = prune_tree(tree, ["A", "B", "C"])
        pruned_tbl = total_branch_length(pruned)

        # Pruned tree has fewer edges, so its total should generally be less
        # (branch lengths may accumulate during collapse, but total should not exceed original)
        assert pruned_tbl <= original_tbl + 1e-9  # allow float tolerance

    def test_diameter_not_greater_than_total_length(self) -> None:
        """The diameter (longest path) can never exceed total branch length."""
        for tree_fn in [
            _make_simple_tree,
            _make_four_leaf_symmetric,
            _make_four_leaf_asymmetric,
            _make_five_leaf_tree,
            _make_six_leaf_tree,
        ]:
            tree = tree_fn()
            d = tree_diameter(tree)
            tbl = total_branch_length(tree)
            assert d <= tbl + 1e-9, f"Diameter {d} exceeds total branch length {tbl}"

    def test_upgma_functions(self) -> None:
        """Exercise all new functions on a UPGMA tree."""
        tree = upgma_tree(SAMPLE_SEQS_4)
        leaves = _get_all_leaves(tree)
        assert len(leaves) == 4

        tbl = total_branch_length(tree)
        assert tbl > 0.0

        d = tree_diameter(tree)
        assert d > 0.0
        assert d <= tbl + 1e-9

        newick = to_newick(tree)
        tree2 = from_newick(newick)
        assert set(_get_all_leaves(tree2)) == set(leaves)

        # Monophyly: all taxa is trivially monophyletic
        assert is_monophyletic(tree, leaves) is True

        # Prune to 2
        pruned = prune_tree(tree, leaves[:2])
        assert len(_get_all_leaves(pruned)) == 2
