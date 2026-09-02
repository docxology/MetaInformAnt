"""Equivalence and behavior tests for the optimized newick parser fast paths.

The perf-lane optimization (tests/perf/BASELINES.md) replaced per-character
scans in ``from_newick`` with paren-aware ``str.find``/run-skipping scans.
These tests pin the observable behavior that must not change: parse results
on nested, sibling-rich, and adversarial trees, plus error discipline.
"""

from __future__ import annotations

import pytest

from metainformant.core.utils.newick import from_newick


def test_deeply_nested_chain_parses() -> None:
    # Right-leaning chain exercises the nested-paren general path.
    newick = "((((a:0.1,b:0.2):0.3,c:0.4):0.5,d:0.6):0.7,e:0.8);"
    tree = from_newick(newick)
    root_children = tree["Internal_0"]
    assert root_children["e"] == 0.8
    assert set(root_children) == {"Internal_1", "e"}


def test_wide_single_level_clade_parses() -> None:
    # Many siblings, no nesting: exercises the no-paren fast split path.
    leaves = [f"sp{i:03d}:{i * 0.001:.3f}" for i in range(50)]
    newick = "(" + ",".join(leaves) + ");"
    tree = from_newick(newick)
    children = tree["Internal_0"]
    assert len(children) == 50
    assert children["sp049"] == 0.049


def test_branch_length_precision_and_missing_lengths() -> None:
    tree = from_newick("(a,b:1.5):2.25;")
    assert tree["Internal_0"] == {"a": 0.0, "b": 1.5}


def test_labels_with_commas_inside_names_are_invalid_newick_but_do_not_crash_splitter() -> None:
    # Labels cannot contain unquoted commas; the parser must raise rather
    # than silently mis-split.
    with pytest.raises(ValueError):
        from_newick("(a,b,c:0.1,);")


def test_unmatched_parenthesis_still_raises_after_optimization() -> None:
    with pytest.raises(ValueError, match="Unmatched parenthesis"):
        from_newick("((a:0.1,b:0.2):0.3;")
    with pytest.raises(ValueError, match="Unmatched parenthesis"):
        from_newick("(a:0.1")
    with pytest.raises(ValueError, match="Unmatched parenthesis"):
        from_newick("((")


def test_trailing_close_paren_tolerance_matches_original_semantics() -> None:
    # Verified against the pre-optimization parser: stray ")" produces the
    # same (permissive) result — the matching-paren scan finds the first ")"
    # identically in both implementations. Pinned here so future strictness
    # is a deliberate change, not a fast-path regression.
    tree = from_newick("(a:0.1)):")
    assert tree == {"a": None, ")": {"a": 0.1}}


def test_single_leaf_and_internal_only_trees() -> None:
    tree = from_newick("solo:0.42;")
    assert tree["solo"] is None
    tree2 = from_newick("((x,y),((p,q),r));")
    # auto-generated internal names are deterministic and unique
    assert "Internal_0" in tree2 and "Internal_1" in tree2


def test_large_random_tree_parse_is_structurally_consistent() -> None:
    import random

    rng = random.Random(20260901)
    nodes = [f"sp{i:04d}:{rng.uniform(0.01, 1.0):.4f}" for i in range(500)]
    while len(nodes) > 1:
        rng.shuffle(nodes)
        pairs = [nodes[i : i + 2] for i in range(0, len(nodes), 2)]
        nodes = [(f"({p[0]},{p[1]}):{rng.uniform(0.01, 1.0):.4f}" if len(p) == 2 else p[0]) for p in pairs]
    tree = from_newick(nodes[0] + ";")
    # every internal node dict maps names -> positive floats; leaves map None
    internals = [v for v in tree.values() if isinstance(v, dict)]
    leaves = [v for v in tree.values() if v is None]
    assert len(internals) == 499  # binary tree: n-1 internal nodes
    assert len(leaves) == 500
