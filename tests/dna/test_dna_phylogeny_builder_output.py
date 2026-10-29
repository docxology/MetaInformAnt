"""``dna.phylogeny`` consumers checked against what the builders really return.

The builders (``neighbor_joining_tree``, ``upgma_tree``, ``bootstrap_support``) return a flat dictionary: every node is
a key, a leaf maps to ``None`` and an internal node maps to ``{child: branch_length}`` (plus a ``"bootstrap"`` entry
once support has been added). Run on that output, four consumers were wrong:

* ``to_ascii`` printed every leaf twice (a leaf line from its parent, then another from its own recursion);
* ``bootstrap_support`` built UPGMA replicates for any ``method`` other than ``"nj"`` (``"NJ"`` or a typo), and
  counted a neighbor-joining branch only if a replicate had the same *rooted* clade although neighbor joining leaves
  the root arbitrary, so the two children of the root, one branch, read very different support (98 and 33);
* ``neighbor_joining_tree`` and ``nj_tree_from_kmer`` joined the last two nodes in set order, which follows the
  per-process hash seed, so the same sequences gave a different Newick string on different runs.

Real trees from real sequences, no test doubles. The sequences are built so their p-distances are those of a known
tree (every branch mutates private sites).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pytest

import metainformant
from metainformant.dna.phylogeny import (
    bootstrap_support,
    neighbor_joining_tree,
    nj_tree_from_kmer,
    to_ascii,
    to_newick,
    upgma_tree,
)
from metainformant.dna.phylogeny.tree_analysis import _find_root


def _sequences(branches: Dict[str, Tuple[str, int]], leaves: List[str], length: int) -> Dict[str, str]:
    """Sequences whose p-distances are the tree metric of ``branches`` (child -> (parent, private sites)) / length."""
    private: Dict[str, List[int]] = {}
    site = 0
    for child, (_, count) in branches.items():
        private[child] = list(range(site, site + count))
        site += count
    assert site <= length

    def sequence(node: str) -> List[str]:
        if node == "root":
            return ["A"] * length
        letters = sequence(branches[node][0])
        for index in private[node]:
            letters[index] = "T"
        return letters

    return {leaf: "".join(sequence(leaf)) for leaf in leaves}


# Six taxa: ((A,B),(C,D)) and (E,F) under the root. Unrooted, that is three cherries round one centre, so the root
# branch is a single split {E,F} | {A,B,C,D} whichever way a replicate is rooted.
SIX_LEAVES = ["A", "B", "C", "D", "E", "F"]
SIX_BRANCHES = {
    "n1": ("n3", 100),
    "n2": ("n3", 100),
    "n3": ("root", 100),
    "n4": ("root", 100),
    "A": ("n1", 60),
    "B": ("n1", 60),
    "C": ("n2", 60),
    "D": ("n2", 60),
    "E": ("n4", 60),
    "F": ("n4", 60),
}
SIX_SEQUENCES = _sequences(SIX_BRANCHES, SIX_LEAVES, 800)

FIVE_SEQUENCES = {
    "A": "ACGTACGTACGTACGTACGT",
    "B": "ACGTACGTACGAACGTACGT",
    "C": "ACGTACCTACGAACGTTCGT",
    "D": "TCGTACCTACGAACGTTCGA",
    "E": "TCGAACCTTCGAACGTTCGA",
}


def _clade(tree: Dict[str, Any], node: str) -> frozenset:
    entry = tree[node]
    if not isinstance(entry, dict):
        return frozenset([node])
    return frozenset().union(*(_clade(tree, child) for child in entry if child != "bootstrap"))


# ---------------------------------------------------------------------------
# to_ascii
# ---------------------------------------------------------------------------


def test_to_ascii_prints_a_hand_built_tree_once_per_node() -> None:
    tree = {
        "A": None,
        "B": None,
        "C": None,
        "n1": {"A": 1.0, "B": 1.0},
        "root": {"n1": 1.0, "C": 2.0},
    }

    assert to_ascii(tree) == "\n".join(["root", "├── n1", "│   ├── A", "│   └── B", "└── C"])


@pytest.mark.parametrize("builder", [neighbor_joining_tree, upgma_tree])
def test_to_ascii_lists_every_node_of_a_builder_tree_exactly_once(builder: Any) -> None:
    tree = builder(FIVE_SEQUENCES)

    lines = to_ascii(tree).splitlines()

    names = [line.lstrip("│├└─ ") for line in lines]
    assert sorted(names) == sorted(tree)  # one line per node: no leaf repeated
    assert lines[0] == _find_root(tree)


def test_to_ascii_of_a_bootstrapped_tree_is_the_text_of_the_plain_tree() -> None:
    tree = neighbor_joining_tree(FIVE_SEQUENCES)
    supported = bootstrap_support(tree, FIVE_SEQUENCES, n_replicates=10)

    assert to_ascii(supported) == to_ascii(tree)
    assert "bootstrap" not in to_ascii(supported)


# ---------------------------------------------------------------------------
# bootstrap_support
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["NJ", "UPGMA", "ml", "", "neighbor_joining"])
def test_bootstrap_support_refuses_an_unknown_method(method: str) -> None:
    """Anything but "nj" used to build UPGMA replicates without a word."""
    tree = neighbor_joining_tree(FIVE_SEQUENCES)

    with pytest.raises(ValueError, match="nj"):
        bootstrap_support(tree, FIVE_SEQUENCES, n_replicates=5, method=method)


def test_bootstrap_support_refuses_an_unknown_method_even_for_empty_sequences() -> None:
    tree = neighbor_joining_tree(FIVE_SEQUENCES)

    with pytest.raises(ValueError, match="upgma"):
        bootstrap_support(tree, {}, method="bogus")


def test_bootstrap_support_refuses_sequences_for_other_taxa() -> None:
    """Names that are not the tree's leaves matched no clade, so every branch silently read 0."""
    tree = neighbor_joining_tree(FIVE_SEQUENCES)
    renamed = {f"{name}x": sequence for name, sequence in FIVE_SEQUENCES.items()}

    with pytest.raises(ValueError, match="same taxa"):
        bootstrap_support(tree, renamed, n_replicates=5)


def test_the_two_children_of_a_neighbor_joining_root_have_the_same_support() -> None:
    """Both children of the root lie on one branch (one split), so a replicate cannot support one and not the other."""
    tree = neighbor_joining_tree(SIX_SEQUENCES)
    root = _find_root(tree)
    first, second = [child for child in tree[root] if child != "bootstrap"]
    assert _clade(tree, first) | _clade(tree, second) == frozenset(SIX_LEAVES)

    supported = bootstrap_support(tree, SIX_SEQUENCES, n_replicates=60, method="nj")

    support = {node: supported[node]["bootstrap"] for node, entry in supported.items() if isinstance(entry, dict)}
    assert support[first] == support[second]
    assert support[first] >= 90  # the signal is strong: private sites on every branch
    assert support[root] == 100


def test_neighbor_joining_support_does_not_depend_on_where_replicates_are_rooted() -> None:
    """Three cherries round one centre: every cherry split is present in (almost) every replicate."""
    tree = neighbor_joining_tree(SIX_SEQUENCES)

    supported = bootstrap_support(tree, SIX_SEQUENCES, n_replicates=60, method="nj")

    internal = {node: supported[node]["bootstrap"] for node, entry in supported.items() if isinstance(entry, dict)}
    assert len(internal) == 5  # three cherries, the centre and the root
    assert all(value >= 90 for value in internal.values()), internal  # a rooted-clade count gave 33 for one of them


def test_neighbor_joining_branch_to_a_single_leaf_is_always_supported() -> None:
    """A root child that is every leaf but one is that leaf's terminal branch, which every tree has."""
    tree = neighbor_joining_tree(FIVE_SEQUENCES)
    root = _find_root(tree)
    sizes = {child: len(_clade(tree, child)) for child in tree[root] if child != "bootstrap"}
    (single,) = [child for child, size in sizes.items() if size == len(FIVE_SEQUENCES) - 1]  # the premise of the test

    supported = bootstrap_support(tree, FIVE_SEQUENCES, n_replicates=20, method="nj")

    assert supported[single]["bootstrap"] == 100


def test_upgma_support_still_counts_rooted_clades() -> None:
    tree = upgma_tree(SIX_SEQUENCES)

    supported = bootstrap_support(tree, SIX_SEQUENCES, n_replicates=40, method="upgma")

    root = _find_root(supported)
    assert supported[root]["bootstrap"] == 100
    for node, entry in supported.items():
        if isinstance(entry, dict):
            assert 0 <= entry["bootstrap"] <= 100


# ---------------------------------------------------------------------------
# Neighbor joining is reproducible across processes
# ---------------------------------------------------------------------------


def _newick_under_hash_seed(builder: str, hash_seed: int) -> str:
    """Build a tree in a fresh interpreter with a given string-hash seed; the child imports the same package."""
    package_root = str(Path(metainformant.__file__).resolve().parent.parent)
    program = (
        "from metainformant.dna.phylogeny import neighbor_joining_tree, nj_tree_from_kmer, to_newick\n"
        "seqs = {'A': 'ACGTACGTACGTACGTACGT', 'B': 'ACGTACGTACGAACGTACGT', 'C': 'ACGTACCTACGAACGTTCGT',"
        " 'D': 'TCGTACCTACGAACGTTCGA', 'E': 'TCGAACCTTCGAACGTTCGA'}\n"
        "tree = neighbor_joining_tree(seqs) if %r == 'nj' else nj_tree_from_kmer(seqs, k=2)\n"
        "print(to_newick(tree))\n" % builder
    )
    env = {"PYTHONHASHSEED": str(hash_seed), "PYTHONPATH": package_root, "PATH": "/usr/bin:/bin"}
    result = subprocess.run([sys.executable, "-c", program], env=env, capture_output=True, text=True, check=True)
    return result.stdout.strip()


@pytest.mark.parametrize("builder", ["nj", "kmer"])
def test_the_newick_of_a_neighbor_joining_tree_does_not_depend_on_the_hash_seed(builder: str) -> None:
    outputs = {_newick_under_hash_seed(builder, seed) for seed in range(6)}

    assert len(outputs) == 1, outputs


def test_neighbor_joining_root_children_are_in_sorted_order() -> None:
    for builder in (neighbor_joining_tree, nj_tree_from_kmer):
        tree = builder(FIVE_SEQUENCES)
        root_children = [child for child in tree[_find_root(tree)] if child != "bootstrap"]

        assert root_children == sorted(root_children)
        assert to_newick(tree).count("(") == len(FIVE_SEQUENCES) - 1
