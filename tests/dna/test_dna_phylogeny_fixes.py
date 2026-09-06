"""Regression tests for phylogeny fixes (NJ selection, UPGMA heights, bootstrap)."""

import random

import pytest

from metainformant.dna.alignment.distances import p_distance
from metainformant.dna.phylogeny.tree_analysis import (
    _extract_clades,
    _find_root,
    _get_all_leaves,
    _get_bipartitions,
    bootstrap_support,
    is_monophyletic,
    to_ascii,
)
from metainformant.dna.phylogeny.tree_construction import neighbor_joining_tree, upgma_tree


def _evolve_four_taxa() -> dict:
    """Deterministically evolve 4 taxa on ((A:0.2,B:0.04):0.04,(C:0.02,D:0.36):0.04).

    The raw closest pair in the resulting p-distance matrix is the NON-cherry
    pair (B, C): correct neighbor joining must recover the {A,B} and {C,D}
    cherries anyway (Q-criterion selection).
    """
    rng = random.Random(42)
    length = 2400
    bases = "ACGT"
    edges_from = {
        "root": [("N1", 0.04), ("N2", 0.04)],
        "N1": [("A", 0.2), ("B", 0.04)],
        "N2": [("C", 0.02), ("D", 0.36)],
    }
    sequences: dict = {}

    def evolve(node: str, seq: list) -> None:
        for child, branch_len in edges_from[node]:
            child_seq = list(seq)
            for i in range(length):
                if rng.random() < branch_len:
                    child_seq[i] = rng.choice([b for b in bases if b != child_seq[i]])
            if child in "ABCD":
                sequences[child] = "".join(child_seq)
            else:
                evolve(child, child_seq)

    evolve("root", [rng.choice(bases) for _ in range(length)])
    return sequences


class TestNeighborJoining:
    def test_recovers_true_topology_despite_deceptive_raw_minimum(self):
        seqs = _evolve_four_taxa()
        distances = {a: {b: p_distance(seqs[a], seqs[b]) for b in seqs} for a in seqs}
        pairs = [(a, b) for i, a in enumerate("ABCD") for b in "ABCD"[i + 1 :]]
        raw_min = min(pairs, key=lambda p: distances[p[0]][p[1]])
        # The setup must be genuinely discriminating: the raw closest pair is
        # NOT a cherry of the true tree.
        assert set(raw_min) == {"B", "C"}

        tree = neighbor_joining_tree(seqs)
        leaves = set(_get_all_leaves(tree))
        splits = set(_get_bipartitions(tree, leaves))
        # The unrooted tree must contain the AB|CD split (either side may be
        # reported as the smaller side).
        assert frozenset({"A", "B"}) in splits or frozenset({"C", "D"}) in splits


class TestUpgma:
    def test_leaf_depths_are_ultrametric(self):
        # Depth-3 merge: {A,B} merge first, then join with C. Every leaf must
        # end up equidistant from the root.
        seqs = {"A": "ACGTACGTAC", "B": "ACGTACGTAT", "C": "TTGATTGATT"}
        tree = upgma_tree(seqs)
        root = _find_root(tree)

        depths: dict = {}

        def walk(node: str, depth: float) -> None:
            data = tree[node]
            if isinstance(data, dict):
                for child, branch_len in data.items():
                    if child == "bootstrap":
                        continue
                    walk(child, depth + branch_len)
            else:
                depths[node] = depth

        walk(root, 0.0)
        assert len(depths) == 3
        assert max(depths.values()) - min(depths.values()) < 1e-9

    def test_nested_cluster_branches_extend_to_merge_height(self):
        seqs = {"A": "ACGTACGTAC", "B": "ACGTACGTAT", "C": "TTGATTGATT"}
        tree = upgma_tree(seqs)
        root = _find_root(tree)
        # The AB cluster merges at d(A,B)/2 = 0.05; the root merges at the
        # averaged distance / 2. A and B must sit at the same depth as C.
        depths: dict = {}

        def walk(node: str, depth: float) -> None:
            data = tree[node]
            if isinstance(data, dict):
                for child, branch_len in data.items():
                    if child == "bootstrap":
                        continue
                    walk(child, depth + branch_len)
            else:
                depths[node] = depth

        walk(root, 0.0)
        assert depths["A"] == pytest.approx(depths["C"])
        assert depths["B"] == pytest.approx(depths["C"])


class TestBootstrapSupport:
    def _setup(self):
        seqs = {
            "A": "AAAAAAAACCCCCCCC",
            "B": "AAAAAAAACCCCCCCC",
            "C": "GGGGGGGGTTTTTTTT",
            "D": "GGGGGGGGTTTTTTTT",
        }
        tree = neighbor_joining_tree(seqs)
        return tree, seqs

    def test_extract_clades_returns_internal_clades(self):
        tree, _ = self._setup()
        clades = _extract_clades(tree)
        assert clades, "_extract_clades must traverse from the tree root"
        assert any(set(clade) == {"A", "B"} for clade in clades)

    def test_internal_nodes_receive_bootstrap_values(self):
        tree, seqs = self._setup()
        supported = bootstrap_support(tree, seqs, n_replicates=12)
        root = _find_root(supported)
        values: dict = {}

        def collect(node: str) -> None:
            data = supported[node]
            if isinstance(data, dict):
                if "bootstrap" in data:
                    values[node] = data["bootstrap"]
                for child in data:
                    if child != "bootstrap":
                        collect(child)

        collect(root)
        assert values, "every internal node must be visited and annotated"
        # Both clusters are internally identical: perfect support expected.
        assert all(v == 100 for v in values.values())

    def test_support_is_deterministic_across_calls(self):
        tree, seqs = self._setup()
        first = bootstrap_support(tree, seqs, n_replicates=12)
        second = bootstrap_support(tree, seqs, n_replicates=12)
        assert first == second

    def test_to_ascii_does_not_print_bootstrap_marker(self):
        tree, seqs = self._setup()
        supported = bootstrap_support(tree, seqs, n_replicates=12)
        art = to_ascii(supported)
        assert "bootstrap" not in art
        for leaf in ("A", "B", "C", "D"):
            assert leaf in art

    def test_monophyletic_clusters(self):
        tree, _ = self._setup()
        assert is_monophyletic(tree, ["A", "B"])
