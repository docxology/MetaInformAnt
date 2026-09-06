from __future__ import annotations

import importlib.util

import matplotlib
import pytest


def test_plot_phylo_tree_smoke():
    """Test that phylogenetic tree plotting works with basic tree structure."""
    matplotlib.use("Agg")
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    if importlib.util.find_spec("networkx") is None:
        pytest.skip("networkx required for phylogenetic tree plotting")

    from metainformant.visualization.genomics.trees import plot_phylo_tree

    # small 3-tip tree
    names = ["A", "B", "C"]
    matrix = DistanceMatrix(names, [[0], [0.1, 0], [0.2, 0.3, 0]])
    tree = DistanceTreeConstructor().nj(matrix)

    ax = plot_phylo_tree(tree)
    assert ax is not None


def test_bio_phylo_tree_conversion_to_networkx():
    """Bio.Phylo trees convert to a weighted directed graph."""
    matplotlib.use("Agg")
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    if importlib.util.find_spec("networkx") is None:
        pytest.skip("networkx required for phylogenetic tree plotting")

    import networkx as nx

    from metainformant.visualization.genomics.trees import _convert_tree_to_networkx

    names = ["A", "B", "C"]
    matrix = DistanceMatrix(names, [[0], [0.1, 0], [0.2, 0.3, 0]])
    tree = DistanceTreeConstructor().nj(matrix)

    G = _convert_tree_to_networkx(tree)

    assert isinstance(G, nx.DiGraph)
    assert {"A", "B", "C"} <= set(G.nodes())
    assert G.number_of_edges() == G.number_of_nodes() - 1  # tree is acyclic
    assert all(data["weight"] > 0 for _, _, data in G.edges(data=True))
