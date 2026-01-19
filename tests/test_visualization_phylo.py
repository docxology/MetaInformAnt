from __future__ import annotations

import matplotlib
import pytest


def test_plot_phylo_tree_smoke():
    """Test that phylogenetic tree plotting works with basic tree structure."""
    matplotlib.use("Agg")
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    try:
        import networkx
    except ImportError:
        pytest.skip("networkx required for phylogenetic tree plotting")

    from metainformant.visualization import plot_phylo_tree

    # small 3-tip tree
    names = ["A", "B", "C"]
    matrix = DistanceMatrix(names, [[0], [0.1, 0], [0.2, 0.3, 0]])
    tree = DistanceTreeConstructor().nj(matrix)

    ax = plot_phylo_tree(tree)
    assert ax is not None
