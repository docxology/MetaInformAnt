from __future__ import annotations

import matplotlib


def test_plot_phylo_tree_smoke():
    matplotlib.use("Agg")
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
    from metainformant.visualization import plot_phylo_tree

    # small 3-tip tree
    names = ["A", "B", "C"]
    matrix = DistanceMatrix(names, [[0], [0.1, 0], [0.2, 0.3, 0]])
    tree = DistanceTreeConstructor().nj(matrix)

    ax = plot_phylo_tree(tree)
    assert ax is not None


