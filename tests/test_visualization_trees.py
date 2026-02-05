"""Tests for phylogenetic tree visualization functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.visualization.genomics.trees import (
    circular_tree_plot,
    plot_phylo_tree,
    tree_annotation_plot,
    tree_comparison_plot,
    unrooted_tree_plot,
)

# Check for optional dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


class TestPlotPhyloTree:
    """Test plot_phylo_tree function."""

    def test_basic_phylogenetic_tree_plot(self):
        """Test basic phylogenetic tree plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create a simple tree
        G = nx.DiGraph()
        G.add_edges_from([("root", "A"), ("root", "B"), ("A", "A1"), ("A", "A2"), ("B", "B1")])

        ax = plot_phylo_tree(G)
        assert ax is not None
        plt.close("all")

    def test_phylogenetic_tree_plot_with_output_path(self, tmp_path: Path):
        """Test phylogenetic tree plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        G.add_edges_from([("root", "species1"), ("root", "species2")])
        output_path = tmp_path / "phylo_tree.png"

        ax = plot_phylo_tree(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_phylogenetic_tree_plot_no_networkx(self):
        """Test phylogenetic tree plot when NetworkX is not available."""
        if HAS_NETWORKX:
            pytest.skip("NetworkX is available")

        fake_tree = "not a tree"

        with pytest.raises(ImportError, match="NetworkX required"):
            plot_phylo_tree(fake_tree)


class TestCircularTreePlot:
    """Test circular_tree_plot function."""

    def test_basic_circular_tree_plot(self):
        """Test basic circular tree plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        G.add_edges_from([("root", "A"), ("root", "B"), ("root", "C")])

        ax = circular_tree_plot(G)
        assert ax is not None
        plt.close("all")

    def test_circular_tree_plot_with_output_path(self, tmp_path: Path):
        """Test circular tree plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        G.add_edges_from([("root", "leaf1"), ("root", "leaf2")])
        output_path = tmp_path / "circular_tree.png"

        ax = circular_tree_plot(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestUnrootedTreePlot:
    """Test unrooted_tree_plot function."""

    def test_basic_unrooted_tree_plot(self):
        """Test basic unrooted tree plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.Graph()  # Undirected graph for unrooted tree
        G.add_edges_from([("A", "B"), ("B", "C"), ("C", "D"), ("C", "E")])

        ax = unrooted_tree_plot(G)
        assert ax is not None
        plt.close("all")

    def test_unrooted_tree_plot_with_output_path(self, tmp_path: Path):
        """Test unrooted tree plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.Graph()
        G.add_edges_from([("X", "Y"), ("Y", "Z")])
        output_path = tmp_path / "unrooted_tree.png"

        ax = unrooted_tree_plot(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestTreeComparisonPlot:
    """Test tree_comparison_plot function."""

    def test_basic_tree_comparison_plot(self):
        """Test basic tree comparison plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create two simple trees
        G1 = nx.DiGraph()
        G1.add_edges_from([("root", "A"), ("root", "B")])

        G2 = nx.DiGraph()
        G2.add_edges_from([("root", "C"), ("root", "D"), ("C", "E")])

        ax = tree_comparison_plot(G1, G2)
        assert ax is not None
        plt.close("all")

    def test_tree_comparison_plot_with_output_path(self, tmp_path: Path):
        """Test tree comparison plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G1 = nx.DiGraph()
        G1.add_edges_from([("r1", "s1")])

        G2 = nx.DiGraph()
        G2.add_edges_from([("r2", "s2")])

        output_path = tmp_path / "tree_comparison.png"

        ax = tree_comparison_plot(G1, G2, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestTreeAnnotationPlot:
    """Test tree_annotation_plot function."""

    def test_basic_tree_annotation_plot(self):
        """Test basic tree annotation plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        G.add_edges_from([("root", "A"), ("root", "B")])

        annotations = {"A": {"color": "red", "label": "Species A"}, "B": {"color": "blue", "label": "Species B"}}

        ax = tree_annotation_plot(G, annotations)
        assert ax is not None
        plt.close("all")

    def test_tree_annotation_plot_with_output_path(self, tmp_path: Path):
        """Test tree annotation plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        G.add_edges_from([("root", "leaf")])

        annotations = {"leaf": {"color": "green", "label": "Annotated leaf"}}

        output_path = tmp_path / "annotated_tree.png"

        ax = tree_annotation_plot(G, annotations, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_tree_annotation_plot_empty_annotations(self):
        """Test tree annotation plot with empty annotations."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for phylogenetic tree plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        G.add_edges_from([("root", "A")])

        ax = tree_annotation_plot(G, {})
        assert ax is not None
        plt.close("all")
