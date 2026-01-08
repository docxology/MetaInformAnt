"""Tests for network visualization functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.visualization.genomics.networks import (
    plot_network_basic,
    plot_network_circular,
    plot_network_hierarchical,
    plot_network_force_directed,
    plot_community_network,
)

# Check for optional dependencies
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


class TestPlotNetworkBasic:
    """Test plot_network_basic function."""

    def test_basic_network_plot(self):
        """Test basic network plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create a simple test graph
        G = nx.erdos_renyi_graph(10, 0.3, seed=42)

        ax = plot_network_basic(G)
        assert ax is not None
        plt.close('all')

    def test_basic_network_plot_with_output_path(self, tmp_path: Path):
        """Test basic network plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.complete_graph(5)
        output_path = tmp_path / "basic_network.png"

        ax = plot_network_basic(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_basic_network_plot_no_networkx(self):
        """Test basic network plot when NetworkX is not available."""
        if HAS_NETWORKX:
            pytest.skip("NetworkX is available")

        # Mock graph object (won't work, but tests import error)
        fake_graph = "not a networkx graph"

        with pytest.raises(ImportError, match="NetworkX required"):
            plot_network_basic(fake_graph)


class TestPlotNetworkCircular:
    """Test plot_network_circular function."""

    def test_circular_network_plot(self):
        """Test circular network plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.cycle_graph(8)

        ax = plot_network_circular(G)
        assert ax is not None
        plt.close('all')

    def test_circular_network_plot_with_output_path(self, tmp_path: Path):
        """Test circular network plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.path_graph(6)
        output_path = tmp_path / "circular_network.png"

        ax = plot_network_circular(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')


class TestPlotNetworkHierarchical:
    """Test plot_network_hierarchical function."""

    def test_hierarchical_network_plot(self):
        """Test hierarchical network plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create a simple directed graph that can be hierarchical
        G = nx.DiGraph()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'), ('C', 'E')])

        ax = plot_network_hierarchical(G)
        assert ax is not None
        plt.close('all')

    def test_hierarchical_network_plot_undirected(self):
        """Test hierarchical network plot with undirected graph."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.star_graph(5)  # Undirected star graph

        ax = plot_network_hierarchical(G)
        assert ax is not None
        plt.close('all')

    def test_hierarchical_network_plot_with_output_path(self, tmp_path: Path):
        """Test hierarchical network plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.balanced_tree(2, 3)
        output_path = tmp_path / "hierarchical_network.png"

        ax = plot_network_hierarchical(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')


class TestPlotNetworkForceDirected:
    """Test plot_network_force_directed function."""

    def test_force_directed_network_plot(self):
        """Test force-directed network plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.erdos_renyi_graph(15, 0.2, seed=42)

        ax = plot_network_force_directed(G)
        assert ax is not None
        plt.close('all')

    def test_force_directed_network_plot_with_output_path(self, tmp_path: Path):
        """Test force-directed network plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.watts_strogatz_graph(10, 3, 0.1, seed=42)
        output_path = tmp_path / "force_directed_network.png"

        ax = plot_network_force_directed(G, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')


class TestPlotCommunityNetwork:
    """Test plot_community_network function."""

    def test_community_network_plot(self):
        """Test community network plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create a test graph with communities
        G = nx.erdos_renyi_graph(12, 0.3, seed=42)

        # Assign communities (simple partitioning)
        communities = {}
        nodes = list(G.nodes())
        for i, node in enumerate(nodes):
            communities[node] = i % 3  # 3 communities

        ax = plot_community_network(G, communities)
        assert ax is not None
        plt.close('all')

    def test_community_network_plot_with_output_path(self, tmp_path: Path):
        """Test community network plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        G = nx.complete_graph(6)
        communities = {0: 0, 1: 0, 2: 1, 3: 1, 4: 2, 5: 2}
        output_path = tmp_path / "community_network.png"

        ax = plot_community_network(G, communities, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_community_network_plot_empty_communities(self):
        """Test community network plot with empty communities."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network plotting")

        G = nx.path_graph(3)
        communities = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_community_network(G, communities)

    def test_community_network_plot_no_networkx(self):
        """Test community network plot when NetworkX is not available."""
        if HAS_NETWORKX:
            pytest.skip("NetworkX is available")

        fake_graph = "not a networkx graph"
        communities = {0: 0, 1: 1}

        with pytest.raises(ImportError, match="NetworkX required"):
            plot_community_network(fake_graph, communities)






