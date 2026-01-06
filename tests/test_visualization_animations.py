"""Tests for animated visualization functions."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest

from metainformant.visualization.animations import (
    animate_time_series,
    animate_evolution,
    animate_clustering,
    animate_network,
    animate_trajectory,
)

# Check for optional dependencies
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


class TestAnimateTimeSeries:
    """Test animate_time_series function."""

    def test_basic_time_series_animation(self):
        """Test basic time series animation creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample time series data
        data = np.random.randn(2, 50)  # 2 series, 50 time points

        fig, anim = animate_time_series(data)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_time_series_animation_single_series(self):
        """Test time series animation with single series."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.sin(np.linspace(0, 4*np.pi, 30))

        fig, anim = animate_time_series(data)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_time_series_animation_with_output_path(self, tmp_path: Path):
        """Test time series animation with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(1, 20)
        output_path = tmp_path / "timeseries.gif"

        fig, anim = animate_time_series(data, output_path=output_path, interval=50)
        assert fig is not None
        assert anim is not None
        # Note: GIF creation might not work in test environment
        plt.close('all')

    def test_time_series_animation_no_numpy(self):
        """Test time series animation when numpy is not available."""
        if hasattr(np, 'array'):  # numpy is available
            pytest.skip("numpy is available")

        data = [[1, 2, 3], [4, 5, 6]]

        with pytest.raises(ImportError, match="numpy required"):
            animate_time_series(data)


class TestAnimateEvolution:
    """Test animate_evolution function."""

    def test_basic_sequence_evolution_animation(self):
        """Test basic sequence evolution animation creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample sequences showing evolution
        sequences = [
            "ATCG",
            "ATCG",  # No change initially
            "ATCG",  # Still no change
            "ATCG"   # Final sequence
        ]

        fig, anim = animate_evolution(sequences)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_sequence_evolution_animation_with_mutations(self):
        """Test sequence evolution animation with actual mutations."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        sequences = [
            "ATCGATCG",
            "ATCGATCG",  # Generation 1
            "ATCGATCA",  # Generation 2: mutation at position 8
            "ATCGATCA",  # Generation 3: same
            "ATCGATGA"   # Generation 4: another mutation
        ]

        fig, anim = animate_evolution(sequences)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_sequence_evolution_animation_empty_sequences(self):
        """Test sequence evolution animation with empty sequences."""
        sequences = []

        with pytest.raises(ValueError, match="cannot be empty"):
            animate_evolution(sequences)


class TestAnimateClustering:
    """Test animate_clustering function."""

    def test_basic_clustering_animation(self):
        """Test basic clustering animation creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample 2D data points
        np.random.seed(42)
        data = np.random.randn(20, 2)

        # Create mock cluster assignments over time
        cluster_labels_over_time = [
            np.random.randint(0, 3, 20) for _ in range(5)  # 5 time steps
        ]

        fig, anim = animate_clustering(data, cluster_labels_over_time)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_clustering_animation_wrong_dimensions(self):
        """Test clustering animation with wrong data dimensions."""
        data = np.random.randn(10, 3)  # 3D data instead of 2D
        cluster_labels = [np.random.randint(0, 2, 10)]

        with pytest.raises(ValueError, match="must be 2D"):
            animate_clustering(data, cluster_labels)

    def test_clustering_animation_empty_labels(self):
        """Test clustering animation with empty cluster labels."""
        data = np.random.randn(10, 2)
        cluster_labels_over_time = []

        with pytest.raises(ValueError, match="cannot be empty"):
            animate_clustering(data, cluster_labels_over_time)


class TestAnimateNetwork:
    """Test animate_network function."""

    def test_basic_network_animation(self):
        """Test basic network animation creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network animations")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create graphs showing network evolution
        graphs = []
        for i in range(3):
            G = nx.erdos_renyi_graph(5 + i, 0.3, seed=42 + i)
            graphs.append(G)

        fig, anim = animate_network(graphs)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_network_animation_empty_graphs(self):
        """Test network animation with empty graphs list."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for network animations")

        graphs = []

        with pytest.raises(ValueError, match="cannot be empty"):
            animate_network(graphs)

    def test_network_animation_no_networkx(self):
        """Test network animation when NetworkX is not available."""
        if HAS_NETWORKX:
            pytest.skip("NetworkX is available")

        graphs = ["fake graph"]

        with pytest.raises(ImportError, match="NetworkX required"):
            animate_network(graphs)


class TestAnimateTrajectory:
    """Test animate_trajectory function."""

    def test_basic_trajectory_animation(self):
        """Test basic trajectory animation creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample trajectories
        np.random.seed(42)
        trajectories = [
            np.random.randn(20, 2).cumsum(axis=0) * 0.1,  # Random walk trajectory 1
            np.random.randn(25, 2).cumsum(axis=0) * 0.1   # Random walk trajectory 2
        ]

        fig, anim = animate_trajectory(trajectories)
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_trajectory_animation_single_trajectory(self):
        """Test trajectory animation with single trajectory."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        trajectory = np.array([[i, i**2] for i in range(10)], dtype=float)

        fig, anim = animate_trajectory([trajectory])
        assert fig is not None
        assert anim is not None
        plt.close('all')

    def test_trajectory_animation_wrong_dimensions(self):
        """Test trajectory animation with wrong trajectory dimensions."""
        # 3D trajectory instead of 2D
        trajectory = np.random.randn(10, 3)

        with pytest.raises(ValueError, match="must be a 2D array with shape"):
            animate_trajectory([trajectory])

    def test_trajectory_animation_empty_trajectories(self):
        """Test trajectory animation with empty trajectories list."""
        trajectories = []

        with pytest.raises(ValueError, match="cannot be empty"):
            animate_trajectory(trajectories)

    def test_trajectory_animation_no_numpy(self):
        """Test trajectory animation when numpy is not available."""
        if hasattr(np, 'array'):  # numpy is available
            pytest.skip("numpy is available")

        trajectories = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]

        with pytest.raises(ImportError, match="numpy required"):
            animate_trajectory(trajectories)






