"""Tests for information theory visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.information import (
    plot_entropy_profile,
    plot_mutual_information_matrix,
    plot_renyi_spectra,
    plot_information_landscape,
    plot_information_network,
)

# Check for optional dependencies
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


class TestPlotEntropyProfile:
    """Test plot_entropy_profile function."""

    def test_basic_entropy_profile_plot(self):
        """Test basic entropy profile plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample entropy data
        entropy_data = {
            'Sequence_A': [1.5, 2.1, 2.8, 3.2, 3.5],
            'Sequence_B': [1.3, 1.9, 2.4, 2.9, 3.1]
        }

        ax = plot_entropy_profile(entropy_data)
        assert ax is not None
        assert len(ax.lines) == 2  # Two sequences
        plt.close('all')

    def test_entropy_profile_with_output_path(self, tmp_path: Path):
        """Test entropy profile plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        entropy_data = {
            'Seq1': [1.0, 1.5, 2.0, 2.3],
            'Seq2': [0.8, 1.2, 1.8, 2.1]
        }
        output_path = tmp_path / "entropy_profile.png"

        ax = plot_entropy_profile(entropy_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_entropy_profile_empty_data(self):
        """Test entropy profile plot with empty data."""
        entropy_data = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_entropy_profile(entropy_data)

    def test_entropy_profile_invalid_data(self):
        """Test entropy profile plot with invalid data types."""
        entropy_data = {
            'Seq1': "not a list"  # Should be list/array
        }

        with pytest.raises(ValueError, match="must be a list or array"):
            plot_entropy_profile(entropy_data)


class TestPlotMutualInformationMatrix:
    """Test plot_mutual_information_matrix function."""

    def test_basic_mutual_information_matrix_plot(self):
        """Test basic mutual information matrix plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample MI matrix
        np.random.seed(42)
        n_vars = 5
        mi_matrix = np.random.rand(n_vars, n_vars)
        # Make it symmetric
        mi_matrix = (mi_matrix + mi_matrix.T) / 2
        # Set diagonal to some value (self-information)
        np.fill_diagonal(mi_matrix, 1.0)

        labels = [f'Var_{i}' for i in range(n_vars)]

        ax = plot_mutual_information_matrix(mi_matrix, labels)
        assert ax is not None
        plt.close('all')

    def test_mutual_information_matrix_without_labels(self):
        """Test mutual information matrix plot without labels."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        mi_matrix = np.array([[1.0, 0.5, 0.3],
                             [0.5, 1.0, 0.2],
                             [0.3, 0.2, 1.0]])

        ax = plot_mutual_information_matrix(mi_matrix)
        assert ax is not None
        plt.close('all')

    def test_mutual_information_matrix_with_output_path(self, tmp_path: Path):
        """Test mutual information matrix plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        mi_matrix = np.array([[1.0, 0.8], [0.8, 1.0]])
        output_path = tmp_path / "mi_matrix.png"

        ax = plot_mutual_information_matrix(mi_matrix, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_mutual_information_matrix_nonsquare(self):
        """Test mutual information matrix plot with non-square matrix."""
        mi_matrix = np.array([[1.0, 0.5, 0.3],
                             [0.5, 1.0, 0.2]])  # Not square

        with pytest.raises(ValueError, match="must be square"):
            plot_mutual_information_matrix(mi_matrix)

    def test_mutual_information_matrix_mismatched_labels(self):
        """Test mutual information matrix plot with mismatched label count."""
        mi_matrix = np.array([[1.0, 0.5], [0.5, 1.0]])
        labels = ['A', 'B', 'C']  # Too many labels

        with pytest.raises(ValueError, match="must match matrix dimensions"):
            plot_mutual_information_matrix(mi_matrix, labels)


class TestPlotRenyiSpectra:
    """Test plot_renyi_spectra function."""

    def test_basic_renyi_spectra_plot(self):
        """Test basic Rényi spectra plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample Rényi data
        alpha_values = [0.5, 1.0, 1.5, 2.0]
        renyi_data = {
            'Distribution_A': [2.1, 2.0, 1.9, 1.8],
            'Distribution_B': [2.3, 2.2, 2.0, 1.9]
        }

        ax = plot_renyi_spectra(renyi_data, alpha_values)
        assert ax is not None
        assert len(ax.lines) == 2  # Two distributions
        plt.close('all')

    def test_renyi_spectra_with_output_path(self, tmp_path: Path):
        """Test Rényi spectra plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        alpha_values = [0.5, 1.0, 2.0]
        renyi_data = {
            'Dist1': [1.8, 1.5, 1.3],
            'Dist2': [2.0, 1.8, 1.5]
        }
        output_path = tmp_path / "renyi_spectra.png"

        ax = plot_renyi_spectra(renyi_data, alpha_values, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_renyi_spectra_empty_data(self):
        """Test Rényi spectra plot with empty data."""
        alpha_values = [0.5, 1.0, 2.0]
        renyi_data = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_renyi_spectra(renyi_data, alpha_values)

    def test_renyi_spectra_mismatched_lengths(self):
        """Test Rényi spectra plot with mismatched data lengths."""
        alpha_values = [0.5, 1.0, 2.0]  # 3 values
        renyi_data = {
            'Dist1': [1.8, 1.5]  # Only 2 values
        }

        with pytest.raises(ValueError, match="must match alpha values length"):
            plot_renyi_spectra(renyi_data, alpha_values)


class TestPlotInformationLandscape:
    """Test plot_information_landscape function."""

    def test_basic_information_landscape_plot(self):
        """Test basic information landscape plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample landscape data
        np.random.seed(42)
        landscape_data = np.random.rand(10, 10)

        ax = plot_information_landscape(landscape_data)
        assert ax is not None
        plt.close('all')

    def test_information_landscape_with_coords(self):
        """Test information landscape plot with coordinate arrays."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        landscape_data = np.random.rand(5, 5)
        x_coords = np.linspace(0, 10, 5)
        y_coords = np.linspace(0, 10, 5)

        ax = plot_information_landscape(landscape_data, x_coords, y_coords)
        assert ax is not None
        plt.close('all')

    def test_information_landscape_with_output_path(self, tmp_path: Path):
        """Test information landscape plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        landscape_data = np.random.rand(4, 4)
        output_path = tmp_path / "info_landscape.png"

        ax = plot_information_landscape(landscape_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_information_landscape_invalid_dimensions(self):
        """Test information landscape plot with invalid dimensions."""
        landscape_data = np.array([1, 2, 3])  # 1D array

        with pytest.raises(ValueError, match="must be 2D"):
            plot_information_landscape(landscape_data)


class TestPlotInformationNetwork:
    """Test plot_information_network function."""

    def test_basic_information_network_plot(self):
        """Test basic information network plot creation."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for information network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample network data
        nodes = ['A', 'B', 'C', 'D', 'E']
        edges = [('A', 'B', 0.8), ('A', 'C', 0.6), ('B', 'D', 0.7),
                ('C', 'D', 0.5), ('D', 'E', 0.9)]

        ax = plot_information_network(nodes, edges)
        assert ax is not None
        plt.close('all')

    def test_information_network_with_output_path(self, tmp_path: Path):
        """Test information network plot with output path."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for information network plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        nodes = ['X', 'Y', 'Z']
        edges = [('X', 'Y', 0.5), ('Y', 'Z', 0.7)]
        output_path = tmp_path / "info_network.png"

        ax = plot_information_network(nodes, edges, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_information_network_empty_nodes(self):
        """Test information network plot with empty nodes."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for information network plotting")

        nodes = []
        edges = [('A', 'B', 0.5)]

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_information_network(nodes, edges)

    def test_information_network_unknown_node(self):
        """Test information network plot with unknown node in edge."""
        if not HAS_NETWORKX:
            pytest.skip("NetworkX required for information network plotting")

        nodes = ['A', 'B', 'C']
        edges = [('A', 'D', 0.5)]  # D is not in nodes

        with pytest.raises(ValueError, match="connects unknown nodes"):
            plot_information_network(nodes, edges)

    def test_information_network_no_networkx(self):
        """Test information network plot when NetworkX is not available."""
        if HAS_NETWORKX:
            pytest.skip("NetworkX is available")

        nodes = ['A', 'B']
        edges = [('A', 'B', 0.5)]

        with pytest.raises(ImportError, match="NetworkX required"):
            plot_information_network(nodes, edges)





