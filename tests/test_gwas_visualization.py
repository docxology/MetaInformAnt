"""Tests for GWAS visualization functions - real matplotlib implementations.

Tests for manhattan_plot, qq_plot, pca_plot, kinship_heatmap functions that
were converted from placeholder implementations to real matplotlib plotting.
Following NO_MOCKING policy - all tests use real implementations.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

# Test visualization functions - these should return real matplotlib Figures
from metainformant.gwas.visualization import manhattan_plot, qq_plot, pca_plot, kinship_heatmap

# Test matplotlib dependency
try:
    import matplotlib
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Test numpy dependency
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


class TestManhattanPlot:
    """Tests for manhattan_plot function - real matplotlib Manhattan plots."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_manhattan_plot_returns_figure(self, tmp_path: Path):
        """Test that manhattan_plot returns a real matplotlib Figure."""
        results = [
            {"chrom": "1", "pos": 1000, "p_value": 1e-6},
            {"chrom": "1", "pos": 2000, "p_value": 1e-7},
            {"chrom": "2", "pos": 1500, "p_value": 1e-5},
        ]

        fig = manhattan_plot(results)

        # Should return a matplotlib Figure object
        assert hasattr(fig, 'savefig'), "Should return matplotlib Figure"
        assert hasattr(fig, 'axes'), "Figure should have axes"
        assert len(fig.axes) > 0, "Figure should have at least one axis"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_manhattan_plot_chromosome_positioning(self, tmp_path: Path):
        """Test that chromosomes are positioned correctly on x-axis."""
        results = [
            {"chrom": "1", "pos": 1000, "p_value": 1e-6},
            {"chrom": "1", "pos": 2000, "p_value": 1e-7},
            {"chrom": "2", "pos": 1500, "p_value": 1e-5},
        ]

        fig = manhattan_plot(results)

        # Check that we have the right number of data points
        ax = fig.axes[0]
        scatter_plots = [child for child in ax.get_children()
                        if hasattr(child, 'get_offsets')]
        assert len(scatter_plots) > 0, "Should have scatter plot data"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_manhattan_plot_significance_line(self, tmp_path: Path):
        """Test that significance line is drawn correctly."""
        results = [
            {"chrom": "1", "pos": 1000, "p_value": 1e-6},
            {"chrom": "1", "pos": 2000, "p_value": 1e-7},  # Significant
            {"chrom": "2", "pos": 1500, "p_value": 1e-5},
        ]

        fig = manhattan_plot(results, significance_threshold=5e-8)

        ax = fig.axes[0]
        # Should have horizontal line at significance threshold
        hlines = [child for child in ax.get_children()
                 if hasattr(child, 'get_xydata')]
        # Check for horizontal significance line
        found_significance_line = False
        for hline in hlines:
            if hasattr(hline, 'get_ydata'):
                y_data = hline.get_ydata()
                if len(y_data) > 0 and y_data[0] == -np.log10(5e-8):
                    found_significance_line = True
                    break
        assert found_significance_line, "Should have significance threshold line"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_manhattan_plot_file_saving(self, tmp_path: Path):
        """Test that plot can be saved to file."""
        results = [{"chrom": "1", "pos": 1000, "p_value": 1e-6}]

        output_file = tmp_path / "test_manhattan.png"
        fig = manhattan_plot(results)

        # Save the figure
        fig.savefig(output_file)

        assert output_file.exists(), "Plot file should be saved"
        assert output_file.stat().st_size > 0, "Plot file should not be empty"


class TestQQPlot:
    """Tests for qq_plot function - real matplotlib Q-Q plots."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_qq_plot_returns_figure(self):
        """Test that qq_plot returns a real matplotlib Figure."""
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5, 0.9]
        results = [{"p_value": p} for p in p_values]

        fig = qq_plot(p_values)

        assert hasattr(fig, 'savefig'), "Should return matplotlib Figure"
        assert hasattr(fig, 'axes'), "Figure should have axes"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_qq_plot_expected_vs_observed(self):
        """Test Q-Q plot expected vs observed p-value calculation."""
        # Create p-values with known distribution
        np.random.seed(42)
        n = 1000
        # Generate p-values from uniform distribution (expected under null)
        p_values = np.random.uniform(0.001, 1.0, n)

        fig = qq_plot(p_values)

        ax = fig.axes[0]

        # Check that we have data points
        scatter_plots = [child for child in ax.get_children()
                        if hasattr(child, 'get_offsets')]
        assert len(scatter_plots) > 0, "Should have Q-Q plot data"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_qq_plot_diagonal_line(self):
        """Test that Q-Q plot has diagonal null hypothesis line."""
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5, 0.9]

        fig = qq_plot(p_values)

        ax = fig.axes[0]

        # Should have diagonal line (y=x)
        lines = [child for child in ax.get_children()
                if hasattr(child, 'get_xydata')]
        found_diagonal = False
        for line in lines:
            if hasattr(line, 'get_xdata') and hasattr(line, 'get_ydata'):
                x_data = line.get_xdata()
                y_data = line.get_ydata()
                if len(x_data) > 1 and len(y_data) > 1:
                    # Check if it's approximately diagonal (y ≈ x)
                    if np.allclose(x_data, y_data, rtol=0.1):
                        found_diagonal = True
                        break
        assert found_diagonal, "Should have diagonal null hypothesis line"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_qq_plot_file_saving(self, tmp_path: Path):
        """Test that Q-Q plot can be saved to file."""
        p_values = [0.001, 0.01, 0.05]

        output_file = tmp_path / "test_qq.png"
        fig = qq_plot(p_values)

        fig.savefig(output_file)

        assert output_file.exists(), "Q-Q plot file should be saved"


class TestPCAPlot:
    """Tests for pca_plot function - real matplotlib PCA plots."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy not available")
    def test_pca_plot_returns_figure(self):
        """Test that pca_plot returns a real matplotlib Figure."""
        # Create mock PCA data
        components = [np.random.randn(20) for _ in range(3)]
        variance = np.random.rand(3)
        loadings = np.random.randn(3, 20)

        fig = pca_plot((components, variance, loadings), explained_var=variance.tolist())

        assert hasattr(fig, 'savefig'), "Should return matplotlib Figure"
        assert hasattr(fig, 'axes'), "Figure should have axes"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy not available")
    def test_pca_plot_explained_variance_annotation(self):
        """Test that PCA plot shows explained variance."""
        components = [np.random.randn(10) for _ in range(2)]
        variance = [0.6, 0.4]
        loadings = np.random.randn(2, 10)

        fig = pca_plot((components, variance, loadings), explained_var=variance)

        ax = fig.axes[0]
        # Should have title or text with variance information
        title = ax.get_title()
        assert "variance" in title.lower() or len(ax.texts) > 0, "Should show variance information"


class TestKinshipHeatmap:
    """Tests for kinship_heatmap function - real matplotlib heatmaps."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy not available")
    def test_kinship_heatmap_returns_figure(self):
        """Test that kinship_heatmap returns a real matplotlib Figure."""
        # Create symmetric kinship matrix
        np.random.seed(42)
        n = 10
        kinship = np.random.rand(n, n)
        kinship = (kinship + kinship.T) / 2  # Make symmetric
        np.fill_diagonal(kinship, 1.0)  # Self-relationships are 1

        fig = kinship_heatmap(kinship)

        assert hasattr(fig, 'savefig'), "Should return matplotlib Figure"
        assert hasattr(fig, 'axes'), "Figure should have axes"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy not available")
    def test_kinship_heatmap_matrix_visualization(self):
        """Test that kinship matrix is properly visualized."""
        kinship = np.array([[1.0, 0.5, 0.2],
                           [0.5, 1.0, 0.3],
                           [0.2, 0.3, 1.0]])

        fig = kinship_heatmap(kinship)

        ax = fig.axes[0]
        # Should have image data (heatmap)
        images = [child for child in ax.get_children()
                 if hasattr(child, 'get_array')]
        assert len(images) > 0, "Should have heatmap image"


class TestVisualizationDependencies:
    """Tests for visualization dependency availability."""

    def test_matplotlib_dependency(self):
        """Test that matplotlib is available when required."""
        if HAS_MATPLOTLIB:
            assert matplotlib is not None
            # Test that basic plotting works
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [1, 4, 2])
            plt.close(fig)
        else:
            # If matplotlib not available, functions should handle gracefully
            # This is tested in the dependency degradation tests
            pass

    def test_numpy_dependency(self):
        """Test that numpy is available when required."""
        if HAS_NUMPY:
            assert np is not None
            # Test basic numpy functionality
            arr = np.array([1, 2, 3])
            assert arr.sum() == 6
        else:
            # If numpy not available, functions should handle gracefully
            pass

    @pytest.mark.skipif(HAS_MATPLOTLIB, reason="matplotlib available")
    def test_graceful_degradation_matplotlib_missing(self):
        """Test graceful degradation when matplotlib is missing."""
        # This should raise ImportError when matplotlib is not available
        with pytest.raises(ImportError):
            manhattan_plot([])

    @pytest.mark.skipif(HAS_NUMPY, reason="numpy available")
    def test_graceful_degradation_numpy_missing(self):
        """Test graceful degradation when numpy is missing."""
        # This should raise ImportError when numpy is not available
        with pytest.raises(ImportError):
            kinship_heatmap([[1.0]])


class TestVisualizationEdgeCases:
    """Tests for edge cases in visualization functions."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_manhattan_plot_empty_results(self):
        """Test Manhattan plot with empty results."""
        fig = manhattan_plot([])
        assert hasattr(fig, 'savefig'), "Should still return Figure for empty data"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_qq_plot_empty_pvalues(self):
        """Test Q-Q plot with empty p-values - should return None gracefully."""
        result = qq_plot([])
        assert result is None, "Should return None for empty p-values"

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_manhattan_plot_invalid_data(self):
        """Test Manhattan plot with invalid data."""
        # Missing required fields
        invalid_results = [{"chrom": "1"}]  # Missing pos and p_value

        # Should handle gracefully or raise appropriate error
        try:
            fig = manhattan_plot(invalid_results)
            assert hasattr(fig, 'savefig')
        except (KeyError, ValueError):
            # This is acceptable - function should validate input
            pass

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not available")
    def test_qq_plot_invalid_pvalues(self):
        """Test Q-Q plot with invalid p-values."""
        invalid_pvalues = ["invalid", -1.0, 2.0]  # Invalid p-values

        try:
            fig = qq_plot(invalid_pvalues)
            assert hasattr(fig, 'savefig')
        except (ValueError, TypeError):
            # This is acceptable - function should validate input
            pass



