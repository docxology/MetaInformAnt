"""Tests for multi-dimensional data visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.plots.multidim import (
    plot_3d_scatter,
    plot_pairwise_relationships,
    plot_parallel_coordinates,
    plot_radar_chart,
)

# Check for optional dependencies
try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


class TestPlotPairwiseRelationships:
    """Test plot_pairwise_relationships function."""

    def test_basic_pairwise_relationships(self):
        """Test basic pairwise relationships plot creation."""
        if not HAS_SEABORN:
            pytest.skip("seaborn required for pairwise relationships plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample multi-dimensional data
        np.random.seed(42)
        data = pd.DataFrame(
            {
                "var1": np.random.randn(50),
                "var2": np.random.randn(50),
                "var3": np.random.randn(50),
                "var4": np.random.randn(50),
            }
        )

        ax = plot_pairwise_relationships(data)
        # Note: pairplot creates its own figure, so ax might be None
        # This is expected behavior for seaborn pairplot
        plt.close("all")

    def test_pairwise_relationships_insufficient_columns(self):
        """Test pairwise relationships with insufficient numeric columns."""
        if not HAS_SEABORN:
            pytest.skip("seaborn required for pairwise relationships plotting")

        data = pd.DataFrame({"var1": [1, 2, 3], "category": ["A", "B", "C"]})  # Only one numeric column

        with pytest.raises(ValueError, match="at least 2 numeric columns"):
            plot_pairwise_relationships(data)

    def test_pairwise_relationships_no_seaborn(self):
        """Test pairwise relationships when seaborn is not available."""
        if HAS_SEABORN:
            pytest.skip("seaborn is available")

        data = pd.DataFrame({"var1": np.random.randn(10), "var2": np.random.randn(10)})

        with pytest.raises(ImportError, match="seaborn required"):
            plot_pairwise_relationships(data)


class TestPlotParallelCoordinates:
    """Test plot_parallel_coordinates function."""

    def test_basic_parallel_coordinates(self):
        """Test basic parallel coordinates plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample multi-dimensional data
        np.random.seed(42)
        data = pd.DataFrame(
            {
                "var1": np.random.randn(20),
                "var2": np.random.randn(20),
                "var3": np.random.randn(20),
                "group": ["A"] * 10 + ["B"] * 10,
            }
        )

        ax = plot_parallel_coordinates(data, color="group")
        assert ax is not None
        plt.close("all")

    def test_parallel_coordinates_no_color_column(self):
        """Test parallel coordinates without color column."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        data = pd.DataFrame({"var1": np.random.randn(10), "var2": np.random.randn(10), "var3": np.random.randn(10)})

        ax = plot_parallel_coordinates(data)
        assert ax is not None
        plt.close("all")

    def test_parallel_coordinates_insufficient_columns(self):
        """Test parallel coordinates with insufficient numeric columns."""
        data = pd.DataFrame({"var1": [1, 2, 3], "category": ["A", "B", "C"]})

        with pytest.raises(ValueError, match="at least 2 numeric columns"):
            plot_parallel_coordinates(data)

    def test_parallel_coordinates_with_output_path(self, tmp_path: Path):
        """Test parallel coordinates plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        data = pd.DataFrame({"var1": np.random.randn(10), "var2": np.random.randn(10), "var3": np.random.randn(10)})
        output_path = tmp_path / "parallel_coords.png"

        ax = plot_parallel_coordinates(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlotRadarChart:
    """Test plot_radar_chart function."""

    def test_basic_radar_chart(self):
        """Test basic radar chart creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample data for radar chart
        data = pd.DataFrame(
            {
                "metric1": [0.8, 0.6],
                "metric2": [0.7, 0.9],
                "metric3": [0.9, 0.5],
                "metric4": [0.6, 0.8],
                "metric5": [0.8, 0.7],
            },
            index=["Group A", "Group B"],
        )

        ax = plot_radar_chart(data)
        assert ax is not None
        plt.close("all")

    def test_radar_chart_insufficient_columns(self):
        """Test radar chart with insufficient columns."""
        data = pd.DataFrame({"metric1": [0.8, 0.6], "metric2": [0.7, 0.9]})

        with pytest.raises(ValueError, match="at least 3 numeric columns"):
            plot_radar_chart(data)

    def test_radar_chart_with_output_path(self, tmp_path: Path):
        """Test radar chart with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        data = pd.DataFrame({"A": [0.8, 0.6, 0.9], "B": [0.7, 0.9, 0.5], "C": [0.9, 0.5, 0.8]})
        output_path = tmp_path / "radar.png"

        ax = plot_radar_chart(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlot3DScatter:
    """Test plot_3d_scatter function."""

    def test_basic_3d_scatter(self):
        """Test basic 3D scatter plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample 3D data
        np.random.seed(42)
        data = pd.DataFrame(
            {
                "x": np.random.randn(50),
                "y": np.random.randn(50),
                "z": np.random.randn(50),
                "color_var": np.random.randn(50),
            }
        )

        ax = plot_3d_scatter(data, c="color_var")
        assert ax is not None
        plt.close("all")

    def test_3d_scatter_explicit_columns(self):
        """Test 3D scatter with explicitly specified columns."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        data = pd.DataFrame(
            {"feature1": np.random.randn(30), "feature2": np.random.randn(30), "feature3": np.random.randn(30)}
        )

        ax = plot_3d_scatter(data, x_col="feature1", y_col="feature2", z_col="feature3")
        assert ax is not None
        plt.close("all")

    def test_3d_scatter_insufficient_columns(self):
        """Test 3D scatter with insufficient numeric columns."""
        data = pd.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})

        with pytest.raises(ValueError, match="at least 3 numeric columns"):
            plot_3d_scatter(data)

    def test_3d_scatter_missing_column(self):
        """Test 3D scatter with missing specified column."""
        data = pd.DataFrame({"x": np.random.randn(10), "y": np.random.randn(10), "z": np.random.randn(10)})

        with pytest.raises(ValueError, match="not found in numeric data"):
            plot_3d_scatter(data, x_col="missing_col")

    def test_3d_scatter_with_output_path(self, tmp_path: Path):
        """Test 3D scatter plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        data = pd.DataFrame({"x": np.random.randn(20), "y": np.random.randn(20), "z": np.random.randn(20)})
        output_path = tmp_path / "3d_scatter.png"

        ax = plot_3d_scatter(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")
