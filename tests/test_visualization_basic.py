"""Tests for basic visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.visualization.basic import (
    lineplot,
    scatter_plot,
    heatmap,
    bar_plot,
    pie_chart,
    area_plot,
    step_plot,
)


class TestLineplot:
    """Test lineplot function."""

    def test_basic_lineplot(self, tmp_path: Path):
        """Test basic line plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])

        ax = lineplot(x, y)
        assert ax is not None
        assert len(ax.lines) == 1
        plt.close('all')

    def test_lineplot_single_array(self, tmp_path: Path):
        """Test line plot with single array (y only)."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        y = np.array([1, 4, 2, 8, 5])

        ax = lineplot(y)
        assert ax is not None
        assert len(ax.lines) == 1
        assert len(ax.lines[0].get_xdata()) == len(y)
        assert len(ax.lines[0].get_ydata()) == len(y)
        plt.close('all')

    def test_lineplot_with_output_path(self, tmp_path: Path):
        """Test line plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])
        output_path = tmp_path / "lineplot.png"

        ax = lineplot(x, y, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_lineplot_mismatched_lengths(self):
        """Test line plot with mismatched array lengths."""
        x = np.array([1, 2, 3])
        y = np.array([1, 4, 2, 8])

        with pytest.raises(ValueError, match="same length"):
            lineplot(x, y)


class TestScatterPlot:
    """Test scatter_plot function."""

    def test_basic_scatter_plot(self):
        """Test basic scatter plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])

        ax = scatter_plot(x, y)
        assert ax is not None
        assert len(ax.collections) == 1  # scatter plots use collections
        plt.close('all')

    def test_scatter_plot_with_output_path(self, tmp_path: Path):
        """Test scatter plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])
        output_path = tmp_path / "scatter.png"

        ax = scatter_plot(x, y, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_scatter_plot_mismatched_lengths(self):
        """Test scatter plot with mismatched array lengths."""
        x = np.array([1, 2, 3])
        y = np.array([1, 4, 2, 8])

        with pytest.raises(ValueError, match="same length"):
            scatter_plot(x, y)


class TestHeatmap:
    """Test heatmap function."""

    def test_basic_heatmap(self):
        """Test basic heatmap creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.rand(10, 10)

        ax = heatmap(data)
        assert ax is not None
        assert len(ax.images) == 1  # heatmaps use images
        plt.close('all')

    def test_heatmap_with_custom_cmap(self):
        """Test heatmap with custom colormap."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.rand(5, 5)

        ax = heatmap(data, cmap='plasma')
        assert ax is not None
        assert ax.images[0].get_cmap().name == 'plasma'
        plt.close('all')

    def test_heatmap_with_output_path(self, tmp_path: Path):
        """Test heatmap with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.rand(5, 5)
        output_path = tmp_path / "heatmap.png"

        ax = heatmap(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_heatmap_invalid_dimensions(self):
        """Test heatmap with non-2D data."""
        data = np.array([1, 2, 3])  # 1D array

        with pytest.raises(ValueError, match="2D array"):
            heatmap(data)


class TestBarPlot:
    """Test bar_plot function."""

    def test_basic_bar_plot(self):
        """Test basic bar plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        height = np.array([1, 4, 2, 8, 5])

        ax = bar_plot(x, height)
        assert ax is not None
        assert len(ax.patches) == 5  # bar plots use patches
        plt.close('all')

    def test_bar_plot_with_strings(self):
        """Test bar plot with string x-values."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = ['A', 'B', 'C', 'D', 'E']
        height = np.array([1, 4, 2, 8, 5])

        ax = bar_plot(x, height)
        assert ax is not None
        assert len(ax.patches) == 5
        plt.close('all')

    def test_bar_plot_with_output_path(self, tmp_path: Path):
        """Test bar plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        height = np.array([1, 4, 2, 8, 5])
        output_path = tmp_path / "barplot.png"

        ax = bar_plot(x, height, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_bar_plot_mismatched_lengths(self):
        """Test bar plot with mismatched array lengths."""
        x = np.array([1, 2, 3])
        height = np.array([1, 4, 2, 8])

        with pytest.raises(ValueError, match="same length"):
            bar_plot(x, height)


class TestPieChart:
    """Test pie_chart function."""

    def test_basic_pie_chart(self):
        """Test basic pie chart creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        sizes = np.array([15, 30, 45, 10])

        ax = pie_chart(sizes)
        assert ax is not None
        plt.close('all')

    def test_pie_chart_with_labels(self):
        """Test pie chart with labels."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        sizes = np.array([15, 30, 45, 10])
        labels = ['A', 'B', 'C', 'D']

        ax = pie_chart(sizes, labels=labels)
        assert ax is not None
        plt.close('all')

    def test_pie_chart_with_output_path(self, tmp_path: Path):
        """Test pie chart with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        sizes = np.array([15, 30, 45, 10])
        output_path = tmp_path / "piechart.png"

        ax = pie_chart(sizes, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_pie_chart_mismatched_labels(self):
        """Test pie chart with mismatched label count."""
        sizes = np.array([15, 30, 45, 10])
        labels = ['A', 'B', 'C']  # Too few labels

        with pytest.raises(ValueError, match="same length"):
            pie_chart(sizes, labels=labels)


class TestAreaPlot:
    """Test area_plot function."""

    def test_basic_area_plot(self):
        """Test basic area plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])

        ax = area_plot(x, y)
        assert ax is not None
        assert len(ax.collections) >= 1  # area plots use fill_between which creates collections
        plt.close('all')

    def test_area_plot_with_output_path(self, tmp_path: Path):
        """Test area plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])
        output_path = tmp_path / "areaplot.png"

        ax = area_plot(x, y, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_area_plot_mismatched_lengths(self):
        """Test area plot with mismatched array lengths."""
        x = np.array([1, 2, 3])
        y = np.array([1, 4, 2, 8])

        with pytest.raises(ValueError, match="same length"):
            area_plot(x, y)


class TestStepPlot:
    """Test step_plot function."""

    def test_basic_step_plot(self):
        """Test basic step plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])

        ax = step_plot(x, y)
        assert ax is not None
        assert len(ax.lines) == 1
        plt.close('all')

    def test_step_plot_with_output_path(self, tmp_path: Path):
        """Test step plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 2, 8, 5])
        output_path = tmp_path / "stepplot.png"

        ax = step_plot(x, y, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_step_plot_mismatched_lengths(self):
        """Test step plot with mismatched array lengths."""
        x = np.array([1, 2, 3])
        y = np.array([1, 4, 2, 8])

        with pytest.raises(ValueError, match="same length"):
            step_plot(x, y)






