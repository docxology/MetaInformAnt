"""Tests for statistical visualization functions."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.statistical import (
    histogram,
    box_plot,
    violin_plot,
    qq_plot,
    correlation_heatmap,
    density_plot,
    ridge_plot,
    roc_curve,
    precision_recall_curve,
    residual_plot,
    leverage_plot,
)

# Check for optional dependencies
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

try:
    from sklearn.linear_model import LinearRegression
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


class TestHistogram:
    """Test histogram function."""

    def test_basic_histogram(self):
        """Test basic histogram creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)

        ax = histogram(data)
        assert ax is not None
        assert len(ax.patches) > 0  # histogram creates patches
        plt.close('all')

    def test_histogram_with_custom_bins(self):
        """Test histogram with custom bin count."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)

        ax = histogram(data, bins=50)
        assert ax is not None
        plt.close('all')

    def test_histogram_with_output_path(self, tmp_path: Path):
        """Test histogram with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)
        output_path = tmp_path / "histogram.png"

        ax = histogram(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_histogram_empty_data(self):
        """Test histogram with empty data."""
        data = np.array([])

        with pytest.raises(ValueError, match="cannot be empty"):
            histogram(data)


class TestBoxPlot:
    """Test box_plot function."""

    def test_basic_box_plot(self):
        """Test basic box plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = [np.random.randn(50) for _ in range(3)]

        ax = box_plot(data)
        assert ax is not None
        assert len(ax.lines) > 0  # boxplot creates lines
        plt.close('all')

    def test_box_plot_with_output_path(self, tmp_path: Path):
        """Test box plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = [np.random.randn(50) for _ in range(3)]
        output_path = tmp_path / "boxplot.png"

        ax = box_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_box_plot_empty_data(self):
        """Test box plot with empty data list."""
        data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            box_plot(data)

    def test_box_plot_empty_array(self):
        """Test box plot with empty array in data."""
        data = [np.array([]), np.random.randn(10)]

        with pytest.raises(ValueError, match="cannot be empty"):
            box_plot(data)


class TestViolinPlot:
    """Test violin_plot function."""

    def test_basic_violin_plot(self):
        """Test basic violin plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = [np.random.randn(50) for _ in range(3)]

        ax = violin_plot(data)
        assert ax is not None
        plt.close('all')

    def test_violin_plot_with_output_path(self, tmp_path: Path):
        """Test violin plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = [np.random.randn(50) for _ in range(3)]
        output_path = tmp_path / "violinplot.png"

        ax = violin_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_violin_plot_empty_data(self):
        """Test violin plot with empty data list."""
        data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            violin_plot(data)


class TestQQPlot:
    """Test qq_plot function."""

    def test_basic_qq_plot(self):
        """Test basic Q-Q plot creation."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for Q-Q plot")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)

        ax = qq_plot(data)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter plot creates collections
        plt.close('all')

    def test_qq_plot_uniform_distribution(self):
        """Test Q-Q plot with uniform distribution."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for Q-Q plot")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.uniform(0, 1, 100)

        ax = qq_plot(data, distribution="uniform")
        assert ax is not None
        plt.close('all')

    def test_qq_plot_with_output_path(self, tmp_path: Path):
        """Test Q-Q plot with output path."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for Q-Q plot")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)
        output_path = tmp_path / "qqplot.png"

        ax = qq_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_qq_plot_empty_data(self):
        """Test Q-Q plot with empty data."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for Q-Q plot")

        data = np.array([])

        with pytest.raises(ValueError, match="cannot be empty"):
            qq_plot(data)

    def test_qq_plot_unsupported_distribution(self):
        """Test Q-Q plot with unsupported distribution."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for Q-Q plot")

        data = np.random.randn(100)

        with pytest.raises(ValueError, match="Unsupported distribution"):
            qq_plot(data, distribution="invalid")

    def test_qq_plot_no_scipy(self):
        """Test Q-Q plot when scipy is not available."""
        if HAS_SCIPY:
            pytest.skip("scipy is available")

        data = np.random.randn(100)

        with pytest.raises(ImportError, match="scipy required"):
            qq_plot(data)


class TestCorrelationHeatmap:
    """Test correlation_heatmap function."""

    def test_basic_correlation_heatmap(self):
        """Test basic correlation heatmap creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'A': np.random.randn(50),
            'B': np.random.randn(50),
            'C': np.random.randn(50),
            'D': np.random.randn(50)
        })

        ax = correlation_heatmap(data)
        assert ax is not None
        plt.close('all')

    def test_correlation_heatmap_with_output_path(self, tmp_path: Path):
        """Test correlation heatmap with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'A': np.random.randn(50),
            'B': np.random.randn(50)
        })
        output_path = tmp_path / "corr_heatmap.png"

        ax = correlation_heatmap(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_correlation_heatmap_no_numeric_columns(self):
        """Test correlation heatmap with no numeric columns."""
        data = pd.DataFrame({
            'A': ['a', 'b', 'c'],
            'B': ['x', 'y', 'z']
        })

        with pytest.raises(ValueError, match="must contain numeric columns"):
            correlation_heatmap(data)


class TestDensityPlot:
    """Test density_plot function."""

    def test_basic_density_plot(self):
        """Test basic density plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)

        ax = density_plot(data)
        assert ax is not None
        plt.close('all')

    def test_density_plot_with_output_path(self, tmp_path: Path):
        """Test density plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = np.random.randn(100)
        output_path = tmp_path / "density.png"

        ax = density_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_density_plot_empty_data(self):
        """Test density plot with empty data."""
        data = np.array([])

        with pytest.raises(ValueError, match="cannot be empty"):
            density_plot(data)


class TestRidgePlot:
    """Test ridge_plot function."""

    def test_basic_ridge_plot(self):
        """Test basic ridge plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = [np.random.randn(50) + i for i in range(3)]

        ax = ridge_plot(data)
        assert ax is not None
        plt.close('all')

    def test_ridge_plot_with_output_path(self, tmp_path: Path):
        """Test ridge plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = [np.random.randn(50) for _ in range(3)]
        output_path = tmp_path / "ridge.png"

        ax = ridge_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_ridge_plot_empty_data(self):
        """Test ridge plot with empty data list."""
        data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            ridge_plot(data)


class TestROCCurve:
    """Test roc_curve function."""

    def test_basic_roc_curve(self):
        """Test basic ROC curve creation."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for ROC curve")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        y_true = np.random.randint(0, 2, 100)
        y_scores = np.random.rand(100)

        ax = roc_curve(y_true, y_scores)
        assert ax is not None
        assert len(ax.lines) > 0  # plot creates lines
        plt.close('all')

    def test_roc_curve_with_output_path(self, tmp_path: Path):
        """Test ROC curve with output path."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for ROC curve")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        y_true = np.random.randint(0, 2, 100)
        y_scores = np.random.rand(100)
        output_path = tmp_path / "roc.png"

        ax = roc_curve(y_true, y_scores, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_roc_curve_mismatched_lengths(self):
        """Test ROC curve with mismatched array lengths."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for ROC curve")

        y_true = np.array([0, 1, 0])
        y_scores = np.array([0.1, 0.5])

        with pytest.raises(ValueError, match="same length"):
            roc_curve(y_true, y_scores)

    def test_roc_curve_no_scipy(self):
        """Test ROC curve when scipy is not available."""
        if HAS_SCIPY:
            pytest.skip("scipy is available")

        y_true = np.array([0, 1, 0])
        y_scores = np.array([0.1, 0.5, 0.3])

        with pytest.raises(ImportError, match="scipy required"):
            roc_curve(y_true, y_scores)


class TestPrecisionRecallCurve:
    """Test precision_recall_curve function."""

    def test_basic_precision_recall_curve(self):
        """Test basic precision-recall curve creation."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for precision-recall curve")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        y_true = np.random.randint(0, 2, 100)
        y_scores = np.random.rand(100)

        ax = precision_recall_curve(y_true, y_scores)
        assert ax is not None
        assert len(ax.lines) > 0
        plt.close('all')

    def test_precision_recall_curve_mismatched_lengths(self):
        """Test precision-recall curve with mismatched array lengths."""
        if not HAS_SCIPY:
            pytest.skip("scipy required for precision-recall curve")

        y_true = np.array([0, 1, 0])
        y_scores = np.array([0.1, 0.5])

        with pytest.raises(ValueError, match="same length"):
            precision_recall_curve(y_true, y_scores)


class TestResidualPlot:
    """Test residual_plot function."""

    def test_basic_residual_plot(self):
        """Test basic residual plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        y_true = np.random.randn(50)
        y_pred = y_true + np.random.randn(50) * 0.1

        ax = residual_plot(y_true, y_pred)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close('all')

    def test_residual_plot_with_output_path(self, tmp_path: Path):
        """Test residual plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        y_true = np.random.randn(50)
        y_pred = y_true + np.random.randn(50) * 0.1
        output_path = tmp_path / "residual.png"

        ax = residual_plot(y_true, y_pred, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_residual_plot_mismatched_lengths(self):
        """Test residual plot with mismatched array lengths."""
        y_true = np.array([1, 2, 3])
        y_pred = np.array([1, 2])

        with pytest.raises(ValueError, match="same length"):
            residual_plot(y_true, y_pred)


class TestLeveragePlot:
    """Test leverage_plot function."""

    def test_basic_leverage_plot(self):
        """Test basic leverage plot creation."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for leverage plot")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        X = np.random.randn(50, 3)
        y = X @ np.array([1, 2, 3]) + np.random.randn(50) * 0.1

        ax = leverage_plot(X, y)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close('all')

    def test_leverage_plot_mismatched_lengths(self):
        """Test leverage plot with mismatched array dimensions."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for leverage plot")

        X = np.random.randn(50, 3)
        y = np.random.randn(30)

        with pytest.raises(ValueError, match="same number of samples"):
            leverage_plot(X, y)

    def test_leverage_plot_no_sklearn(self):
        """Test leverage plot when scikit-learn is not available."""
        if HAS_SKLEARN:
            pytest.skip("scikit-learn is available")

        X = np.random.randn(50, 3)
        y = np.random.randn(50)

        with pytest.raises(ImportError, match="scikit-learn required"):
            leverage_plot(X, y)
