"""Tests for dimensionality reduction visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.analysis.dimred import (
    biplot,
    plot_pca,
    plot_pca_loadings,
    plot_tsne,
    plot_umap,
)

# Check for optional dependencies
try:
    from sklearn.decomposition import PCA

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import umap

    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False


class TestPlotPCA:
    """Test plot_pca function."""

    def test_basic_pca_plot(self):
        """Test basic PCA plot creation."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(50, 10)

        ax = plot_pca(data)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close("all")

    def test_pca_plot_with_dataframe(self):
        """Test PCA plot with pandas DataFrame."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = pd.DataFrame(np.random.randn(30, 5), columns=["feat1", "feat2", "feat3", "feat4", "feat5"])

        ax = plot_pca(data)
        assert ax is not None
        plt.close("all")

    def test_pca_plot_3d(self):
        """Test 3D PCA plot."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(50, 10)

        ax = plot_pca(data, n_components=3)
        assert ax is not None
        plt.close("all")

    def test_pca_plot_with_output_path(self, tmp_path: Path):
        """Test PCA plot with output path."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(30, 8)
        output_path = tmp_path / "pca.png"

        ax = plot_pca(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_pca_plot_invalid_n_components(self):
        """Test PCA plot with invalid n_components."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA plotting")

        data = np.random.randn(20, 5)

        with pytest.raises(ValueError, match="n_components must be 2 or 3"):
            plot_pca(data, n_components=4)

    def test_pca_plot_no_sklearn(self):
        """Test PCA plot when scikit-learn is not available."""
        if HAS_SKLEARN:
            pytest.skip("scikit-learn is available")

        data = np.random.randn(20, 5)

        with pytest.raises(ImportError, match="scikit-learn required"):
            plot_pca(data)


class TestPlotUMAP:
    """Test plot_umap function."""

    def test_basic_umap_plot(self):
        """Test basic UMAP plot creation."""
        if not HAS_UMAP:
            pytest.skip("umap-learn required for UMAP plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(50, 10)

        ax = plot_umap(data)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close("all")

    def test_umap_plot_3d(self):
        """Test 3D UMAP plot."""
        if not HAS_UMAP:
            pytest.skip("umap-learn required for UMAP plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(50, 10)

        ax = plot_umap(data, n_components=3)
        assert ax is not None
        plt.close("all")

    def test_umap_plot_no_umap(self):
        """Test UMAP plot when umap-learn is not available."""
        if HAS_UMAP:
            pytest.skip("umap-learn is available")

        data = np.random.randn(20, 5)

        with pytest.raises(ImportError, match="umap-learn required"):
            plot_umap(data)


class TestPlotTSNE:
    """Test plot_tsne function."""

    def test_basic_tsne_plot(self):
        """Test basic t-SNE plot creation."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for t-SNE plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(30, 10)  # Smaller dataset for speed

        ax = plot_tsne(data)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close("all")

    def test_tsne_plot_3d(self):
        """Test 3D t-SNE plot."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for t-SNE plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(30, 10)

        ax = plot_tsne(data, n_components=3)
        assert ax is not None
        plt.close("all")

    def test_tsne_plot_no_sklearn(self):
        """Test t-SNE plot when scikit-learn is not available."""
        if HAS_SKLEARN:
            pytest.skip("scikit-learn is available")

        data = np.random.randn(20, 5)

        with pytest.raises(ImportError, match="scikit-learn required"):
            plot_tsne(data)


class TestPlotPCALoadings:
    """Test plot_pca_loadings function."""

    def test_basic_pca_loadings_plot(self):
        """Test basic PCA loadings plot creation."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA loadings plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create and fit PCA model
        np.random.seed(42)
        data = np.random.randn(50, 10)
        pca = PCA(n_components=2)
        pca.fit(data)

        ax = plot_pca_loadings(pca)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close("all")

    def test_pca_loadings_plot_with_output_path(self, tmp_path: Path):
        """Test PCA loadings plot with output path."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA loadings plotting")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(40, 8)
        pca = PCA(n_components=2)
        pca.fit(data)
        output_path = tmp_path / "pca_loadings.png"

        ax = plot_pca_loadings(pca, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_pca_loadings_plot_invalid_n_components(self):
        """Test PCA loadings plot with invalid n_components."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA loadings plotting")

        np.random.seed(42)
        data = np.random.randn(30, 5)
        pca = PCA(n_components=3)
        pca.fit(data)

        with pytest.raises(ValueError, match="only supports 2 components"):
            plot_pca_loadings(pca, n_components=3)

    def test_pca_loadings_plot_invalid_model(self):
        """Test PCA loadings plot with invalid model."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for PCA loadings plotting")

        # Pass non-PCA object
        invalid_model = "not a pca model"

        with pytest.raises(ValueError, match="must be a fitted sklearn PCA model"):
            plot_pca_loadings(invalid_model)


class TestBiplot:
    """Test biplot function."""

    def test_basic_biplot(self):
        """Test basic biplot creation."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for biplot")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(50, 6)
        pca = PCA(n_components=2)
        pca.fit(data)

        ax = biplot(data, pca)
        assert ax is not None
        # Biplot should have both scatter points and arrows
        plt.close("all")

    def test_biplot_with_dataframe(self):
        """Test biplot with pandas DataFrame."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for biplot")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = pd.DataFrame(np.random.randn(30, 4), columns=["gene1", "gene2", "gene3", "gene4"])
        pca = PCA(n_components=2)
        pca.fit(data.values)

        ax = biplot(data, pca)
        assert ax is not None
        plt.close("all")

    def test_biplot_with_output_path(self, tmp_path: Path):
        """Test biplot with output path."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for biplot")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        np.random.seed(42)
        data = np.random.randn(40, 5)
        pca = PCA(n_components=2)
        pca.fit(data)
        output_path = tmp_path / "biplot.png"

        ax = biplot(data, pca, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_biplot_invalid_model(self):
        """Test biplot with invalid PCA model."""
        if not HAS_SKLEARN:
            pytest.skip("scikit-learn required for biplot")

        data = np.random.randn(20, 3)
        invalid_model = "not a pca model"

        with pytest.raises(ValueError, match="must be a fitted sklearn PCA model"):
            biplot(data, invalid_model)

    def test_biplot_no_sklearn(self):
        """Test biplot when scikit-learn is not available."""
        if HAS_SKLEARN:
            pytest.skip("scikit-learn is available")

        data = np.random.randn(20, 3)
        fake_model = type("MockPCA", (), {"components_": np.random.randn(3, 2)})()

        with pytest.raises(ImportError, match="scikit-learn required"):
            biplot(data, fake_model)
