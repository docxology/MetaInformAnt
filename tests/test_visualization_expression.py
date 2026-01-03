"""Tests for expression analysis visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.expression import (
    plot_expression_heatmap,
    plot_enrichment_barplot,
    plot_differential_expression,
)


class TestPlotExpressionHeatmap:
    """Test plot_expression_heatmap function."""

    def test_basic_expression_heatmap(self):
        """Test basic expression heatmap creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample expression data
        np.random.seed(42)
        data = pd.DataFrame(
            np.random.randn(20, 10),
            index=[f'Gene_{i}' for i in range(20)],
            columns=[f'Sample_{i}' for i in range(10)]
        )

        ax = plot_expression_heatmap(data)
        assert ax is not None
        plt.close('all')

    def test_expression_heatmap_with_output_path(self, tmp_path: Path):
        """Test expression heatmap with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame(
            np.random.randn(5, 5),
            index=['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
            columns=['S1', 'S2', 'S3', 'S4', 'S5']
        )
        output_path = tmp_path / "expression_heatmap.png"

        ax = plot_expression_heatmap(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_expression_heatmap_empty_data(self):
        """Test expression heatmap with empty data."""
        data = pd.DataFrame()

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_expression_heatmap(data)

    def test_expression_heatmap_no_numeric_data(self):
        """Test expression heatmap with no numeric columns."""
        data = pd.DataFrame({
            'Gene': ['Gene1', 'Gene2'],
            'Sample': ['S1', 'S2']
        })

        with pytest.raises(ValueError, match="must contain numeric columns"):
            plot_expression_heatmap(data)


class TestPlotEnrichmentBarplot:
    """Test plot_enrichment_barplot function."""

    def test_basic_enrichment_barplot(self):
        """Test basic enrichment barplot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample enrichment results
        enrichment_data = pd.DataFrame({
            'Term': ['Pathway A', 'Pathway B', 'Pathway C', 'Pathway D', 'Pathway E'],
            'padj': [1e-10, 1e-8, 1e-6, 1e-4, 0.01],
            'Count': [15, 12, 8, 6, 4]
        })

        ax = plot_enrichment_barplot(enrichment_data)
        assert ax is not None
        plt.close('all')

    def test_enrichment_barplot_with_custom_columns(self):
        """Test enrichment barplot with custom column names."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        enrichment_data = pd.DataFrame({
            'pathway': ['Path A', 'Path B', 'Path C'],
            'p.adjust': [1e-5, 1e-4, 1e-3],
            'gene_count': [10, 8, 5]
        })

        ax = plot_enrichment_barplot(enrichment_data)
        assert ax is not None
        plt.close('all')

    def test_enrichment_barplot_with_output_path(self, tmp_path: Path):
        """Test enrichment barplot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        enrichment_data = pd.DataFrame({
            'Term': ['Pathway 1', 'Pathway 2'],
            'padj': [1e-5, 1e-4],
            'Count': [8, 6]
        })
        output_path = tmp_path / "enrichment_barplot.png"

        ax = plot_enrichment_barplot(enrichment_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_enrichment_barplot_empty_data(self):
        """Test enrichment barplot with empty data."""
        enrichment_data = pd.DataFrame()

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_enrichment_barplot(enrichment_data)

    def test_enrichment_barplot_missing_term_column(self):
        """Test enrichment barplot with missing term column."""
        enrichment_data = pd.DataFrame({
            'padj': [1e-5, 1e-4],
            'Count': [8, 6]
        })

        with pytest.raises(ValueError, match="must contain a term/pathway column"):
            plot_enrichment_barplot(enrichment_data)

    def test_enrichment_barplot_missing_pval_column(self):
        """Test enrichment barplot with missing p-value column."""
        enrichment_data = pd.DataFrame({
            'Term': ['Path A', 'Path B'],
            'Count': [8, 6]
        })

        with pytest.raises(ValueError, match="must contain an adjusted p-value column"):
            plot_enrichment_barplot(enrichment_data)


class TestPlotDifferentialExpression:
    """Test plot_differential_expression function."""

    def test_basic_differential_expression_plot(self):
        """Test basic differential expression plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample DE results
        np.random.seed(42)
        de_data = pd.DataFrame({
            'log2FoldChange': np.random.randn(100),
            'padj': np.random.uniform(1e-10, 1, 100),
            'baseMean': np.random.exponential(100, 100)
        })

        ax = plot_differential_expression(de_data)
        assert ax is not None
        plt.close('all')

    def test_differential_expression_plot_with_custom_columns(self):
        """Test differential expression plot with custom column names."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        de_data = pd.DataFrame({
            'logFC': np.random.randn(20),
            'adj.P.Val': np.random.uniform(1e-8, 1, 20),
            'AveExpr': np.random.exponential(50, 20)
        })

        ax = plot_differential_expression(de_data)
        assert ax is not None
        plt.close('all')

    def test_differential_expression_plot_with_output_path(self, tmp_path: Path):
        """Test differential expression plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        de_data = pd.DataFrame({
            'log2FoldChange': [1, -2, 0.5, -0.8],
            'padj': [1e-6, 1e-7, 0.01, 0.05],
            'baseMean': [100, 200, 50, 75]
        })
        output_path = tmp_path / "de_plot.png"

        ax = plot_differential_expression(de_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_differential_expression_plot_empty_data(self):
        """Test differential expression plot with empty data."""
        de_data = pd.DataFrame()

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_differential_expression(de_data)

    def test_differential_expression_plot_missing_logfc_column(self):
        """Test differential expression plot with missing log fold change column."""
        de_data = pd.DataFrame({
            'padj': [1e-5, 1e-4],
            'baseMean': [100, 200]
        })

        with pytest.raises(ValueError, match="must contain a log fold change column"):
            plot_differential_expression(de_data)

    def test_differential_expression_plot_missing_pval_column(self):
        """Test differential expression plot with missing p-value column."""
        de_data = pd.DataFrame({
            'log2FoldChange': [1, -2],
            'baseMean': [100, 200]
        })

        with pytest.raises(ValueError, match="must contain an adjusted p-value column"):
            plot_differential_expression(de_data)





