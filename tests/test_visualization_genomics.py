"""Tests for genomic visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.genomics import (
    manhattan_plot,
    volcano_plot,
    regional_plot,
    circular_manhattan_plot,
    chromosome_ideogram,
    coverage_plot,
    variant_plot,
)


class TestManhattanPlot:
    """Test manhattan_plot function."""

    def test_basic_manhattan_plot(self):
        """Test basic Manhattan plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample GWAS data
        np.random.seed(42)
        data = pd.DataFrame({
            'CHR': np.random.choice([1, 2, 3], 100),
            'BP': np.random.randint(1000000, 50000000, 100),
            'P': np.random.uniform(1e-10, 1, 100)
        })

        ax = manhattan_plot(data)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close('all')

    def test_manhattan_plot_with_custom_columns(self):
        """Test Manhattan plot with custom column names."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'chromosome': [1, 1, 2],
            'position': [1000, 2000, 1500],
            'pvalue': [1e-5, 1e-6, 1e-4]
        })

        ax = manhattan_plot(data, chr_col='chromosome', pos_col='position', pval_col='pvalue')
        assert ax is not None
        plt.close('all')

    def test_manhattan_plot_with_output_path(self, tmp_path: Path):
        """Test Manhattan plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'CHR': [1, 2, 3],
            'BP': [1000, 2000, 3000],
            'P': [0.01, 0.001, 0.05]
        })
        output_path = tmp_path / "manhattan.png"

        ax = manhattan_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_manhattan_plot_missing_columns(self):
        """Test Manhattan plot with missing required columns."""
        data = pd.DataFrame({
            'CHR': [1, 2],
            'BP': [1000, 2000]
            # Missing P column
        })

        with pytest.raises(ValueError, match="Missing required columns"):
            manhattan_plot(data)


class TestVolcanoPlot:
    """Test volcano_plot function."""

    def test_basic_volcano_plot(self):
        """Test basic volcano plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample DE analysis data
        np.random.seed(42)
        data = pd.DataFrame({
            'log2FoldChange': np.random.randn(100),
            'padj': np.random.uniform(1e-10, 1, 100)
        })

        ax = volcano_plot(data)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close('all')

    def test_volcano_plot_with_custom_columns(self):
        """Test volcano plot with custom column names."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'logFC': np.random.randn(10),
            'p_adj': np.random.uniform(1e-5, 1, 10)
        })

        ax = volcano_plot(data, log2fc_col='logFC', pval_col='p_adj')
        assert ax is not None
        plt.close('all')

    def test_volcano_plot_with_output_path(self, tmp_path: Path):
        """Test volcano plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'log2FoldChange': [1, -2, 0.5],
            'padj': [0.001, 0.01, 0.05]
        })
        output_path = tmp_path / "volcano.png"

        ax = volcano_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_volcano_plot_missing_columns(self):
        """Test volcano plot with missing required columns."""
        data = pd.DataFrame({
            'log2FoldChange': [1, 2, 3]
            # Missing padj column
        })

        with pytest.raises(ValueError, match="Missing required columns"):
            volcano_plot(data)


class TestRegionalPlot:
    """Test regional_plot function."""

    def test_basic_regional_plot(self):
        """Test basic regional plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample regional data
        data = pd.DataFrame({
            'CHR': ['1'] * 10,
            'BP': list(range(1000000, 1010000, 100000)),
            'P': np.random.uniform(1e-8, 1, 10)
        })

        ax = regional_plot(data, chr='1', start=900000, end=1100000)
        assert ax is not None
        plt.close('all')

    def test_regional_plot_empty_region(self):
        """Test regional plot with no data in region."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'CHR': ['2'] * 5,
            'BP': [2000000, 2001000, 2002000, 2003000, 2004000],
            'P': [0.1, 0.2, 0.3, 0.4, 0.5]
        })

        # Query region with no data
        ax = regional_plot(data, chr='1', start=1000000, end=2000000)
        assert ax is not None
        plt.close('all')

    def test_regional_plot_with_output_path(self, tmp_path: Path):
        """Test regional plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'CHR': ['1'] * 3,
            'BP': [1000000, 1010000, 1020000],
            'P': [0.01, 0.001, 0.05]
        })
        output_path = tmp_path / "regional.png"

        ax = regional_plot(data, chr='1', start=900000, end=1100000, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_regional_plot_invalid_region(self):
        """Test regional plot with invalid region parameters."""
        data = pd.DataFrame({
            'CHR': ['1'] * 3,
            'BP': [1000000, 1010000, 1020000],
            'P': [0.01, 0.001, 0.05]
        })

        with pytest.raises(ValueError, match="Start position must be less than end position"):
            regional_plot(data, chr='1', start=2000000, end=1000000)


class TestCircularManhattanPlot:
    """Test circular_manhattan_plot function."""

    def test_basic_circular_manhattan_plot(self):
        """Test basic circular Manhattan plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample GWAS data
        data = pd.DataFrame({
            'CHR': np.random.choice([1, 2, 3], 50),
            'BP': np.random.randint(1000000, 50000000, 50),
            'P': np.random.uniform(1e-8, 1, 50)
        })

        ax = circular_manhattan_plot(data)
        assert ax is not None
        plt.close('all')

    def test_circular_manhattan_plot_with_output_path(self, tmp_path: Path):
        """Test circular Manhattan plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.DataFrame({
            'CHR': [1, 2, 1],
            'BP': [1000000, 2000000, 3000000],
            'P': [0.01, 0.001, 0.05]
        })
        output_path = tmp_path / "circular_manhattan.png"

        ax = circular_manhattan_plot(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')


class TestChromosomeIdeogram:
    """Test chromosome_ideogram function."""

    def test_basic_chromosome_ideogram(self):
        """Test basic chromosome ideogram creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        ax = chromosome_ideogram()
        assert ax is not None
        assert len(ax.patches) > 0  # barh creates patches
        plt.close('all')

    def test_chromosome_ideogram_with_output_path(self, tmp_path: Path):
        """Test chromosome ideogram with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        output_path = tmp_path / "ideogram.png"

        ax = chromosome_ideogram(output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')


class TestCoveragePlot:
    """Test coverage_plot function."""

    def test_basic_coverage_plot(self):
        """Test basic coverage plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        positions = np.arange(1000000, 1010000, 1000)
        coverage = np.random.poisson(30, len(positions))

        ax = coverage_plot(coverage, positions)
        assert ax is not None
        assert len(ax.lines) > 0  # plot creates lines
        plt.close('all')

    def test_coverage_plot_with_output_path(self, tmp_path: Path):
        """Test coverage plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        positions = np.array([1000000, 1001000, 1002000])
        coverage = np.array([10, 20, 15])
        output_path = tmp_path / "coverage.png"

        ax = coverage_plot(coverage, positions, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_coverage_plot_mismatched_lengths(self):
        """Test coverage plot with mismatched array lengths."""
        positions = np.array([1000000, 1001000])
        coverage = np.array([10, 20, 30])

        with pytest.raises(ValueError, match="same length"):
            coverage_plot(coverage, positions)


class TestVariantPlot:
    """Test variant_plot function."""

    def test_basic_variant_plot(self):
        """Test basic variant plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        variants = pd.DataFrame({
            'POS': [1000000, 1005000, 1010000, 1015000, 1020000]
        })

        ax = variant_plot(variants)
        assert ax is not None
        assert len(ax.collections) > 0  # scatter creates collections
        plt.close('all')

    def test_variant_plot_missing_pos_column(self):
        """Test variant plot with missing POS column."""
        variants = pd.DataFrame({
            'CHR': ['1', '1', '1'],
            'REF': ['A', 'T', 'G']
        })

        with pytest.raises(ValueError, match="must contain 'POS' column"):
            variant_plot(variants)

    def test_variant_plot_with_output_path(self, tmp_path: Path):
        """Test variant plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        variants = pd.DataFrame({
            'POS': [1000000, 1010000, 1020000]
        })
        output_path = tmp_path / "variants.png"

        ax = variant_plot(variants, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')
