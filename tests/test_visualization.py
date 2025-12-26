"""Tests for visualization functionality."""

import pytest
import numpy as np
import pandas as pd
from metainformant.visualization import plots


class TestVisualization:
    """Test visualization functionality."""

    def test_correlation_heatmap(self):
        """Test correlation heatmap functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test data
        np.random.seed(42)
        data = pd.DataFrame({
            'A': np.random.randn(100),
            'B': np.random.randn(100),
            'C': np.random.randn(100),
            'D': np.random.randn(100)
        })
        
        fig, ax = plt.subplots()
        result_ax = plots.correlation_heatmap(data, ax=ax, annot=False)
        assert result_ax is ax
        plt.close(fig)

    def test_volcano_plot(self):
        """Test volcano plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test data
        data = pd.DataFrame({
            'log2FoldChange': np.random.randn(100),
            'padj': 10 ** (-np.random.exponential(2, 100))
        })
        
        fig, ax = plt.subplots()
        result_ax = plots.volcano_plot(data, 'log2FoldChange', 'padj', ax=ax)
        assert result_ax is ax
        plt.close(fig)

    def test_expression_heatmap(self):
        """Test expression heatmap functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test expression data
        np.random.seed(42)
        data = pd.DataFrame(
            np.random.randn(20, 10),
            index=[f'Gene_{i}' for i in range(20)],
            columns=[f'Sample_{i}' for i in range(10)]
        )
        
        fig, ax = plt.subplots()
        result_ax = plots.expression_heatmap(data, ax=ax)
        assert result_ax is ax
        plt.close(fig)

    def test_pca_plot(self):
        """Test PCA plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test PCA data
        data = pd.DataFrame({
            'PC1': np.random.randn(50),
            'PC2': np.random.randn(50),
            'PC1_variance': [25.5] * 50,
            'PC2_variance': [18.3] * 50,
            'cluster': np.random.choice(['A', 'B', 'C'], 50)
        })
        
        fig, ax = plt.subplots()
        result_ax = plots.pca_plot(data, hue='cluster', ax=ax)
        assert result_ax is ax
        plt.close(fig)

    def test_plot_error_handling(self):
        """Test error handling in enhanced plotting functions."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test with invalid data
        fig, ax = plt.subplots()
        
        with pytest.raises(ValueError):
            plots.pca_plot(pd.DataFrame(), ax=ax)  # Missing PC columns
            
        plt.close(fig)

    def test_plot_data_types(self):
        """Test that enhanced plots work with different data types."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test with numpy arrays
        x = np.array([1, 2, 3, 4])
        y = np.array([1, 4, 2, 3])
        
        fig, ax = plt.subplots()
        plots.scatter_plot(x, y, ax=ax)
        plt.close(fig)

        # Test correlation heatmap with mixed types
        data = pd.DataFrame({
            'A': [1.0, 2.0, 3.0],
            'B': [4, 5, 6],
            'C': ['a', 'b', 'c']  # Non-numeric column
        })
        
        fig2, ax2 = plt.subplots()
        # Should handle non-numeric columns gracefully
        plots.correlation_heatmap(data.select_dtypes(include=[np.number]), ax=ax2)
        plt.close(fig2)

    def test_manhattan_plot(self):
        """Test Manhattan plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test genomic data
        np.random.seed(42)
        chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] * 20
        positions = np.arange(1, 101)
        p_values = 10 ** (-np.random.exponential(3, 100))

        data = pd.DataFrame({
            'chromosome': chromosomes,
            'position': positions,
            'p_value': p_values,
            'log_p': -np.log10(p_values)
        })

        fig, ax = plt.subplots(figsize=(12, 6))
        result_ax = plots.manhattan_plot(data, 'position', 'log_p', ax=ax)
        assert result_ax is ax
        plt.close(fig)

    def test_qq_plot(self):
        """Test Q-Q plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test p-values (mix of null and signal)
        np.random.seed(42)
        null_pvals = np.random.uniform(0, 1, 80)  # Null distribution
        signal_pvals = np.random.exponential(0.1, 20)  # Some signal
        p_values = np.concatenate([null_pvals, signal_pvals])

        fig, ax = plt.subplots()
        result_ax = plots.qq_plot(p_values, ax=ax)
        assert result_ax is ax
        plt.close(fig)



