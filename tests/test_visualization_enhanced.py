"""Enhanced tests for visualization functionality."""

import pytest
import numpy as np
from metainformant.visualization import plots


class TestVisualizationEnhanced:
    """Test enhanced visualization functionality."""

    def test_scatter_plot_functionality(self):
        """Test enhanced scatter plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test basic scatter plot
        fig, ax = plt.subplots()
        x = [1, 2, 3, 4, 5]
        y = [1, 4, 2, 3, 5]

        result_ax = plots.scatter_plot(x, y, ax=ax)
        assert result_ax is ax

        # Test with custom styling
        fig2, ax2 = plt.subplots()
        result_ax2 = plots.scatter_plot(
            x, y,
            ax=ax2,
            color='red',
            size=50,
            alpha=0.8,
            xlabel='X Axis',
            ylabel='Y Axis',
            title='Test Scatter Plot'
        )

        # Verify labels and title were set
        assert ax2.get_xlabel() == 'X Axis'
        assert ax2.get_ylabel() == 'Y Axis'
        assert ax2.get_title() == 'Test Scatter Plot'

        plt.close('all')

    def test_histogram_functionality(self):
        """Test enhanced histogram functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test basic histogram
        fig, ax = plt.subplots()
        data = [1, 2, 2, 3, 3, 3, 4, 4, 5]

        result_ax = plots.histogram(data, ax=ax)
        assert result_ax is ax

        # Test with custom options
        fig2, ax2 = plt.subplots()
        result_ax2 = plots.histogram(
            data,
            ax=ax2,
            bins=10,
            density=True,
            alpha=0.6,
            color='blue',
            xlabel='Values',
            ylabel='Frequency',
            title='Test Histogram'
        )

        # Verify settings
        assert ax2.get_xlabel() == 'Values'
        assert ax2.get_ylabel() == 'Frequency'
        assert ax2.get_title() == 'Test Histogram'

        plt.close('all')

    def test_box_plot_functionality(self):
        """Test enhanced box plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test with single dataset
        fig, ax = plt.subplots()
        data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        result_ax = plots.box_plot(data, ax=ax)
        assert result_ax is ax

        # Test with multiple datasets
        fig2, ax2 = plt.subplots()
        data_multi = [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]]

        result_ax2 = plots.box_plot(
            data_multi,
            ax=ax2,
            positions=[1, 2, 3],
            labels=['Group A', 'Group B', 'Group C'],
            xlabel='Groups',
            ylabel='Values',
            title='Test Box Plot'
        )

        # Verify labels
        assert ax2.get_xlabel() == 'Groups'
        assert ax2.get_ylabel() == 'Values'
        assert ax2.get_title() == 'Test Box Plot'

        plt.close('all')

    def test_violin_plot_functionality(self):
        """Test enhanced violin plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test with single dataset
        fig, ax = plt.subplots()
        data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        result_ax = plots.violin_plot(data, ax=ax)
        assert result_ax is ax

        # Test with multiple datasets
        fig2, ax2 = plt.subplots()
        data_multi = [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]]

        result_ax2 = plots.violin_plot(
            data_multi,
            ax=ax2,
            positions=[1, 2, 3],
            labels=['A', 'B', 'C'],
            xlabel='Categories',
            ylabel='Values',
            title='Test Violin Plot'
        )

        # Verify labels
        assert ax2.get_xlabel() == 'Categories'
        assert ax2.get_ylabel() == 'Values'
        assert ax2.get_title() == 'Test Violin Plot'

        plt.close('all')

    def test_plot_error_handling(self):
        """Test error handling in plotting functions."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test empty data
        fig, ax = plt.subplots()

        with pytest.raises((ValueError, TypeError)):
            plots.scatter_plot([], [], ax=ax)

        with pytest.raises((ValueError, TypeError)):
            plots.histogram([], ax=ax)

        plt.close('all')

    def test_plot_with_matplotlib_backend(self):
        """Test that plots work with different matplotlib backends."""
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        # Test that all plot types work with Agg backend
        x = [1, 2, 3, 4]
        y = [1, 4, 2, 3]

        plots.scatter_plot(x, y, ax=ax)
        plots.histogram(x, ax=ax)
        plots.box_plot(x, ax=ax)
        plots.violin_plot(x, ax=ax)

        plt.close('all')

        # Test with different backends if available
        try:
            matplotlib.use('SVG')
            fig, ax = plt.subplots()
            plots.scatter_plot(x, y, ax=ax)
            plt.close('all')
        except ImportError:
            pass  # SVG backend not available

    def test_plot_return_values(self):
        """Test that plot functions return the correct axes objects."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        # Test that each function returns the axes object
        result = plots.scatter_plot([1, 2, 3], [1, 2, 3], ax=ax)
        assert result is ax

        fig2, ax2 = plt.subplots()
        result = plots.histogram([1, 2, 3], ax=ax2)
        assert result is ax2

        fig3, ax3 = plt.subplots()
        result = plots.box_plot([1, 2, 3], ax=ax3)
        assert result is ax3

        fig4, ax4 = plt.subplots()
        result = plots.violin_plot([1, 2, 3], ax=ax4)
        assert result is ax4

        plt.close('all')

    def test_plot_data_types(self):
        """Test that plots work with different data types."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test with numpy arrays
        import numpy as np
        x = np.array([1, 2, 3, 4])
        y = np.array([1, 4, 2, 3])

        fig, ax = plt.subplots()
        plots.scatter_plot(x, y, ax=ax)
        plt.close(fig)

        # Test with lists
        fig2, ax2 = plt.subplots()
        plots.histogram([1, 2, 2, 3, 3, 3], ax=ax2)
        plt.close(fig2)

        # Test with tuples
        fig3, ax3 = plt.subplots()
        plots.box_plot(([1, 2, 3], [2, 3, 4]), ax=ax3)
        plt.close(fig3)

import pytest
import numpy as np
import pandas as pd
from metainformant.visualization import plots


class TestVisualizationEnhanced:
    """Test enhanced visualization functionality."""

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

import pytest
import numpy as np
import pandas as pd
from metainformant.visualization import plots


class TestVisualizationAdvanced:
    """Test advanced visualization functionality."""

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
        np.random.seed(42)
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
        np.random.seed(42)
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

    def test_enrichment_plot(self):
        """Test enrichment plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test enrichment data
        np.random.seed(42)
        data = pd.DataFrame({
            'pathway': [f'Pathway_{i}' for i in range(50)],
            'odds_ratio': np.random.exponential(1, 50),
            'p_value': 10 ** (-np.random.exponential(2, 50)),
            'log_p': -np.log10(10 ** (-np.random.exponential(2, 50)))
        })
        
        fig, ax = plt.subplots(figsize=(10, 8))
        result_ax = plots.enrichment_plot(data, 'odds_ratio', 'log_p', ax=ax)
        assert result_ax is ax
        plt.close(fig)

    def test_network_plot(self):
        """Test network plot functionality."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create test network data
        nodes = ['A', 'B', 'C', 'D', 'E']
        edges = [('A', 'B'), ('A', 'C'), ('B', 'D'), ('C', 'E'), ('D', 'E')]
        
        fig, ax = plt.subplots()
        result_ax = plots.network_plot(nodes, edges, ax=ax)
        assert result_ax is ax
        plt.close(fig)

    def test_plot_error_handling(self):
        """Test error handling in enhanced plotting functions."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test empty data for QQ plot
        fig, ax = plt.subplots()
        result_ax = plots.qq_plot([], ax=ax)
        assert result_ax is ax
        plt.close(fig)

        # Test PCA plot with missing columns
        fig, ax = plt.subplots()
        data = pd.DataFrame({'A': [1, 2, 3]})
        
        with pytest.raises(ValueError):
            plots.pca_plot(data, ax=ax)
            
        plt.close(fig)

    def test_plot_data_types(self):
        """Test that enhanced plots work with different data types."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Test with numpy arrays for correlation heatmap
        np.random.seed(42)
        data = np.random.randn(10, 5)
        df = pd.DataFrame(data, columns=['A', 'B', 'C', 'D', 'E'])
        
        fig, ax = plt.subplots()
        plots.correlation_heatmap(df, ax=ax)
        plt.close(fig)

        # Test with mixed data types for volcano plot
        data = pd.DataFrame({
            'log2FoldChange': [1.5, -2.1, 0.8, -1.2],
            'padj': [0.001, 0.01, 0.05, 0.1]
        })
        
        fig, ax = plt.subplots()
        plots.volcano_plot(data, 'log2FoldChange', 'padj', ax=ax)
        plt.close(fig)
