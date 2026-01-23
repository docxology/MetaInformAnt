"""Tests for quality control visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.analysis.quality import (
    plot_quality_metrics,
    plot_adapter_content,
    plot_gc_distribution,
    plot_length_distribution,
)


class TestPlotQualityMetrics:
    """Test plot_quality_metrics function."""

    def test_basic_quality_metrics_plot(self):
        """Test basic quality metrics plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample QC data
        qc_data = {
            "per_base_quality": {
                "positions": list(range(100)),
                "mean_qualities": np.random.uniform(20, 40, 100),
                "read_counts_at_position": np.random.randint(1000, 2000, 100),
            },
            "gc_content_distribution": {
                "bins": np.linspace(0, 100, 21),
                "counts": np.random.poisson(50, 20),
                "mean_gc_content": 45.0,
                "median_gc_content": 44.5,
            },
            "sequence_length_distribution": {"lengths": list(range(50, 151)), "counts": np.random.poisson(100, 101)},
            "basic_statistics": {
                "num_reads": 10000,
                "total_bases": 1500000,
                "min_length": 50,
                "max_length": 150,
                "mean_length": 150.0,
            },
        }

        ax = plot_quality_metrics(qc_data)
        assert ax is not None
        plt.close("all")

    def test_quality_metrics_partial_data(self):
        """Test quality metrics plot with partial data."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Only some QC metrics
        qc_data = {
            "per_base_quality": {"positions": list(range(50)), "mean_qualities": np.random.uniform(25, 35, 50)},
            "basic_statistics": {"num_reads": 5000, "total_bases": 750000, "mean_length": 150.0},
        }

        ax = plot_quality_metrics(qc_data)
        assert ax is not None
        plt.close("all")

    def test_quality_metrics_empty_data(self):
        """Test quality metrics plot with empty data."""
        qc_data = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_quality_metrics(qc_data)


class TestPlotAdapterContent:
    """Test plot_adapter_content function."""

    def test_basic_adapter_content_plot(self):
        """Test basic adapter content plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample adapter data
        adapter_data = {
            "Illumina_Universal_Adapter": [0.5, 1.2, 0.8, 2.1, 1.5],
            "Illumina_Small_RNA_3p_Adapter": [0.2, 0.5, 0.3, 0.8, 0.4],
            "TruSeq_Adapter_Index_1": [0.1, 0.3, 0.2, 0.4, 0.2],
        }

        ax = plot_adapter_content(adapter_data)
        assert ax is not None
        plt.close("all")

    def test_adapter_content_with_output_path(self, tmp_path: Path):
        """Test adapter content plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        adapter_data = {"Adapter1": [1.0, 1.5], "Adapter2": [0.5, 0.8]}
        output_path = tmp_path / "adapter_content.png"

        ax = plot_adapter_content(adapter_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_adapter_content_empty_data(self):
        """Test adapter content plot with empty data."""
        adapter_data = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_adapter_content(adapter_data)


class TestPlotGcDistribution:
    """Test plot_gc_distribution function."""

    def test_basic_gc_distribution_plot(self):
        """Test basic GC distribution plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample GC content data
        gc_data = np.random.normal(45, 10, 1000)
        gc_data = np.clip(gc_data, 0, 100)  # Ensure valid range

        ax = plot_gc_distribution(gc_data)
        assert ax is not None
        plt.close("all")

    def test_gc_distribution_with_output_path(self, tmp_path: Path):
        """Test GC distribution plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        gc_data = [40, 45, 50, 35, 55, 42, 48, 38, 52, 46]
        output_path = tmp_path / "gc_dist.png"

        ax = plot_gc_distribution(gc_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_gc_distribution_empty_data(self):
        """Test GC distribution plot with empty data."""
        gc_data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_gc_distribution(gc_data)


class TestPlotLengthDistribution:
    """Test plot_length_distribution function."""

    def test_basic_length_distribution_plot(self):
        """Test basic length distribution plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample length data
        length_data = np.random.normal(150, 20, 1000).astype(int)
        length_data = np.clip(length_data, 50, 250)  # Ensure reasonable range

        ax = plot_length_distribution(length_data)
        assert ax is not None
        plt.close("all")

    def test_length_distribution_with_output_path(self, tmp_path: Path):
        """Test length distribution plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        length_data = [100, 150, 120, 140, 160, 130, 145, 155, 125, 135]
        output_path = tmp_path / "length_dist.png"

        ax = plot_length_distribution(length_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_length_distribution_empty_data(self):
        """Test length distribution plot with empty data."""
        length_data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_length_distribution(length_data)
