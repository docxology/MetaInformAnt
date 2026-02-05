"""Comprehensive tests for the visualization module.

Tests all new functionality: themes, palettes, composite dashboards,
interactive plots, and quality control submodules.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
from pathlib import Path

import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Themes tests
# ---------------------------------------------------------------------------


class TestThemes:
    """Tests for themes module."""

    def test_list_themes(self):
        from metainformant.visualization.themes import list_themes

        themes = list_themes()
        assert isinstance(themes, list)
        assert "publication" in themes
        assert "presentation" in themes
        assert "poster" in themes
        assert "dark" in themes
        assert "minimal" in themes
        assert "colorblind" in themes

    def test_get_theme(self):
        from metainformant.visualization.themes import get_theme

        params = get_theme("publication")
        assert isinstance(params, dict)
        assert "figure.dpi" in params
        assert params["figure.dpi"] == 300

    def test_get_theme_unknown_raises(self):
        from metainformant.visualization.themes import get_theme

        with pytest.raises(KeyError, match="Unknown theme"):
            get_theme("nonexistent_theme")

    def test_apply_and_reset_theme(self):
        from metainformant.visualization.themes import apply_theme, reset_theme

        original_dpi = matplotlib.rcParams["figure.dpi"]
        apply_theme("publication")
        assert matplotlib.rcParams["figure.dpi"] == 300
        reset_theme()
        # After reset, DPI reverts to default (not necessarily original if it was changed)
        assert matplotlib.rcParams["figure.dpi"] != 300 or original_dpi == 300

    def test_theme_context_manager(self):
        from metainformant.visualization.themes import theme

        original_dpi = matplotlib.rcParams["figure.dpi"]
        with theme("poster"):
            assert matplotlib.rcParams["figure.dpi"] == 150
        # Context restored
        assert matplotlib.rcParams["figure.dpi"] == original_dpi

    def test_register_custom_theme(self):
        from metainformant.visualization.themes import get_theme, register_theme

        register_theme("custom_test", {"figure.dpi": 999})
        params = get_theme("custom_test")
        assert params["figure.dpi"] == 999


# ---------------------------------------------------------------------------
# Palettes tests
# ---------------------------------------------------------------------------


class TestPalettes:
    """Tests for palettes module."""

    def test_chromosome_palette_default(self):
        from metainformant.visualization.palettes import chromosome_palette

        colors = chromosome_palette()
        assert len(colors) == 25  # 22 autosomes + X, Y, MT

    def test_chromosome_palette_specific(self):
        from metainformant.visualization.palettes import chromosome_palette

        colors = chromosome_palette(["1", "2", "X"])
        assert len(colors) == 3
        assert all(c.startswith("#") for c in colors)

    def test_categorical_wong(self):
        from metainformant.visualization.palettes import categorical

        colors = categorical(5, palette="wong")
        assert len(colors) == 5
        assert all(c.startswith("#") for c in colors)

    def test_categorical_fallback_large_n(self):
        from metainformant.visualization.palettes import categorical

        colors = categorical(15, palette="wong")  # more than 8
        assert len(colors) == 15

    def test_expression_gradient(self):
        from metainformant.visualization.palettes import expression_gradient

        cmap = expression_gradient(256)
        assert hasattr(cmap, "__call__")

    def test_significance_palette(self):
        from metainformant.visualization.palettes import significance_palette

        palette = significance_palette()
        assert "highly_significant" in palette
        assert "not_significant" in palette

    def test_significance_color(self):
        from metainformant.visualization.palettes import significance_color

        assert significance_color(0.0001) == "#D32F2F"
        assert significance_color(0.005) == "#FF9800"
        assert significance_color(0.03) == "#FFC107"
        assert significance_color(0.1) == "#9E9E9E"

    def test_heatmap_cmap_valid(self):
        from metainformant.visualization.palettes import heatmap_cmap

        cmap = heatmap_cmap("expression")
        assert cmap is not None

    def test_heatmap_cmap_invalid_raises(self):
        from metainformant.visualization.palettes import heatmap_cmap

        with pytest.raises(KeyError):
            heatmap_cmap("invalid_name")

    def test_alternating_pair(self):
        from metainformant.visualization.palettes import alternating_pair

        colors = alternating_pair(2)
        assert len(colors) == 4  # 2 pairs


# ---------------------------------------------------------------------------
# Composite dashboards tests
# ---------------------------------------------------------------------------


class TestComposite:
    """Tests for composite dashboard functions."""

    def test_multi_panel(self, tmp_path: Path):
        from metainformant.visualization.composite import multi_panel

        def plot_func(ax, data=None):
            ax.plot([1, 2, 3], [1, 2, 3])

        panels = [
            {"func": plot_func, "kwargs": {}, "title": "Panel A"},
            {"func": plot_func, "kwargs": {}, "title": "Panel B"},
        ]
        output = tmp_path / "multi_panel.png"
        fig = multi_panel(panels, ncols=2, output_path=output)
        assert fig is not None
        assert output.exists()
        plt.close(fig)

    def test_genomic_overview(self, tmp_path: Path):
        from metainformant.visualization.composite import genomic_overview

        data = {
            "pca_data": np.random.randn(50, 2),
            "pvalues": np.random.uniform(0.0001, 1, 100),
            "gc_content": np.random.normal(50, 10, 200),
            "quality_scores": np.random.normal(30, 5, 200),
        }
        output = tmp_path / "genomic_overview.png"
        fig = genomic_overview(data, output_path=output)
        assert fig is not None
        assert output.exists()
        plt.close(fig)

    def test_qc_summary(self, tmp_path: Path):
        from metainformant.visualization.composite import qc_summary

        qc_data = {
            "read_lengths": np.random.normal(150, 10, 500).astype(int),
            "gc_content": np.random.normal(50, 8, 500),
            "quality_scores": np.random.normal(35, 3, 500),
            "mapping_rates": {"SampleA": 95.2, "SampleB": 92.1},
        }
        output = tmp_path / "qc_summary.png"
        fig = qc_summary(qc_data, output_path=output)
        assert fig is not None
        assert output.exists()
        plt.close(fig)


# ---------------------------------------------------------------------------
# Interactive plots tests
# ---------------------------------------------------------------------------


class TestInteractive:
    """Tests for interactive plot wrappers."""

    def test_interactive_scatter_fallback(self, tmp_path: Path):
        from metainformant.visualization.interactive import interactive_scatter

        df = pd.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
        output = tmp_path / "scatter.png"
        result = interactive_scatter(df, "x", "y", output_path=output)
        # May return Plotly Figure or Axes
        assert result is not None

    def test_interactive_heatmap(self, tmp_path: Path):
        from metainformant.visualization.interactive import interactive_heatmap

        data = np.random.randn(5, 5)
        output = tmp_path / "heatmap.png"
        result = interactive_heatmap(data, output_path=output)
        assert result is not None

    def test_interactive_volcano(self, tmp_path: Path):
        from metainformant.visualization.interactive import interactive_volcano

        df = pd.DataFrame(
            {
                "log2FoldChange": np.random.uniform(-3, 3, 100),
                "pvalue": np.random.uniform(0.0001, 1, 100),
                "gene": [f"Gene{i}" for i in range(100)],
            }
        )
        output = tmp_path / "volcano.png"
        result = interactive_volcano(df, output_path=output)
        assert result is not None

    def test_interactive_manhattan(self, tmp_path: Path):
        from metainformant.visualization.interactive import interactive_manhattan

        df = pd.DataFrame(
            {
                "chromosome": np.repeat([str(i) for i in range(1, 6)], 20),
                "position": np.tile(np.arange(1000, 21000, 1000), 5),
                "pvalue": np.random.uniform(1e-8, 1, 100),
            }
        )
        output = tmp_path / "manhattan.png"
        result = interactive_manhattan(df, output_path=output)
        assert result is not None


# ---------------------------------------------------------------------------
# Quality submodules tests
# ---------------------------------------------------------------------------


class TestQualitySequencing:
    """Tests for quality_sequencing module."""

    def test_plot_quality_metrics(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_sequencing import plot_quality_metrics

        qc_data = {
            "per_base_quality": {
                "positions": list(range(150)),
                "mean_qualities": np.random.normal(35, 2, 150).tolist(),
            },
            "gc_content_distribution": {
                "bins": list(range(0, 101, 5)),
                "counts": np.random.randint(10, 100, 20).tolist(),
            },
            "sequence_length_distribution": {
                "lengths": list(range(140, 160)),
                "counts": np.random.randint(1, 50, 20).tolist(),
            },
            "basic_statistics": {
                "num_reads": 1000,
                "total_bases": 150000,
                "min_length": 140,
                "max_length": 155,
                "mean_length": 150,
            },
        }
        output = tmp_path / "quality_metrics.png"
        ax = plot_quality_metrics(qc_data, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")

    def test_plot_gc_distribution(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_sequencing import plot_gc_distribution

        gc_data = list(np.random.normal(50, 10, 200))
        output = tmp_path / "gc_dist.png"
        ax = plot_gc_distribution(gc_data, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")

    def test_plot_length_distribution(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_sequencing import plot_length_distribution

        lengths = list(np.random.normal(150, 5, 300).astype(int))
        output = tmp_path / "length_dist.png"
        ax = plot_length_distribution(lengths, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")

    def test_plot_kmer_profiles(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_sequencing import plot_kmer_profiles

        kmer_counts = {f"ACGT{i}": np.random.randint(100, 1000) for i in range(30)}
        output = tmp_path / "kmer.png"
        ax = plot_kmer_profiles(kmer_counts, top_n=15, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")


class TestQualityOmics:
    """Tests for quality_omics module."""

    def test_plot_singlecell_qc_metrics(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_omics import plot_singlecell_qc_metrics

        qc_metrics = {
            "n_counts": np.random.exponential(5000, 500),
            "n_genes": np.random.exponential(2000, 500),
            "percent_mito": np.random.beta(2, 20, 500) * 100,
        }
        output = tmp_path / "singlecell_qc.png"
        ax = plot_singlecell_qc_metrics(qc_metrics, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")

    def test_plot_multiomics_quality_overview(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_omics import plot_multiomics_quality_overview

        quality_reports = {
            "Transcriptomics": {"overall_quality": 0.92},
            "Proteomics": {"overall_quality": 0.85},
            "Metabolomics": {"overall_quality": 0.78},
        }
        output = tmp_path / "multiomics_qc.png"
        ax = plot_multiomics_quality_overview(quality_reports, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")


class TestQualityAssessment:
    """Tests for quality_assessment module."""

    def test_plot_coverage_uniformity(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_assessment import plot_coverage_uniformity

        coverage = np.random.poisson(30, 500).astype(float)
        output = tmp_path / "coverage.png"
        ax = plot_coverage_uniformity(coverage, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")

    def test_plot_error_profiles(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_assessment import plot_error_profiles

        profiles = {
            "mismatch": np.random.exponential(0.001, 150),
            "insertion": np.random.exponential(0.0005, 150),
            "deletion": np.random.exponential(0.0005, 150),
        }
        output = tmp_path / "errors.png"
        ax = plot_error_profiles(profiles, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")

    def test_plot_data_integrity_metrics(self, tmp_path: Path):
        from metainformant.visualization.analysis.quality_assessment import plot_data_integrity_metrics

        metrics = {
            "completeness": 0.95,
            "missing_rate": 0.02,
            "valid_entries": 0.99,
            "error_rate": 0.005,
        }
        output = tmp_path / "integrity.png"
        ax = plot_data_integrity_metrics(metrics, output_path=output)
        assert ax is not None
        assert output.exists()
        plt.close("all")


# ---------------------------------------------------------------------------
# Analysis __init__ exports test
# ---------------------------------------------------------------------------


class TestAnalysisExports:
    """Test that analysis subpackage exports all expected functions."""

    def test_all_exports_present(self):
        from metainformant.visualization import analysis

        # Dimensionality reduction
        assert hasattr(analysis, "plot_pca")
        assert hasattr(analysis, "plot_umap")
        assert hasattr(analysis, "plot_tsne")
        assert hasattr(analysis, "biplot")

        # Statistical
        assert hasattr(analysis, "histogram")
        assert hasattr(analysis, "box_plot")
        assert hasattr(analysis, "violin_plot")
        assert hasattr(analysis, "qq_plot")
        assert hasattr(analysis, "correlation_heatmap")
        assert hasattr(analysis, "roc_curve")

        # Time series
        assert hasattr(analysis, "plot_time_series")
        assert hasattr(analysis, "plot_autocorrelation")
        assert hasattr(analysis, "plot_forecast")

        # Information
        assert hasattr(analysis, "plot_entropy_profile")
        assert hasattr(analysis, "plot_mutual_information_matrix")

        # Quality - sequencing
        assert hasattr(analysis, "plot_quality_metrics")
        assert hasattr(analysis, "plot_gc_distribution")
        assert hasattr(analysis, "plot_kmer_profiles")

        # Quality - omics
        assert hasattr(analysis, "plot_vcf_quality_metrics")
        assert hasattr(analysis, "plot_singlecell_qc_metrics")

        # Quality - assessment
        assert hasattr(analysis, "plot_coverage_uniformity")
        assert hasattr(analysis, "plot_error_profiles")


# ---------------------------------------------------------------------------
# Main visualization __init__ exports test
# ---------------------------------------------------------------------------


class TestMainExports:
    """Test that main visualization package exports new modules."""

    def test_themes_exported(self):
        from metainformant.visualization import apply_theme, list_themes, theme, themes

        assert themes is not None
        assert callable(list_themes)
        assert callable(apply_theme)
        assert callable(theme)

    def test_palettes_exported(self):
        from metainformant.visualization import WONG, categorical, chromosome_palette, palettes

        assert palettes is not None
        assert callable(chromosome_palette)
        assert callable(categorical)
        assert isinstance(WONG, list)

    def test_composite_exported(self):
        from metainformant.visualization import composite, genomic_overview, multi_panel, qc_summary

        assert composite is not None
        assert callable(multi_panel)
        assert callable(genomic_overview)
        assert callable(qc_summary)

    def test_interactive_exported(self):
        from metainformant.visualization import interactive, interactive_scatter, interactive_volcano

        assert interactive is not None
        assert callable(interactive_scatter)
        assert callable(interactive_volcano)
