"""Tests for population genetics visualization functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna.population_viz import (
    plot_demographic_comparison,
    plot_diversity_comparison,
    plot_fst_comparison,
    plot_kinship_matrix,
    plot_linkage_disequilibrium_decay,
    plot_neutrality_test_summary,
    plot_pca_results,
    plot_site_frequency_spectrum,
    plot_summary_statistics_grid,
    plot_tajimas_d_comparison,
)


class TestPlotDiversityComparison:
    """Test diversity comparison plotting."""

    def test_basic_plot(self, tmp_path: Path):
        """Test basic diversity comparison plot."""
        diversity = {
            "Neutral": 0.01,
            "High": 0.05,
            "Low": 0.001,
        }
        
        output_path = tmp_path / "diversity.png"
        fig = plot_diversity_comparison(
            diversity,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()
        assert len(fig.axes) == 1

    def test_plot_without_saving(self):
        """Test plotting without saving."""
        diversity = {"Scenario1": 0.01, "Scenario2": 0.02}
        fig = plot_diversity_comparison(diversity)
        assert fig is not None


class TestPlotTajimasDComparison:
    """Test Tajima's D comparison plotting."""

    def test_basic_plot(self, tmp_path: Path):
        """Test basic Tajima's D plot."""
        tajimas_d = {
            "Neutral": 0.1,
            "Negative": -1.5,
            "Positive": 1.5,
        }
        
        output_path = tmp_path / "tajimas_d.png"
        fig = plot_tajimas_d_comparison(
            tajimas_d,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()


class TestPlotFstComparison:
    """Test Fst comparison plotting."""

    def test_basic_plot(self, tmp_path: Path):
        """Test basic Fst plot."""
        fst_values = {
            "Low": 0.05,
            "Moderate": 0.15,
            "High": 0.3,
        }
        
        output_path = tmp_path / "fst.png"
        fig = plot_fst_comparison(
            fst_values,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()


class TestPlotPCAResults:
    """Test PCA results plotting."""

    def test_pca_plot(self, tmp_path: Path):
        """Test PCA plotting."""
        pca_result = {
            "status": "success",
            "pcs": [[0.1, 0.2, 0.3] for _ in range(50)],
            "explained_variance_ratio": [0.3, 0.2, 0.1, 0.05, 0.03],
            "n_components": 5,
        }
        
        output_path = tmp_path / "pca.png"
        fig = plot_pca_results(
            pca_result,
            output_path=str(output_path),
            n_components=5
        )
        
        assert fig is not None
        assert output_path.exists()
        assert len(fig.axes) == 3  # Three subplots

    def test_pca_plot_failure(self):
        """Test PCA plot with failed status."""
        pca_result = {"status": "failed"}
        
        with pytest.raises(ValueError, match="status is not 'success'"):
            plot_pca_results(pca_result)


class TestPlotKinshipMatrix:
    """Test kinship matrix plotting."""

    def test_kinship_plot(self, tmp_path: Path):
        """Test kinship matrix plotting."""
        kinship_result = {
            "status": "success",
            "kinship_matrix": [[1.0 if i == j else 0.1 for j in range(20)] for i in range(20)],
            "method": "vanraden",
        }
        
        output_path = tmp_path / "kinship.png"
        fig = plot_kinship_matrix(
            kinship_result,
            output_path=str(output_path),
            max_samples=20
        )
        
        assert fig is not None
        assert output_path.exists()

    def test_kinship_plot_failure(self):
        """Test kinship plot with failed status."""
        kinship_result = {"status": "failed"}
        
        with pytest.raises(ValueError, match="status is not 'success'"):
            plot_kinship_matrix(kinship_result)


class TestPlotSiteFrequencySpectrum:
    """Test site frequency spectrum plotting."""

    def test_sfs_plot(self, tmp_path: Path):
        """Test SFS plotting."""
        sfs = [50, 30, 20, 10, 5]  # Counts per frequency bin
        
        output_path = tmp_path / "sfs.png"
        fig = plot_site_frequency_spectrum(
            sfs,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()


class TestPlotNeutralityTestSummary:
    """Test neutrality test summary plotting."""

    def test_neutrality_summary_plot(self, tmp_path: Path):
        """Test neutrality test summary plotting."""
        neutrality_data = {
            "Neutral": {
                "tajimas_d": 0.1,
                "pi_theta_ratio": 1.0,
                "nucleotide_diversity": 0.01,
                "segregating_sites": 100,
            },
            "Bottleneck": {
                "tajimas_d": -1.5,
                "pi_theta_ratio": 0.8,
                "nucleotide_diversity": 0.005,
                "segregating_sites": 50,
            },
        }
        
        output_path = tmp_path / "neutrality_summary.png"
        fig = plot_neutrality_test_summary(
            neutrality_data,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()
        assert len(fig.axes) == 4  # Four subplots


class TestPlotDemographicComparison:
    """Test demographic comparison plotting."""

    def test_demographic_plot(self, tmp_path: Path):
        """Test demographic comparison plotting."""
        demographic_results = {
            "bottleneck": {
                "estimated_ne": 15.0,
                "observed_diversity": 0.025,
            },
            "expansion": {
                "estimated_ne": 250.0,
                "observed_diversity": 0.005,
            },
        }
        
        output_path = tmp_path / "demographic.png"
        fig = plot_demographic_comparison(
            demographic_results,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()
        assert len(fig.axes) == 2  # Two subplots


class TestPlotSummaryStatisticsGrid:
    """Test summary statistics grid plotting."""

    def test_summary_grid_plot(self, tmp_path: Path):
        """Test summary statistics grid plotting."""
        summary_stats = {
            "Neutral": {
                "nucleotide_diversity": 0.01,
                "segregating_sites": 100,
                "wattersons_theta": 0.01,
                "tajimas_d": 0.1,
                "sample_size": 30,
                "sequence_length": 1000,
            },
            "High": {
                "nucleotide_diversity": 0.05,
                "segregating_sites": 500,
                "wattersons_theta": 0.05,
                "tajimas_d": -0.5,
                "sample_size": 30,
                "sequence_length": 1000,
            },
        }
        
        output_path = tmp_path / "summary_grid.png"
        fig = plot_summary_statistics_grid(
            summary_stats,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()
        assert len(fig.axes) == 6  # Six subplots


class TestPlotLinkageDisequilibriumDecay:
    """Test LD decay plotting."""

    def test_ld_decay_plot(self, tmp_path: Path):
        """Test LD decay plotting."""
        ld_values = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
        distances = [0, 1, 2, 3, 4, 5]
        
        output_path = tmp_path / "ld_decay.png"
        fig = plot_linkage_disequilibrium_decay(
            ld_values,
            distances=distances,
            output_path=str(output_path)
        )
        
        assert fig is not None
        assert output_path.exists()

    def test_ld_decay_plot_no_distances(self, tmp_path: Path):
        """Test LD decay plot without distances."""
        ld_values = [0.5, 0.4, 0.3, 0.2]
        
        fig = plot_linkage_disequilibrium_decay(
            ld_values,
            output_path=str(tmp_path / "ld_decay2.png")
        )
        
        assert fig is not None

