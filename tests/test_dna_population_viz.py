"""Tests for population genetics visualization functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna.population.visualization import (
    plot_allele_frequency_spectrum,
    plot_bootstrap_distribution,
    plot_demographic_comparison,
    plot_diversity_comparison,
    plot_fst_comparison,
    plot_fst_matrix,
    plot_fst_matrix,
    plot_hardy_weinberg_test,
    plot_heterozygosity_distribution,
    plot_kinship_matrix,
    plot_linkage_disequilibrium_decay,
    plot_neutrality_test_suite,
    plot_neutrality_test_summary,
    plot_outlier_detection,
    plot_pca_results,
    plot_permutation_test,
    plot_pi_vs_theta,
    plot_site_frequency_spectrum,
    plot_statistic_correlation_matrix,
    plot_statistic_distribution,
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
        
        from metainformant.core.utils.errors import ValidationError
        with pytest.raises(ValidationError, match="PCA result status is not 'success'"):
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
        
        from metainformant.core.utils.errors import ValidationError
        with pytest.raises(ValidationError, match="Kinship result status is not 'success'"):
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


def test_plot_allele_frequency_spectrum():
    """Test allele frequency spectrum plot."""
    sfs = [10, 5, 3, 2, 1]
    fig = plot_allele_frequency_spectrum(sfs)

    assert fig is not None


def test_plot_heterozygosity_distribution():
    """Test heterozygosity distribution plot."""
    het_values = [0.1, 0.2, 0.3, 0.4, 0.5]
    fig = plot_heterozygosity_distribution(het_values)

    assert fig is not None


def test_plot_statistic_distribution():
    """Test statistic distribution plot."""
    stats = {
        "scenario1": [1.0, 2.0, 3.0],
        "scenario2": [4.0, 5.0, 6.0],
    }
    fig = plot_statistic_distribution(stats, plot_type="histogram")

    assert fig is not None


def test_plot_pi_vs_theta():
    """Test π vs θ plot."""
    pi_values = [0.01, 0.02, 0.03]
    theta_values = [0.01, 0.02, 0.03]
    fig = plot_pi_vs_theta(pi_values, theta_values)

    assert fig is not None


def test_plot_statistic_correlation_matrix():
    """Test statistic correlation matrix plot."""
    stats = {
        "pi": [0.01, 0.02, 0.03],
        "theta": [0.01, 0.02, 0.03],
        "d": [-0.5, 0.0, 0.5],
    }
    fig = plot_statistic_correlation_matrix(stats)

    assert fig is not None


def test_plot_fst_matrix():
    """Test Fst matrix plot."""
    fst_matrix = [[0.0, 0.1, 0.2], [0.1, 0.0, 0.15], [0.2, 0.15, 0.0]]
    fig = plot_fst_matrix(fst_matrix)

    assert fig is not None


def test_plot_neutrality_test_suite():
    """Test neutrality test suite plot."""
    test_results = {
        "tajimas_d": {"statistic": -0.5},
        "fu_and_li_d": {"statistic": -0.3},
        "fay_wu_h": {"statistic": -0.2},
    }
    fig = plot_neutrality_test_suite(test_results)

    assert fig is not None


def test_plot_hardy_weinberg_test():
    """Test Hardy-Weinberg test plot."""
    hwe_results = [{
        "locus": "Test_Locus",
        "chi_square": 2.5,
        "p_value": 0.1,
        "degrees_of_freedom": 1,
        "hwe_deviated": False,
    }]
    fig = plot_hardy_weinberg_test(hwe_results)

    assert fig is not None


def test_plot_bootstrap_distribution():
    """Test bootstrap distribution plot."""
    bootstrap_values = [1.0, 1.5, 2.0, 2.5, 3.0]
    fig = plot_bootstrap_distribution(bootstrap_values, observed_value=2.0)

    assert fig is not None


def test_plot_permutation_test():
    """Test permutation test plot."""
    permuted_values = [0.1, 0.2, 0.3, 0.4, 0.5]
    observed_value = 0.8
    fig = plot_permutation_test(permuted_values, observed_value, p_value=0.01)

    assert fig is not None


def test_plot_outlier_detection():
    """Test outlier detection plot."""
    statistic_values = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0]
    outlier_indices = [3]
    fig = plot_outlier_detection(statistic_values, outlier_indices)

    assert fig is not None

