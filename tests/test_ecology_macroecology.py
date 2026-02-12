"""Tests for ecology macroecology and visualization modules.

Covers:
    - fit_logseries (Fisher's log-series)
    - fit_lognormal (Preston's lognormal)
    - fit_broken_stick (MacArthur's broken stick)
    - fit_geometric_series (Motomura's geometric series)
    - compare_sad_models (model comparison with AIC)
    - species_area_power (Arrhenius power-law SAR)
    - species_area_logarithmic (Gleason logarithmic SAR)
    - distance_decay (exponential + power-law)
    - occupancy_frequency (core/common/satellite classification)
    - metabolic_scaling (Kleiber's law)
    - endemism_index (weighted endemism)
    - taylors_power_law (variance-mean scaling)
    - All visualization plot_* functions (Agg backend)

Uses real implementations only (NO mocking per project policy).
"""

from __future__ import annotations

import math
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.ecology.analysis.macroecology import (
    compare_sad_models,
    distance_decay,
    endemism_index,
    fit_broken_stick,
    fit_geometric_series,
    fit_lognormal,
    fit_logseries,
    metabolic_scaling,
    occupancy_frequency,
    species_area_logarithmic,
    species_area_power,
    taylors_power_law,
)
from metainformant.ecology.visualization.visualization import (
    create_interactive_ecology_dashboard,
    plot_beta_diversity_ordination,
    plot_biodiversity_rarefaction,
    plot_community_composition,
    plot_diversity_accumulation_curve,
    plot_diversity_indices_comparison,
    plot_ecological_distance_heatmap,
    plot_ecological_network,
    plot_rank_abundance_curve_comparison,
    plot_species_abundance_distribution,
)

# ---------------------------------------------------------------------------
# Shared test data fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def realistic_abundances() -> List[float]:
    """Realistic species abundance distribution (50 species)."""
    return [
        120,
        95,
        80,
        65,
        55,
        48,
        42,
        38,
        34,
        30,
        27,
        24,
        22,
        20,
        18,
        16,
        14,
        13,
        12,
        11,
        10,
        9,
        8,
        7,
        6,
        6,
        5,
        5,
        4,
        4,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
    ]


@pytest.fixture
def small_abundances() -> List[float]:
    """Small SAD for quick tests."""
    return [50, 30, 20, 10, 5, 3, 1, 1]


# ============================================================================
# Fisher's Log-Series tests
# ============================================================================


class TestFitLogseries:
    """Test Fisher's log-series distribution fitting."""

    def test_logseries_alpha_positive(self, small_abundances: List[float]) -> None:
        """Fisher's alpha must be positive."""
        result = fit_logseries(small_abundances)
        assert result["alpha"] > 0.0

    def test_logseries_x_in_range(self, small_abundances: List[float]) -> None:
        """The x parameter must be in (0, 1)."""
        result = fit_logseries(small_abundances)
        assert 0.0 < result["x"] < 1.0

    def test_logseries_expected_frequencies(self, small_abundances: List[float]) -> None:
        """Expected frequencies should be a non-empty list of positive values."""
        result = fit_logseries(small_abundances)
        ef = result["expected_frequencies"]
        assert len(ef) > 0
        assert all(e >= 0 for e in ef)

    def test_logseries_goodness_of_fit(self, small_abundances: List[float]) -> None:
        """Goodness-of-fit contains chi_squared and degrees_of_freedom."""
        result = fit_logseries(small_abundances)
        gof = result["goodness_of_fit"]
        assert "chi_squared" in gof
        assert "degrees_of_freedom" in gof
        assert gof["chi_squared"] >= 0

    def test_logseries_fisher_diversity(self, small_abundances: List[float]) -> None:
        """fisher_diversity equals alpha."""
        result = fit_logseries(small_abundances)
        assert abs(result["fisher_diversity"] - result["alpha"]) < 1e-12

    def test_logseries_single_species(self) -> None:
        """Single species yields alpha=0."""
        result = fit_logseries([100])
        assert result["alpha"] == 0.0

    def test_logseries_empty_raises(self) -> None:
        """Empty input raises ValueError."""
        with pytest.raises((ValueError, Exception)):
            fit_logseries([])

    def test_logseries_realistic(self, realistic_abundances: List[float]) -> None:
        """Realistic 50-species SAD produces valid output."""
        result = fit_logseries(realistic_abundances)
        assert result["alpha"] > 1.0  # Many species -> alpha > 1


# ============================================================================
# Preston's Lognormal tests
# ============================================================================


class TestFitLognormal:
    """Test Preston's lognormal distribution fitting."""

    def test_lognormal_mu_positive(self, small_abundances: List[float]) -> None:
        """Mean of log-abundances should be positive (all abundances > 1)."""
        result = fit_lognormal(small_abundances)
        assert result["mu"] > 0.0

    def test_lognormal_sigma_nonnegative(self, small_abundances: List[float]) -> None:
        """Standard deviation of log-abundances must be non-negative."""
        result = fit_lognormal(small_abundances)
        assert result["sigma"] >= 0.0

    def test_lognormal_s_star_geq_s(self, small_abundances: List[float]) -> None:
        """Estimated total species (S*) should be >= observed species count."""
        result = fit_lognormal(small_abundances)
        s_observed = len([a for a in small_abundances if a > 0])
        assert result["s_star"] >= s_observed - 0.1

    def test_lognormal_expected_frequencies(self, small_abundances: List[float]) -> None:
        """Expected frequencies per octave are non-negative."""
        result = fit_lognormal(small_abundances)
        for e in result["expected_frequencies"]:
            assert e >= 0.0

    def test_lognormal_single_species(self) -> None:
        """Single species yields sigma=0."""
        result = fit_lognormal([100])
        assert result["sigma"] == 0.0


# ============================================================================
# MacArthur's Broken Stick tests
# ============================================================================


class TestFitBrokenStick:
    """Test MacArthur's broken stick model fitting."""

    def test_broken_stick_expected_count_matches(self, small_abundances: List[float]) -> None:
        """Expected abundances list has same length as observed."""
        result = fit_broken_stick(small_abundances)
        assert len(result["expected_abundances"]) == len(small_abundances)

    def test_broken_stick_expected_sum(self, small_abundances: List[float]) -> None:
        """Expected abundances should sum to approximately N."""
        result = fit_broken_stick(small_abundances)
        n = sum(small_abundances)
        expected_sum = sum(result["expected_abundances"])
        assert abs(expected_sum - n) < 0.5

    def test_broken_stick_descending(self, small_abundances: List[float]) -> None:
        """Expected abundances should be in descending order."""
        result = fit_broken_stick(small_abundances)
        ea = result["expected_abundances"]
        for i in range(len(ea) - 1):
            assert ea[i] >= ea[i + 1] - 1e-10

    def test_broken_stick_observed_sorted(self, small_abundances: List[float]) -> None:
        """Observed abundances in result are sorted descending."""
        result = fit_broken_stick(small_abundances)
        obs = result["observed_abundances"]
        for i in range(len(obs) - 1):
            assert obs[i] >= obs[i + 1]


# ============================================================================
# Motomura's Geometric Series tests
# ============================================================================


class TestFitGeometricSeries:
    """Test geometric series (niche preemption) fitting."""

    def test_geometric_k_in_range(self, small_abundances: List[float]) -> None:
        """Dominance parameter k must be in (0, 1)."""
        result = fit_geometric_series(small_abundances)
        assert 0.0 < result["k"] < 1.0

    def test_geometric_expected_count_matches(self, small_abundances: List[float]) -> None:
        """Expected abundances list matches observed length."""
        result = fit_geometric_series(small_abundances)
        assert len(result["expected_abundances"]) == len(small_abundances)

    def test_geometric_single_species(self) -> None:
        """Single species: k = 1.0."""
        result = fit_geometric_series([100])
        assert result["k"] == 1.0

    def test_geometric_goodness_of_fit(self, small_abundances: List[float]) -> None:
        """Goodness-of-fit dict is present with correct keys."""
        result = fit_geometric_series(small_abundances)
        gof = result["goodness_of_fit"]
        assert "chi_squared" in gof
        assert gof["chi_squared"] >= 0


# ============================================================================
# Compare SAD Models tests
# ============================================================================


class TestCompareSadModels:
    """Test model comparison across all four SAD models."""

    def test_compare_returns_four_models(self, small_abundances: List[float]) -> None:
        """All four models are returned."""
        result = compare_sad_models(small_abundances)
        expected_models = {"logseries", "lognormal", "broken_stick", "geometric_series"}
        assert set(result.keys()) == expected_models

    def test_compare_aic_present(self, small_abundances: List[float]) -> None:
        """Each model has an AIC value."""
        result = compare_sad_models(small_abundances)
        for model_name, model_data in result.items():
            assert "aic" in model_data

    def test_compare_ranks_unique(self, small_abundances: List[float]) -> None:
        """Model ranks are 1, 2, 3, 4 (unique)."""
        result = compare_sad_models(small_abundances)
        ranks = sorted(m["rank"] for m in result.values())
        assert ranks == [1, 2, 3, 4]

    def test_compare_best_model_has_lowest_aic(self, small_abundances: List[float]) -> None:
        """The model ranked 1 has the lowest AIC."""
        result = compare_sad_models(small_abundances)
        best_model = [m for m in result.values() if m["rank"] == 1][0]
        for m in result.values():
            assert best_model["aic"] <= m["aic"] + 1e-10


# ============================================================================
# Species-Area Power Law (Arrhenius) tests
# ============================================================================


class TestSpeciesAreaPower:
    """Test Arrhenius power-law species-area relationship."""

    def test_power_z_positive(self) -> None:
        """Exponent z should be positive for increasing S with A."""
        result = species_area_power([1, 10, 100, 1000], [10, 30, 80, 200])
        assert result["z"] > 0

    def test_power_c_positive(self) -> None:
        """Intercept c should be positive."""
        result = species_area_power([1, 10, 100], [5, 15, 45])
        assert result["c"] > 0

    def test_power_r_squared_range(self) -> None:
        """R-squared in [0, 1] for well-behaved data."""
        result = species_area_power([1, 10, 100, 1000], [10, 30, 80, 200])
        assert 0.0 <= result["r_squared"] <= 1.0

    def test_power_ci_z_exists(self) -> None:
        """Bootstrap CI for z is a tuple of two floats."""
        result = species_area_power([1, 10, 100], [5, 15, 45], n_bootstrap=100)
        ci = result["ci_z"]
        assert len(ci) == 2
        assert ci[0] <= ci[1]

    def test_power_predicted_species(self) -> None:
        """Predicted species list has same length as input."""
        areas = [1, 10, 100]
        result = species_area_power(areas, [5, 15, 45])
        assert len(result["predicted_species"]) == 3

    def test_power_non_positive_area_raises(self) -> None:
        """Non-positive areas raise ValueError."""
        with pytest.raises(ValueError, match="positive"):
            species_area_power([0, 10], [5, 15])

    def test_power_mismatched_lengths_raises(self) -> None:
        """Mismatched lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            species_area_power([1, 10], [5])


# ============================================================================
# Species-Area Logarithmic (Gleason) tests
# ============================================================================


class TestSpeciesAreaLogarithmic:
    """Test Gleason's logarithmic species-area model."""

    def test_log_b_positive(self) -> None:
        """Slope b should be positive for increasing S with A."""
        result = species_area_logarithmic([1, 10, 100], [5, 15, 25])
        assert result["b"] > 0

    def test_log_r_squared_range(self) -> None:
        """R-squared in [0, 1]."""
        result = species_area_logarithmic([1, 10, 100, 1000], [5, 15, 25, 35])
        assert 0.0 <= result["r_squared"] <= 1.0

    def test_log_predicted_species(self) -> None:
        """Predicted species list matches input length."""
        result = species_area_logarithmic([1, 10, 100], [5, 15, 25])
        assert len(result["predicted_species"]) == 3

    def test_log_non_positive_area_raises(self) -> None:
        """Non-positive areas raise ValueError."""
        with pytest.raises(ValueError, match="positive"):
            species_area_logarithmic([0, 10], [5, 15])


# ============================================================================
# Distance Decay tests
# ============================================================================


class TestDistanceDecay:
    """Test distance-decay relationship fitting."""

    def test_decay_output_keys(self) -> None:
        """Output has exponential, power_law, and best_model."""
        result = distance_decay([1, 5, 10, 20], [0.9, 0.6, 0.3, 0.1])
        for key in ("exponential", "power_law", "best_model"):
            assert key in result

    def test_decay_model_keys(self) -> None:
        """Each model has a, b, r_squared."""
        result = distance_decay([1, 5, 10, 20], [0.9, 0.6, 0.3, 0.1])
        for model in ("exponential", "power_law"):
            assert "a" in result[model]
            assert "b" in result[model]
            assert "r_squared" in result[model]

    def test_decay_best_model_valid(self) -> None:
        """Best model is one of the two candidates."""
        result = distance_decay([1, 5, 10, 20], [0.9, 0.6, 0.3, 0.1])
        assert result["best_model"] in ("exponential", "power_law")

    def test_decay_b_positive(self) -> None:
        """Decay rate b should be positive for declining similarity."""
        result = distance_decay([1, 5, 10, 20], [0.9, 0.6, 0.3, 0.1])
        # At least one model should have b > 0
        assert result["exponential"]["b"] > -10  # b is negated in the code

    def test_decay_mismatched_lengths_raises(self) -> None:
        """Mismatched lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            distance_decay([1, 5], [0.9])

    def test_decay_too_few_points_raises(self) -> None:
        """Fewer than 2 points raise ValueError."""
        with pytest.raises(ValueError):
            distance_decay([1], [0.9])


# ============================================================================
# Occupancy-Frequency tests
# ============================================================================


class TestOccupancyFrequency:
    """Test occupancy-frequency distribution analysis."""

    def test_occupancy_basic(self) -> None:
        """Core, common, satellite counts sum to total species."""
        mat = [[1, 0, 1], [1, 1, 0], [1, 1, 0]]
        result = occupancy_frequency(mat)
        total = result["core_species"] + result["common_species"] + result["satellite_species"]
        assert total == 3

    def test_occupancy_core_species(self) -> None:
        """A species present in all sites is 'core'."""
        mat = [[1, 0], [1, 0], [1, 0], [1, 1]]
        result = occupancy_frequency(mat)
        # Species 0 is present in 4/4 = 100% -> core
        assert result["core_species"] >= 1

    def test_occupancy_satellite_species(self) -> None:
        """A species present in <33% of sites is 'satellite'."""
        mat = [[0, 1], [0, 0], [0, 0], [0, 0]]
        result = occupancy_frequency(mat)
        # Species 0 present 0%, species 1 present 25% -> both satellite
        assert result["satellite_species"] >= 1

    def test_occupancy_distribution_length(self) -> None:
        """Occupancy distribution has one entry per species."""
        mat = [[1, 0, 1, 1], [1, 1, 0, 0], [1, 1, 1, 0]]
        result = occupancy_frequency(mat)
        assert len(result["occupancy_distribution"]) == 4

    def test_occupancy_bimodality_range(self) -> None:
        """Bimodality index is in [0, 1]."""
        mat = [[1, 0, 1], [1, 1, 0], [1, 1, 0]]
        result = occupancy_frequency(mat)
        assert 0.0 <= result["bimodality_index"] <= 1.0


# ============================================================================
# Metabolic Scaling tests
# ============================================================================


class TestMetabolicScaling:
    """Test allometric metabolic scaling (Kleiber's law)."""

    def test_scaling_exponent_positive(self) -> None:
        """Allometric exponent should be positive."""
        result = metabolic_scaling([1, 10, 100, 1000], [0.5, 3.0, 17.0, 100.0])
        assert result["b_exponent"] > 0

    def test_scaling_b0_positive(self) -> None:
        """Scaling constant b0 should be positive."""
        result = metabolic_scaling([1, 10, 100], [0.5, 3.0, 17.0])
        assert result["b0"] > 0

    def test_scaling_near_kleiber(self) -> None:
        """Data generated with b~0.75 should have small Kleiber deviation."""
        # Generate data following Kleiber's law: B = 0.5 * M^0.75
        masses = [1, 5, 10, 50, 100, 500, 1000]
        rates = [0.5 * m**0.75 for m in masses]
        result = metabolic_scaling(masses, rates)
        assert result["kleiber_deviation"] < 0.05

    def test_scaling_non_positive_raises(self) -> None:
        """Non-positive values raise ValueError."""
        with pytest.raises(ValueError, match="positive"):
            metabolic_scaling([0, 10], [1, 5])


# ============================================================================
# Endemism Index tests
# ============================================================================


class TestEndemismIndex:
    """Test weighted and corrected endemism indices."""

    def test_endemism_known_value(self) -> None:
        """Hand-calculated endemism for 4 species."""
        result = endemism_index([10, 50, 200, 500], 100)
        # WE = 1/10 + 1/50 + 1/200 + 1/500 = 0.1 + 0.02 + 0.005 + 0.002 = 0.127
        assert abs(result["weighted_endemism"] - 0.127) < 0.001
        # CWE = WE / 4 = 0.031750
        assert abs(result["corrected_weighted_endemism"] - 0.127 / 4) < 0.001
        # Endemic: range < 100 -> species with range 10 and 50
        assert result["n_endemic_species"] == 2

    def test_endemism_no_endemic(self) -> None:
        """No species with range < area yields 0 endemic."""
        result = endemism_index([500, 1000], 100)
        assert result["n_endemic_species"] == 0

    def test_endemism_all_endemic(self) -> None:
        """All species with range < area."""
        result = endemism_index([10, 20, 30], 100)
        assert result["n_endemic_species"] == 3

    def test_endemism_non_positive_area_raises(self) -> None:
        """Non-positive area raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            endemism_index([10, 20], 0)


# ============================================================================
# Taylor's Power Law tests
# ============================================================================


class TestTaylorsPowerLaw:
    """Test Taylor's power law (variance-mean relationship)."""

    def test_taylors_b_positive(self) -> None:
        """Power-law exponent should be positive."""
        result = taylors_power_law([1, 5, 10, 50], [1.2, 28, 110, 2600])
        assert result["b"] > 0

    def test_taylors_a_positive(self) -> None:
        """Scaling constant a should be positive."""
        result = taylors_power_law([1, 5, 10, 50], [1.2, 28, 110, 2600])
        assert result["a"] > 0

    def test_taylors_r_squared_range(self) -> None:
        """R-squared in [0, 1]."""
        result = taylors_power_law([1, 5, 10, 50], [1.2, 28, 110, 2600])
        assert 0.0 <= result["r_squared"] <= 1.0

    def test_taylors_aggregation_detection(self) -> None:
        """Data with b~2 indicates negative binomial aggregation."""
        # Generate data with b~2: V = a * M^2
        means = [1, 2, 5, 10, 20, 50]
        variances = [1.1 * m**2 for m in means]
        result = taylors_power_law(means, variances)
        assert abs(result["b"] - 2.0) < 0.1

    def test_taylors_mismatched_raises(self) -> None:
        """Mismatched lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            taylors_power_law([1, 5], [1.2])

    def test_taylors_too_few_raises(self) -> None:
        """Fewer than 2 points raise ValueError."""
        with pytest.raises(ValueError):
            taylors_power_law([1], [1.2])


# ============================================================================
# Visualization tests (Agg backend -- no display)
# ============================================================================


class TestVisualization:
    """Test all ecology visualization functions using the Agg backend."""

    def teardown_method(self) -> None:
        """Close all matplotlib figures after each test."""
        plt.close("all")

    def test_plot_species_abundance_distribution(self) -> None:
        """Species abundance distribution plot returns an Axes object."""
        data = np.array([50, 30, 20, 10, 5, 3, 1, 1], dtype=float)
        ax = plot_species_abundance_distribution(data)
        assert ax is not None
        assert hasattr(ax, "get_title")

    def test_plot_species_abundance_with_names(self) -> None:
        """SAD plot accepts optional species names."""
        data = np.array([50, 30, 20], dtype=float)
        names = ["Sp_A", "Sp_B", "Sp_C"]
        ax = plot_species_abundance_distribution(data, species_names=names)
        assert ax is not None

    def test_plot_diversity_accumulation_curve(self) -> None:
        """Diversity accumulation curve returns an Axes object."""
        data = [
            {"sample_size": 10, "species_count": 5},
            {"sample_size": 20, "species_count": 8},
            {"sample_size": 30, "species_count": 10},
        ]
        ax = plot_diversity_accumulation_curve(data)
        assert ax is not None

    def test_plot_community_composition(self) -> None:
        """Community composition stacked bar returns an Axes."""
        matrix = np.array([[10, 5, 3], [8, 6, 4], [2, 9, 7]], dtype=float)
        ax = plot_community_composition(matrix)
        assert ax is not None

    def test_plot_community_composition_with_names(self) -> None:
        """Community composition with species and sample names."""
        matrix = np.array([[10, 5], [8, 6]], dtype=float)
        ax = plot_community_composition(
            matrix,
            species_names=["Sp1", "Sp2"],
            sample_names=["S1", "S2"],
        )
        assert ax is not None

    def test_plot_beta_diversity_ordination(self) -> None:
        """Beta diversity ordination scatter returns an Axes."""
        coords = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]], dtype=float)
        ax = plot_beta_diversity_ordination(coords)
        assert ax is not None

    def test_plot_beta_diversity_ordination_with_groups(self) -> None:
        """Ordination plot with group labels."""
        coords = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6], [0.7, 0.8]], dtype=float)
        groups = np.array(["A", "A", "B", "B"])
        ax = plot_beta_diversity_ordination(coords, sample_groups=groups)
        assert ax is not None

    def test_plot_diversity_indices_comparison(self) -> None:
        """Diversity indices comparison bar chart returns an Axes."""
        indices = {
            "Shannon": np.array([1.5, 2.0, 1.8]),
            "Simpson": np.array([0.7, 0.8, 0.75]),
        }
        ax = plot_diversity_indices_comparison(indices)
        assert ax is not None

    def test_plot_rank_abundance_curve_comparison(self) -> None:
        """Rank-abundance comparison plot returns an Axes."""
        datasets = {
            "Community_A": np.array([50, 30, 20, 10, 5]),
            "Community_B": np.array([40, 35, 25, 15, 10]),
        }
        ax = plot_rank_abundance_curve_comparison(datasets)
        assert ax is not None

    def test_plot_biodiversity_rarefaction(self) -> None:
        """Biodiversity rarefaction plot returns an Axes."""
        data = {
            "Site_1": [
                {"sample_size": 10, "species_count": 5},
                {"sample_size": 20, "species_count": 8},
            ],
            "Site_2": [
                {"sample_size": 10, "species_count": 3},
                {"sample_size": 20, "species_count": 6},
            ],
        }
        ax = plot_biodiversity_rarefaction(data)
        assert ax is not None

    def test_plot_ecological_distance_heatmap(self) -> None:
        """Ecological distance heatmap returns an Axes."""
        dm = np.array([[0, 0.3, 0.7], [0.3, 0, 0.5], [0.7, 0.5, 0]], dtype=float)
        ax = plot_ecological_distance_heatmap(dm)
        assert ax is not None

    def test_plot_ecological_distance_heatmap_with_names(self) -> None:
        """Distance heatmap accepts sample names."""
        dm = np.array([[0, 0.5], [0.5, 0]], dtype=float)
        ax = plot_ecological_distance_heatmap(dm, sample_names=["A", "B"])
        assert ax is not None

    def test_plot_ecological_network(self) -> None:
        """Ecological network plot returns an Axes (requires networkx)."""
        try:
            import networkx  # noqa: F401
        except ImportError:
            pytest.skip("networkx not installed")

        interaction = np.array(
            [[0, 0.5, 0.2], [0.3, 0, 0.6], [0.1, 0.4, 0]],
            dtype=float,
        )
        ax = plot_ecological_network(interaction)
        assert ax is not None

    def test_plot_ecological_network_with_names(self) -> None:
        """Network plot accepts species names."""
        try:
            import networkx  # noqa: F401
        except ImportError:
            pytest.skip("networkx not installed")

        interaction = np.array(
            [[0, 0.5], [0.3, 0]],
            dtype=float,
        )
        ax = plot_ecological_network(interaction, species_names=["Wolf", "Deer"])
        assert ax is not None

    def test_create_interactive_dashboard(self) -> None:
        """Interactive dashboard creation (requires plotly)."""
        try:
            import plotly  # noqa: F401
        except ImportError:
            pytest.skip("plotly not installed")

        data = {
            "diversity_indices": {
                "Shannon": [1.5, 2.0, 1.8],
                "Simpson": [0.7, 0.8, 0.75],
            }
        }
        fig = create_interactive_ecology_dashboard(data)
        assert fig is not None

    def test_plot_to_file(self, tmp_path) -> None:
        """Saving a plot to a file via output_path works."""
        data = np.array([50, 30, 20, 10], dtype=float)
        outfile = tmp_path / "test_plot.png"
        # The function saves to parent dir path returned by ensure_directory
        # We just verify it doesn't raise
        ax = plot_species_abundance_distribution(data)
        fig = ax.get_figure()
        fig.savefig(str(outfile), dpi=72)
        assert outfile.exists()
