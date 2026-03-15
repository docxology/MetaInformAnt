"""Comprehensive tests for the ecology module.

Tests all analysis submodules: community, ordination, indicators, functional, macroecology.
Uses real data only (NO mocking per project policy).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import List

import pytest

# ──────────────────────────────────────────────────────────────────────────────
# Test data fixtures
# ──────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def simple_abundances() -> List[float]:
    """Simple species abundance vector."""
    return [10, 8, 6, 4, 2, 1, 1, 1]


@pytest.fixture
def equal_abundances() -> List[float]:
    """Perfectly even community."""
    return [10, 10, 10, 10]


@pytest.fixture
def single_species() -> List[float]:
    """Single dominant species."""
    return [100, 0, 0, 0]


@pytest.fixture
def community_matrix() -> List[List[float]]:
    """Multiple communities for beta diversity testing."""
    return [
        [10, 8, 6, 4, 2, 0, 0],  # Community A
        [0, 6, 8, 10, 4, 2, 1],  # Community B
        [10, 8, 6, 4, 2, 0, 0],  # Community C (same as A)
        [1, 1, 1, 1, 1, 1, 1],  # Community D (even)
    ]


@pytest.fixture
def presence_absence() -> List[List[int]]:
    """Presence-absence matrix for nestedness."""
    return [
        [1, 1, 1, 1, 1],  # Richest site
        [1, 1, 1, 1, 0],  # Nested subset
        [1, 1, 1, 0, 0],  # Further nested
        [1, 1, 0, 0, 0],  # Poorest
    ]


@pytest.fixture
def trait_matrix() -> List[List[float]]:
    """Species trait matrix (species x traits)."""
    return [
        [1.0, 2.0],
        [3.0, 4.0],
        [5.0, 6.0],
        [2.0, 3.0],
        [4.0, 5.0],
    ]


@pytest.fixture
def group_labels() -> List[str]:
    """Group labels for indicator/ANOSIM testing."""
    return ["A", "A", "B", "B"]


# ──────────────────────────────────────────────────────────────────────────────
# Community module tests
# ──────────────────────────────────────────────────────────────────────────────


class TestAlphaDiversity:
    """Tests for alpha diversity indices."""

    def test_shannon_equal_abundances(self, equal_abundances):
        from metainformant.ecology.analysis.community import shannon_diversity

        result = shannon_diversity(equal_abundances)
        expected = math.log(4)  # ln(S) for perfect evenness
        assert abs(result - expected) < 1e-10

    def test_shannon_single_species(self):
        from metainformant.ecology.analysis.community import shannon_diversity

        assert shannon_diversity([100]) == 0.0

    def test_shannon_empty(self):
        from metainformant.ecology.analysis.community import shannon_diversity

        assert shannon_diversity([]) == 0.0

    def test_shannon_increases_with_richness(self):
        from metainformant.ecology.analysis.community import shannon_diversity

        h2 = shannon_diversity([1, 1])
        h4 = shannon_diversity([1, 1, 1, 1])
        h8 = shannon_diversity([1, 1, 1, 1, 1, 1, 1, 1])
        assert h2 < h4 < h8

    def test_simpson_equal_abundances(self, equal_abundances):
        from metainformant.ecology.analysis.community import simpson_diversity

        result = simpson_diversity(equal_abundances)
        expected = 1.0 - 4 * (0.25**2)
        assert abs(result - expected) < 1e-10

    def test_simpson_range(self, simple_abundances):
        from metainformant.ecology.analysis.community import simpson_diversity

        result = simpson_diversity(simple_abundances)
        assert 0.0 <= result <= 1.0

    def test_simpson_single_species(self):
        from metainformant.ecology.analysis.community import simpson_diversity

        assert simpson_diversity([100]) == 0.0

    def test_species_richness_simple_list(self, simple_abundances):
        from metainformant.ecology.analysis.community import species_richness

        result = species_richness(simple_abundances)
        assert result == 8
        assert isinstance(result, int)

    def test_species_richness_with_zeros(self):
        from metainformant.ecology.analysis.community import species_richness

        assert species_richness([1, 0, 3, 0, 5]) == 3

    def test_species_richness_nested_list(self, community_matrix):
        from metainformant.ecology.analysis.community import species_richness

        result = species_richness(community_matrix)
        assert isinstance(result, list)
        assert len(result) == 4
        assert result[0] == 5  # Community A has 5 non-zero

    def test_species_richness_dict(self):
        from metainformant.ecology.analysis.community import species_richness

        data = {"siteA": [10, 5, 0], "siteB": [0, 0, 10]}
        result = species_richness(data)
        assert isinstance(result, dict)
        assert result["siteA"] == 2
        assert result["siteB"] == 1

    def test_species_richness_empty(self):
        from metainformant.ecology.analysis.community import species_richness

        assert species_richness([]) == 0

    def test_pielou_perfect_evenness(self, equal_abundances):
        from metainformant.ecology.analysis.community import pielou_evenness

        result = pielou_evenness(equal_abundances)
        assert abs(result - 1.0) < 1e-10

    def test_pielou_range(self, simple_abundances):
        from metainformant.ecology.analysis.community import pielou_evenness

        result = pielou_evenness(simple_abundances)
        assert 0.0 <= result <= 1.0

    def test_chao1_geq_observed(self, simple_abundances):
        from metainformant.ecology.analysis.community import chao1_estimator, species_richness_simple

        chao1 = chao1_estimator(simple_abundances)
        observed = species_richness_simple(simple_abundances)
        assert chao1 >= observed

    def test_chao1_no_singletons(self):
        from metainformant.ecology.analysis.community import chao1_estimator

        # No singletons = Chao1 equals observed
        assert chao1_estimator([10, 8, 6, 4]) == 4

    def test_chao1_all_singletons(self):
        from metainformant.ecology.analysis.community import chao1_estimator

        # f1=5, f2=0: S_obs + f1*(f1-1)/2 = 5 + 10 = 15
        assert chao1_estimator([1, 1, 1, 1, 1]) == 15.0

    def test_community_metrics_keys(self, simple_abundances):
        from metainformant.ecology.analysis.community import community_metrics

        metrics = community_metrics(simple_abundances)
        assert set(metrics.keys()) == {"shannon", "simpson", "richness", "pielou", "chao1"}

    def test_community_metrics_consistency(self, equal_abundances):
        from metainformant.ecology.analysis.community import community_metrics

        metrics = community_metrics(equal_abundances)
        assert abs(metrics["pielou"] - 1.0) < 1e-10
        assert metrics["richness"] == 4


class TestBetaDiversity:
    """Tests for beta diversity measures."""

    def test_bray_curtis_identical(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm = [10, 8, 6]
        assert beta_diversity(comm, comm, "bray_curtis") == 0.0

    def test_bray_curtis_different(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm1 = [10, 0, 0]
        comm2 = [0, 0, 10]
        result = beta_diversity(comm1, comm2, "bray_curtis")
        assert result == 1.0  # Completely different

    def test_bray_curtis_range(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm1 = [10, 8, 6, 4]
        comm2 = [8, 6, 4, 2]
        result = beta_diversity(comm1, comm2, "bray_curtis")
        assert 0.0 <= result <= 1.0

    def test_jaccard_identical(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm = [10, 8, 6]
        assert beta_diversity(comm, comm, "jaccard") == 0.0

    def test_jaccard_no_overlap(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm1 = [10, 0, 0]
        comm2 = [0, 10, 0]
        assert beta_diversity(comm1, comm2, "jaccard") == 1.0

    def test_sorensen_identical(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm = [5, 5, 5]
        assert beta_diversity(comm, comm, "sorensen") == 0.0

    def test_sorensen_range(self):
        from metainformant.ecology.analysis.community import beta_diversity

        comm1 = [10, 5, 0, 3]
        comm2 = [0, 5, 10, 0]
        result = beta_diversity(comm1, comm2, "sorensen")
        assert 0.0 <= result <= 1.0

    def test_invalid_method(self):
        from metainformant.ecology.analysis.community import beta_diversity

        with pytest.raises(ValueError, match="Unsupported"):
            beta_diversity([1], [1], "invalid")

    def test_similarity_matrix_symmetric(self, community_matrix):
        from metainformant.ecology.analysis.community import community_similarity_matrix

        sim = community_similarity_matrix(community_matrix)
        n = len(community_matrix)
        for i in range(n):
            for j in range(n):
                assert abs(sim[i][j] - sim[j][i]) < 1e-10

    def test_similarity_matrix_diagonal(self, community_matrix):
        from metainformant.ecology.analysis.community import community_similarity_matrix

        sim = community_similarity_matrix(community_matrix)
        for i in range(len(community_matrix)):
            assert abs(sim[i][i] - 1.0) < 1e-10

    def test_alpha_beta_gamma(self, community_matrix):
        from metainformant.ecology.analysis.community import alpha_beta_gamma_diversity

        result = alpha_beta_gamma_diversity(community_matrix)
        assert "alpha" in result
        assert "beta" in result
        assert "gamma" in result
        assert result["gamma"] >= result["alpha"]

    def test_alpha_beta_gamma_empty(self):
        from metainformant.ecology.analysis.community import alpha_beta_gamma_diversity

        result = alpha_beta_gamma_diversity([])
        assert result["alpha"] == 0.0


class TestCurves:
    """Tests for rarefaction, rank-abundance, and related curves."""

    def test_rarefaction_increasing(self, simple_abundances):
        from metainformant.ecology.analysis.community import rarefaction_curve

        curve = rarefaction_curve(simple_abundances, max_samples=20)
        richness_values = [r for _, r in curve]
        # Rarefaction curve should be monotonically non-decreasing
        for i in range(1, len(richness_values)):
            assert richness_values[i] >= richness_values[i - 1]

    def test_rarefaction_empty(self):
        from metainformant.ecology.analysis.community import rarefaction_curve

        curve = rarefaction_curve([])
        assert len(curve) == 1
        assert curve[0] == (0, 0.0)

    def test_rank_abundance_descending(self, simple_abundances):
        from metainformant.ecology.analysis.community import rank_abundance_curve

        curve = rank_abundance_curve(simple_abundances)
        abundances = [a for _, a in curve]
        for i in range(1, len(abundances)):
            assert abundances[i] <= abundances[i - 1]

    def test_rank_abundance_preserves_count(self, simple_abundances):
        from metainformant.ecology.analysis.community import rank_abundance_curve

        curve = rank_abundance_curve(simple_abundances)
        assert len(curve) == len(simple_abundances)

    def test_dominance_curve_reaches_one(self, simple_abundances):
        from metainformant.ecology.analysis.community import dominance_diversity_curve

        curve = dominance_diversity_curve(simple_abundances)
        last_dominance = curve[-1][0]
        assert abs(last_dominance - 1.0) < 1e-10

    def test_dominance_curve_empty(self):
        from metainformant.ecology.analysis.community import dominance_diversity_curve

        curve = dominance_diversity_curve([])
        assert curve == [(0.0, 0.0)]

    def test_species_accumulation(self):
        from metainformant.ecology.analysis.community import species_accumulation_curve

        effort = [5, 10, 14, 17, 19, 20]
        curve = species_accumulation_curve(effort)
        assert len(curve) == 6
        assert curve[0] == (1, 5.0)


class TestNestedness:
    """Tests for nestedness temperature."""

    def test_perfectly_nested(self, presence_absence):
        from metainformant.ecology.analysis.community import nestedness_temperature_calculator

        temp = nestedness_temperature_calculator(presence_absence)
        assert temp == 0.0  # Perfectly nested

    def test_random_matrix(self):
        from metainformant.ecology.analysis.community import nestedness_temperature_calculator

        # Anti-nested pattern
        matrix = [
            [1, 0, 0, 1, 0],
            [0, 1, 0, 0, 1],
            [0, 0, 1, 1, 0],
            [1, 0, 1, 0, 0],
        ]
        temp = nestedness_temperature_calculator(matrix)
        assert temp > 0.0  # Not perfectly nested

    def test_empty_matrix(self):
        from metainformant.ecology.analysis.community import nestedness_temperature_calculator

        assert nestedness_temperature_calculator([]) == 100.0


class TestSpeciesArea:
    """Tests for species-area relationship."""

    def test_power_law_fit(self):
        from metainformant.ecology.analysis.community import species_area_relationship

        areas = [1, 10, 100, 1000]
        # Perfect power law: S = 10 * A^0.3
        species = [int(10 * a**0.3) for a in areas]
        result = species_area_relationship(species, areas)
        assert result["model"] == "power_law"
        assert result["r_squared"] is not None
        assert result["r_squared"] > 0.9

    def test_too_few_points(self):
        from metainformant.ecology.analysis.community import species_area_relationship

        with pytest.raises(ValueError):
            species_area_relationship([10], [1])


class TestBiodiversityIndices:
    """Tests for calculate_biodiversity_indices."""

    def test_all_indices(self, community_matrix):
        from metainformant.ecology.analysis.community import calculate_biodiversity_indices

        result = calculate_biodiversity_indices(community_matrix)
        assert "shannon" in result
        assert "simpson" in result
        assert len(result["shannon"]) == len(community_matrix)

    def test_custom_indices(self, community_matrix):
        from metainformant.ecology.analysis.community import calculate_biodiversity_indices

        result = calculate_biodiversity_indices(community_matrix, indices=["shannon", "richness"])
        assert len(result) == 2
        assert "simpson" not in result


class TestEcologyReport:
    """Tests for report generation."""

    def test_report_generation(self, community_matrix, tmp_path):
        from metainformant.ecology.analysis.community import generate_ecology_report

        report = generate_ecology_report(community_matrix, output_path=tmp_path / "report.txt")
        assert "ECOLOGICAL COMMUNITY ANALYSIS REPORT" in report
        assert "Shannon" in report or "shannon" in report

    def test_report_to_file(self, community_matrix, tmp_path):
        from metainformant.ecology.analysis.community import generate_ecology_report

        output = tmp_path / "report.txt"
        generate_ecology_report(community_matrix, output_path=output)
        assert output.exists()
        content = output.read_text()
        assert len(content) > 0


class TestCalculateDiversity:
    """Tests for the main calculate_diversity dispatcher."""

    def test_dict_input(self):
        from metainformant.ecology.analysis.community import calculate_diversity

        data = {"siteA": [10, 5, 3], "siteB": [8, 8, 8]}
        result = calculate_diversity(data, "shannon")
        assert isinstance(result, dict)
        assert "siteA" in result

    def test_list_input(self, community_matrix):
        from metainformant.ecology.analysis.community import calculate_diversity

        result = calculate_diversity(community_matrix, "shannon")
        assert isinstance(result, list)
        assert len(result) == len(community_matrix)

    def test_invalid_method(self, community_matrix):
        from metainformant.ecology.analysis.community import calculate_diversity

        with pytest.raises(ValueError, match="Unsupported"):
            calculate_diversity(community_matrix, "invalid")

    def test_invsimpson(self, equal_abundances):
        from metainformant.ecology.analysis.community import calculate_single_diversity

        result = calculate_single_diversity(equal_abundances, "invsimpson")
        # Inverse Simpson for 4 equal species = 4.0
        assert abs(result - 4.0) < 1e-10


class TestEvenness:
    """Tests for calculate_evenness."""

    def test_pielou_method(self, equal_abundances):
        from metainformant.ecology.analysis.community import calculate_evenness

        result = calculate_evenness(equal_abundances, "pielou")
        assert abs(result - 1.0) < 1e-10

    def test_simpson_method(self, simple_abundances):
        from metainformant.ecology.analysis.community import calculate_evenness

        result = calculate_evenness(simple_abundances, "simpson")
        assert result > 0.0

    def test_invalid_method(self, simple_abundances):
        from metainformant.ecology.analysis.community import calculate_evenness

        with pytest.raises(ValueError, match="Unsupported"):
            calculate_evenness(simple_abundances, "invalid")


# ──────────────────────────────────────────────────────────────────────────────
# Ordination module tests
# ──────────────────────────────────────────────────────────────────────────────


class TestDistanceMatrix:
    """Tests for pairwise distance matrix computation."""

    def test_bray_curtis_distance(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix

        dist = distance_matrix(community_matrix, method="bray_curtis")
        n = len(community_matrix)
        assert len(dist) == n
        assert len(dist[0]) == n
        # Diagonal should be 0
        for i in range(n):
            assert abs(dist[i][i]) < 1e-10
        # Symmetric
        for i in range(n):
            for j in range(n):
                assert abs(dist[i][j] - dist[j][i]) < 1e-10

    def test_euclidean_distance(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix

        dist = distance_matrix(community_matrix, method="euclidean")
        assert dist[0][2] < 1e-10  # Community A == Community C

    def test_jaccard_distance(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix

        dist = distance_matrix(community_matrix, method="jaccard")
        for i in range(len(community_matrix)):
            assert abs(dist[i][i]) < 1e-10


class TestPCoA:
    """Tests for Principal Coordinates Analysis."""

    def test_pcoa_basic(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix, pcoa

        dist = distance_matrix(community_matrix)
        result = pcoa(dist, n_components=2)
        assert "coordinates" in result
        assert "eigenvalues" in result
        assert "variance_explained" in result
        assert len(result["coordinates"]) == 4
        assert len(result["coordinates"][0]) == 2

    def test_pcoa_variance_explained(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix, pcoa

        dist = distance_matrix(community_matrix)
        result = pcoa(dist, n_components=2)
        # Variance explained should be between 0 and 1
        for v in result["variance_explained"]:
            assert 0.0 <= v <= 1.0 + 1e-10

    def test_pcoa_identical_communities_close(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix, pcoa

        dist = distance_matrix(community_matrix)
        result = pcoa(dist, n_components=2)
        coords = result["coordinates"]
        # Community 0 and 2 are identical, should have same coordinates
        for dim in range(2):
            assert abs(coords[0][dim] - coords[2][dim]) < 1e-6


class TestNMDS:
    """Tests for Non-metric Multidimensional Scaling."""

    def test_nmds_basic(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix, nmds

        dist = distance_matrix(community_matrix)
        result = nmds(dist, n_components=2, max_iter=50, n_init=2)
        assert "coordinates" in result
        assert "stress" in result
        assert len(result["coordinates"]) == 4

    def test_nmds_stress_range(self, community_matrix):
        from metainformant.ecology.analysis.ordination import distance_matrix, nmds

        dist = distance_matrix(community_matrix)
        result = nmds(dist, n_components=2, max_iter=100, n_init=2)
        assert result["stress"] >= 0.0


class TestProcrustes:
    """Tests for Procrustes analysis."""

    def test_identical_coords(self):
        from metainformant.ecology.analysis.ordination import procrustes

        coords = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        result = procrustes(coords, coords)
        assert result["m2"] < 1e-6


# ──────────────────────────────────────────────────────────────────────────────
# Indicators module tests
# ──────────────────────────────────────────────────────────────────────────────


class TestIndVal:
    """Tests for Indicator Value analysis."""

    def test_indval_basic(self):
        from metainformant.ecology.analysis.indicators import indval

        # Species 1 is abundant only in group A, species 3 only in group B
        matrix = [
            [10, 0, 0, 5],  # Site 1 (group A)
            [8, 1, 0, 4],  # Site 2 (group A)
            [0, 5, 10, 3],  # Site 3 (group B)
            [0, 4, 8, 2],  # Site 4 (group B)
        ]
        groups = ["A", "A", "B", "B"]
        result = indval(matrix, groups)
        assert isinstance(result, list)
        assert len(result) > 0
        # Each result should have indval, specificity, fidelity
        for item in result:
            assert "indval" in item
            assert "best_group" in item

    def test_indval_perfect_indicator(self):
        from metainformant.ecology.analysis.indicators import indval

        # Species 0 only in group A, species 1 only in group B
        matrix = [
            [10, 0],
            [10, 0],
            [0, 10],
            [0, 10],
        ]
        groups = ["A", "A", "B", "B"]
        result = indval(matrix, groups)
        # Both species should have high indval (100)
        indvals = {r["species_idx"]: r["indval"] for r in result}
        assert indvals[0] == 100.0
        assert indvals[1] == 100.0


class TestANOSIM:
    """Tests for Analysis of Similarities."""

    def test_anosim_basic(self):
        from metainformant.ecology.analysis.indicators import anosim

        # Clear group structure
        dist = [
            [0.0, 0.1, 0.8, 0.9],
            [0.1, 0.0, 0.7, 0.8],
            [0.8, 0.7, 0.0, 0.1],
            [0.9, 0.8, 0.1, 0.0],
        ]
        groups = ["A", "A", "B", "B"]
        result = anosim(dist, groups, n_permutations=99)
        assert "r_statistic" in result
        assert "p_value" in result
        assert -1.0 <= result["r_statistic"] <= 1.0

    def test_anosim_no_structure(self):
        from metainformant.ecology.analysis.indicators import anosim

        # No group structure (all distances equal)
        dist = [
            [0.0, 0.5, 0.5, 0.5],
            [0.5, 0.0, 0.5, 0.5],
            [0.5, 0.5, 0.0, 0.5],
            [0.5, 0.5, 0.5, 0.0],
        ]
        groups = ["A", "A", "B", "B"]
        result = anosim(dist, groups, n_permutations=99)
        assert abs(result["r_statistic"]) < 0.3  # Should be near 0


class TestPERMANOVA:
    """Tests for Permutational MANOVA."""

    def test_permanova_basic(self):
        from metainformant.ecology.analysis.indicators import permanova

        dist = [
            [0.0, 0.1, 0.8, 0.9],
            [0.1, 0.0, 0.7, 0.8],
            [0.8, 0.7, 0.0, 0.1],
            [0.9, 0.8, 0.1, 0.0],
        ]
        groups = ["A", "A", "B", "B"]
        result = permanova(dist, groups, n_permutations=99)
        assert "f_statistic" in result
        assert "p_value" in result
        assert "r_squared" in result
        assert 0.0 <= result["r_squared"] <= 1.0


class TestClustering:
    """Tests for hierarchical community clustering."""

    def test_upgma_basic(self):
        from metainformant.ecology.analysis.indicators import cluster_communities

        dist = [
            [0.0, 0.2, 0.8, 0.9],
            [0.2, 0.0, 0.7, 0.8],
            [0.8, 0.7, 0.0, 0.1],
            [0.9, 0.8, 0.1, 0.0],
        ]
        result = cluster_communities(dist, method="upgma", n_clusters=2)
        assert "cluster_labels" in result
        assert "dendrogram" in result
        labels = result["cluster_labels"]
        # Sites 0,1 should be in same cluster, sites 2,3 in another
        assert labels[0] == labels[1]
        assert labels[2] == labels[3]
        assert labels[0] != labels[2]


class TestSIMPER:
    """Tests for Similarity Percentages."""

    def test_simper_basic(self):
        from metainformant.ecology.analysis.indicators import simper

        matrix = [
            [10, 0, 5],
            [8, 1, 6],
            [0, 10, 2],
            [1, 8, 3],
        ]
        groups = ["A", "A", "B", "B"]
        result = simper(matrix, groups)
        assert isinstance(result, list)
        assert len(result) > 0
        # Contributions should sum to ~100%
        total = sum(item["contribution_pct"] for item in result)
        assert abs(total - 100.0) < 1.0


# ──────────────────────────────────────────────────────────────────────────────
# Functional ecology module tests
# ──────────────────────────────────────────────────────────────────────────────


class TestFunctionalDiversity:
    """Tests for functional diversity indices."""

    def test_raos_q(self, trait_matrix, simple_abundances):
        from metainformant.ecology.analysis.functional import raos_quadratic_entropy

        abundances = simple_abundances[: len(trait_matrix)]
        result = raos_quadratic_entropy(trait_matrix, abundances)
        assert result >= 0.0

    def test_raos_q_single_species(self):
        from metainformant.ecology.analysis.functional import raos_quadratic_entropy

        traits = [[1.0, 2.0]]
        abundances = [10]
        result = raos_quadratic_entropy(traits, abundances)
        assert result == 0.0  # Single species, no functional diversity

    def test_cwm_basic(self, trait_matrix):
        from metainformant.ecology.analysis.functional import community_weighted_mean

        abundances = [1, 1, 1, 1, 1]  # Equal weights
        result = community_weighted_mean(trait_matrix, abundances)
        assert len(result) == 2  # Two traits
        # Mean of [1,3,5,2,4] = 3.0, mean of [2,4,6,3,5] = 4.0
        assert abs(result[0] - 3.0) < 1e-10
        assert abs(result[1] - 4.0) < 1e-10

    def test_functional_dispersion(self, trait_matrix):
        from metainformant.ecology.analysis.functional import functional_dispersion

        abundances = [10, 8, 6, 4, 2]
        result = functional_dispersion(trait_matrix, abundances)
        assert result >= 0.0

    def test_functional_redundancy(self, trait_matrix):
        from metainformant.ecology.analysis.functional import functional_redundancy

        abundances = [10, 8, 6, 4, 2]
        result = functional_redundancy(trait_matrix, abundances)
        # Redundancy can be negative (functional diversity > species diversity)
        assert isinstance(result, float)

    def test_functional_richness(self):
        from metainformant.ecology.analysis.functional import functional_richness

        # Use non-collinear points so convex hull has non-zero area
        non_collinear = [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.5, 0.5],
        ]
        result = functional_richness(non_collinear)
        assert result > 0.0

    def test_functional_richness_single(self):
        from metainformant.ecology.analysis.functional import functional_richness

        result = functional_richness([[1.0, 2.0]])
        assert result == 0.0  # Single point has no volume

    def test_trait_distance_euclidean(self, trait_matrix):
        from metainformant.ecology.analysis.functional import trait_distance_matrix

        dist = trait_distance_matrix(trait_matrix, method="euclidean")
        assert len(dist) == len(trait_matrix)
        # Diagonal should be 0
        for i in range(len(trait_matrix)):
            assert abs(dist[i][i]) < 1e-10
        # Symmetric
        for i in range(len(trait_matrix)):
            for j in range(len(trait_matrix)):
                assert abs(dist[i][j] - dist[j][i]) < 1e-10


# ──────────────────────────────────────────────────────────────────────────────
# Macroecology module tests
# ──────────────────────────────────────────────────────────────────────────────


class TestSADModels:
    """Tests for Species Abundance Distribution models."""

    def test_logseries_fit(self):
        from metainformant.ecology.analysis.macroecology import fit_logseries

        # Typical SAD data
        abundances = [100, 50, 25, 12, 6, 3, 2, 1, 1, 1]
        result = fit_logseries(abundances)
        assert "alpha" in result
        assert result["alpha"] > 0

    def test_lognormal_fit(self):
        from metainformant.ecology.analysis.macroecology import fit_lognormal

        abundances = [100, 80, 60, 40, 20, 10, 5, 2, 1, 1]
        result = fit_lognormal(abundances)
        assert "mu" in result
        assert "sigma" in result

    def test_broken_stick(self):
        from metainformant.ecology.analysis.macroecology import fit_broken_stick

        abundances = [50, 30, 15, 5]
        result = fit_broken_stick(abundances)
        assert "expected_abundances" in result
        assert len(result["expected_abundances"]) == len(abundances)

    def test_geometric_series(self):
        from metainformant.ecology.analysis.macroecology import fit_geometric_series

        abundances = [100, 50, 25, 12, 6, 3]
        result = fit_geometric_series(abundances)
        assert "k" in result
        assert 0.0 < result["k"] < 1.0

    def test_compare_sad_models(self):
        from metainformant.ecology.analysis.macroecology import compare_sad_models

        abundances = [100, 80, 60, 40, 20, 10, 5, 2, 1, 1]
        result = compare_sad_models(abundances)
        assert isinstance(result, dict)
        assert len(result) >= 3  # At least 3 models compared


class TestSAR:
    """Tests for Species-Area Relationships."""

    def test_power_law_sar(self):
        from metainformant.ecology.analysis.macroecology import species_area_power

        areas = [1, 10, 100, 1000, 10000]
        species = [10, 20, 40, 80, 160]
        result = species_area_power(areas, species)
        assert "c" in result
        assert "z" in result
        assert "r_squared" in result
        assert result["r_squared"] > 0.5

    def test_logarithmic_sar(self):
        from metainformant.ecology.analysis.macroecology import species_area_logarithmic

        areas = [1, 10, 100, 1000]
        species = [5, 10, 15, 20]
        result = species_area_logarithmic(areas, species)
        assert "a" in result
        assert "b" in result


class TestDistanceDecay:
    """Tests for distance-decay of similarity."""

    def test_distance_decay_basic(self):
        from metainformant.ecology.analysis.macroecology import distance_decay

        distances = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        similarities = [0.9, 0.7, 0.5, 0.3, 0.15, 0.05]
        result = distance_decay(distances, similarities)
        assert "exponential" in result or "power_law" in result
        assert "best_model" in result


class TestOccupancyFrequency:
    """Tests for occupancy-frequency distribution."""

    def test_occupancy_basic(self, presence_absence):
        from metainformant.ecology.analysis.macroecology import occupancy_frequency

        result = occupancy_frequency(presence_absence)
        assert "core_species" in result
        assert "satellite_species" in result
        total = result["core_species"] + result["common_species"] + result["satellite_species"]
        assert total == len(presence_absence[0])


class TestTaylorsLaw:
    """Tests for Taylor's Power Law."""

    def test_taylors_basic(self):
        from metainformant.ecology.analysis.macroecology import taylors_power_law

        means = [1.0, 5.0, 10.0, 50.0, 100.0]
        variances = [1.0, 25.0, 100.0, 2500.0, 10000.0]
        result = taylors_power_law(means, variances)
        assert "a" in result
        assert "b" in result
        # b should be approximately 2 for this quadratic relationship
        assert abs(result["b"] - 2.0) < 0.5


class TestMetabolicScaling:
    """Tests for Metabolic Theory of Ecology."""

    def test_metabolic_scaling_basic(self):
        from metainformant.ecology.analysis.macroecology import metabolic_scaling

        # Kleiber's law: B = B0 * M^0.75
        masses = [1, 10, 100, 1000, 10000]
        rates = [m**0.75 for m in masses]
        result = metabolic_scaling(masses, rates)
        assert "b_exponent" in result
        assert abs(result["b_exponent"] - 0.75) < 0.1


# ──────────────────────────────────────────────────────────────────────────────
# Package-level import tests
# ──────────────────────────────────────────────────────────────────────────────


class TestPackageImports:
    """Test that all modules and functions are importable."""

    def test_import_ecology(self):
        import metainformant.ecology

    def test_import_analysis(self):
        from metainformant.ecology import analysis

    def test_import_community(self):
        from metainformant.ecology.analysis import community

    def test_import_ordination(self):
        from metainformant.ecology.analysis import ordination

    def test_import_indicators(self):
        from metainformant.ecology.analysis import indicators

    def test_import_functional(self):
        from metainformant.ecology.analysis import functional

    def test_import_macroecology(self):
        from metainformant.ecology.analysis import macroecology

    def test_import_visualization(self):
        from metainformant.ecology import visualization

    def test_top_level_re_exports(self):
        """Test that key functions are available at the ecology package level."""
        from metainformant.ecology import (
            beta_diversity,
            community_metrics,
            shannon_diversity,
            simpson_diversity,
            species_richness,
        )
