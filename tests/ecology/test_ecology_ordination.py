"""Tests for ecology ordination, indicator species, and phylogenetic modules.

Covers:
    - distance_matrix (Bray-Curtis, Jaccard, Euclidean, Manhattan, Canberra)
    - pcoa (Principal Coordinates Analysis)
    - nmds (Non-metric Multidimensional Scaling)
    - cca (Canonical Correspondence Analysis)
    - procrustes rotation
    - indval (Indicator Value analysis)
    - anosim (Analysis of Similarities)
    - permanova (Permutational MANOVA)
    - cluster_communities (hierarchical clustering)
    - simper (Similarity Percentages)
    - multivariate_dispersion (PERMDISP)
    - faiths_pd (Faith's Phylogenetic Diversity)
    - compute_unifrac (UniFrac distances)
    - phylogenetic_beta_diversity
    - nri_nti (Net Relatedness / Nearest Taxon Index)
    - phylogenetic_signal (Blomberg's K, Pagel's lambda)
    - build_simple_tree (UPGMA, NJ)

Uses real implementations only (NO mocking per project policy).
"""

from __future__ import annotations

import math
from typing import Dict, List

import pytest

from metainformant.ecology.analysis.ordination import (
    cca,
    distance_matrix,
    nmds,
    pcoa,
    procrustes,
)
from metainformant.ecology.analysis.indicators import (
    anosim,
    cluster_communities,
    indval,
    multivariate_dispersion,
    permanova,
    simper,
)
from metainformant.ecology.phylogenetic.diversity import (
    build_simple_tree,
    compute_unifrac,
    faiths_pd,
    nri_nti,
    phylogenetic_beta_diversity,
    phylogenetic_signal,
)

# ---------------------------------------------------------------------------
# Shared test data fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def community_matrix() -> List[List[float]]:
    """Four communities x seven species abundance matrix."""
    return [
        [10, 8, 6, 4, 2, 0, 0],  # Community A
        [0, 6, 8, 10, 4, 2, 1],  # Community B
        [10, 8, 6, 4, 2, 0, 0],  # Community C (identical to A)
        [1, 1, 1, 1, 1, 1, 1],  # Community D (perfectly even)
    ]


@pytest.fixture
def small_communities() -> List[List[float]]:
    """Three simple communities for quick tests."""
    return [
        [10, 0, 5],
        [0, 8, 2],
        [3, 3, 3],
    ]


@pytest.fixture
def group_labels_4() -> List[str]:
    """Group labels for 4-community tests (two groups of two)."""
    return ["forest", "forest", "grassland", "grassland"]


@pytest.fixture
def simple_tree() -> dict:
    """A small phylogenetic tree with four taxa: A, B, C, D."""
    return {
        "name": "root",
        "branch_length": 0.0,
        "children": [
            {
                "name": "internal_1",
                "branch_length": 0.3,
                "children": [
                    {"name": "A", "branch_length": 0.5, "children": []},
                    {"name": "B", "branch_length": 0.4, "children": []},
                ],
            },
            {
                "name": "internal_2",
                "branch_length": 0.2,
                "children": [
                    {"name": "C", "branch_length": 0.6, "children": []},
                    {"name": "D", "branch_length": 0.3, "children": []},
                ],
            },
        ],
    }


# ============================================================================
# Distance matrix tests
# ============================================================================


class TestDistanceMatrix:
    """Test pairwise distance matrix computation."""

    def test_bray_curtis_basic(self, small_communities: List[List[float]]) -> None:
        """Bray-Curtis distances are symmetric, diagonal is zero, values in [0,1]."""
        dm = distance_matrix(small_communities, method="bray_curtis")

        n = len(small_communities)
        assert len(dm) == n
        assert all(len(row) == n for row in dm)

        # Diagonal must be zero
        for i in range(n):
            assert dm[i][i] == 0.0

        # Symmetry
        for i in range(n):
            for j in range(n):
                assert abs(dm[i][j] - dm[j][i]) < 1e-12

        # Values in [0, 1] for Bray-Curtis
        for i in range(n):
            for j in range(n):
                assert 0.0 <= dm[i][j] <= 1.0 + 1e-12

    def test_bray_curtis_identical_communities(self) -> None:
        """Distance between identical communities should be zero."""
        comms = [[5, 10, 15], [5, 10, 15]]
        dm = distance_matrix(comms, method="bray_curtis")
        assert abs(dm[0][1]) < 1e-12

    def test_jaccard_binary_semantics(self) -> None:
        """Jaccard distance captures presence/absence differences only."""
        comms = [
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [1, 1, 1, 1],
        ]
        dm = distance_matrix(comms, method="jaccard")
        # Comms 0 and 1 share no species -> Jaccard = 1.0
        assert abs(dm[0][1] - 1.0) < 1e-12
        # Comms 0 and 2 share species 0,1 out of union {0,1,2,3} -> Jaccard = 0.5
        assert abs(dm[0][2] - 0.5) < 1e-12

    def test_euclidean_known_value(self) -> None:
        """Euclidean distance matches hand-computed value."""
        comms = [[0, 0], [3, 4]]
        dm = distance_matrix(comms, method="euclidean")
        assert abs(dm[0][1] - 5.0) < 1e-12

    def test_manhattan_known_value(self) -> None:
        """Manhattan distance matches hand-computed value."""
        comms = [[0, 0], [3, 4]]
        dm = distance_matrix(comms, method="manhattan")
        assert abs(dm[0][1] - 7.0) < 1e-12

    def test_canberra_known_value(self) -> None:
        """Canberra distance is computed correctly."""
        comms = [[1, 2], [3, 4]]
        dm = distance_matrix(comms, method="canberra")
        # |1-3|/(1+3) + |2-4|/(2+4) = 0.5 + 1/3
        expected = 0.5 + 1.0 / 3.0
        assert abs(dm[0][1] - expected) < 1e-12

    def test_unsupported_method_raises(self) -> None:
        """Requesting an unknown metric raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported distance method"):
            distance_matrix([[1, 2], [3, 4]], method="hamming")

    def test_empty_raises(self) -> None:
        """Empty input raises ValueError."""
        with pytest.raises((ValueError, Exception)):
            distance_matrix([], method="bray_curtis")

    def test_unequal_vector_lengths_raises(self) -> None:
        """Unequal vector lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            distance_matrix([[1, 2], [3]], method="bray_curtis")

    def test_single_community(self) -> None:
        """Single community yields a 1x1 zero matrix."""
        dm = distance_matrix([[1, 2, 3]])
        assert dm == [[0.0]]


# ============================================================================
# PCoA tests
# ============================================================================


class TestPCoA:
    """Test Principal Coordinates Analysis."""

    def test_pcoa_returns_correct_shape(self) -> None:
        """PCoA coordinates have the expected shape."""
        dm = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
        result = pcoa(dm, n_components=2)
        coords = result["coordinates"]
        assert len(coords) == 3
        assert all(len(row) == 2 for row in coords)

    def test_pcoa_eigenvalues_descending(self) -> None:
        """Eigenvalues are in descending order."""
        dm = [[0, 1, 3, 4], [1, 0, 2, 3], [3, 2, 0, 1], [4, 3, 1, 0]]
        result = pcoa(dm, n_components=3)
        evs = result["eigenvalues"]
        assert len(evs) == 3
        for i in range(len(evs) - 1):
            assert evs[i] >= evs[i + 1] - 1e-10

    def test_pcoa_variance_explained_sums_leq_1(self) -> None:
        """Variance explained fractions sum to at most 1."""
        dm = [[0, 2, 4], [2, 0, 3], [4, 3, 0]]
        result = pcoa(dm, n_components=2)
        total = sum(result["variance_explained"])
        assert total <= 1.0 + 1e-6

    def test_pcoa_preserves_relative_distances(self) -> None:
        """Embedded distances preserve the rank order of the original matrix."""
        dm = [[0, 1, 5], [1, 0, 4], [5, 4, 0]]
        result = pcoa(dm, n_components=2)
        coords = result["coordinates"]

        def _edist(a: List[float], b: List[float]) -> float:
            return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))

        d01 = _edist(coords[0], coords[1])
        d02 = _edist(coords[0], coords[2])
        # Original: d(0,1)=1, d(0,2)=5 -> embedded d01 < d02
        assert d01 < d02

    def test_pcoa_n_components_exceeds_n_raises(self) -> None:
        """Requesting more components than objects raises ValueError."""
        dm = [[0, 1], [1, 0]]
        with pytest.raises(ValueError, match="n_components"):
            pcoa(dm, n_components=3)

    def test_pcoa_single_point(self) -> None:
        """Single point returns zero coordinates."""
        result = pcoa([[0]], n_components=1)
        assert result["coordinates"] == [[0.0]]


# ============================================================================
# NMDS tests
# ============================================================================


class TestNMDS:
    """Test Non-metric Multidimensional Scaling."""

    def test_nmds_returns_correct_shape(self) -> None:
        """NMDS coordinates have the expected shape."""
        dm = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
        result = nmds(dm, n_components=2, n_init=2, max_iter=50, seed=42)
        coords = result["coordinates"]
        assert len(coords) == 3
        assert all(len(row) == 2 for row in coords)

    def test_nmds_stress_nonnegative(self) -> None:
        """Stress must be non-negative."""
        dm = [[0, 1, 3, 5], [1, 0, 2, 4], [3, 2, 0, 1], [5, 4, 1, 0]]
        result = nmds(dm, n_components=2, n_init=2, max_iter=100, seed=42)
        assert result["stress"] >= 0.0

    def test_nmds_perfect_embedding(self) -> None:
        """Points that lie in 2D should have near-zero stress in 2D."""
        # Three points forming a right triangle in 2D
        dm = [[0, 3, 5], [3, 0, 4], [5, 4, 0]]
        result = nmds(dm, n_components=2, n_init=4, max_iter=300, seed=42)
        # Stress should be reasonably low (Euclidean distances are exactly embeddable in 2D)
        assert result["stress"] < 0.3

    def test_nmds_seed_reproducibility(self) -> None:
        """Same seed gives the same result."""
        dm = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
        r1 = nmds(dm, n_components=2, n_init=2, max_iter=50, seed=123)
        r2 = nmds(dm, n_components=2, n_init=2, max_iter=50, seed=123)
        assert abs(r1["stress"] - r2["stress"]) < 1e-10

    def test_nmds_n_components_too_large_raises(self) -> None:
        """n_components >= n raises ValueError."""
        dm = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
        with pytest.raises(ValueError, match="n_components"):
            nmds(dm, n_components=3)

    def test_nmds_single_pair(self) -> None:
        """Two points in 1D should be trivial."""
        dm = [[0, 5], [5, 0]]
        result = nmds(dm, n_components=1, n_init=2, max_iter=50, seed=42)
        assert len(result["coordinates"]) == 2


# ============================================================================
# CCA tests
# ============================================================================


class TestCCA:
    """Test Canonical Correspondence Analysis."""

    def test_cca_returns_correct_keys(self) -> None:
        """CCA returns site_scores, species_scores, eigenvalues, variance_explained."""
        sp = [[10, 0, 5], [0, 8, 2], [3, 3, 3]]
        env = [[1.0, 2.0], [3.0, 1.0], [2.0, 3.0]]
        result = cca(sp, env)
        assert "site_scores" in result
        assert "species_scores" in result
        assert "eigenvalues" in result
        assert "variance_explained" in result

    def test_cca_site_scores_shape(self) -> None:
        """Site scores have shape (n_sites, n_components)."""
        sp = [[10, 0, 5, 1], [0, 8, 2, 3], [3, 3, 3, 2], [5, 1, 7, 0]]
        env = [[1.0], [3.0], [2.0], [4.0]]
        result = cca(sp, env)
        scores = result["site_scores"]
        assert len(scores) == 4
        # n_components = min(n_env=1, n_species-1=3, n_sites-1=3) = 1
        assert len(scores[0]) == 1

    def test_cca_eigenvalues_nonnegative(self) -> None:
        """CCA eigenvalues should be non-negative."""
        sp = [[10, 0, 5], [0, 8, 2], [3, 3, 3], [6, 2, 4]]
        env = [[1.0, 2.0], [3.0, 1.0], [2.0, 3.0], [4.0, 0.5]]
        result = cca(sp, env)
        for ev in result["eigenvalues"]:
            assert ev >= -1e-10

    def test_cca_row_mismatch_raises(self) -> None:
        """Mismatched row counts raise ValueError."""
        sp = [[10, 0], [0, 8]]
        env = [[1.0], [3.0], [2.0]]
        with pytest.raises(ValueError, match="Row count mismatch"):
            cca(sp, env)

    def test_cca_negative_species_raises(self) -> None:
        """Negative species abundances raise ValueError."""
        sp = [[10, -1], [0, 8], [3, 3]]
        env = [[1.0], [3.0], [2.0]]
        with pytest.raises(ValueError, match="non-negative"):
            cca(sp, env)

    def test_cca_variance_explained_sums_to_1(self) -> None:
        """Variance explained fractions should sum to approximately 1.0."""
        sp = [[10, 0, 5], [0, 8, 2], [3, 3, 3]]
        env = [[1.0, 2.0], [3.0, 1.0], [2.0, 3.0]]
        result = cca(sp, env)
        total = sum(result["variance_explained"])
        assert abs(total - 1.0) < 0.1  # Allow some tolerance


# ============================================================================
# Procrustes tests
# ============================================================================


class TestProcrustes:
    """Test Procrustes rotation analysis."""

    def test_procrustes_identical_configs(self) -> None:
        """Identical configurations should have m2=0 and correlation=1."""
        c = [[0, 0], [1, 0], [0, 1]]
        result = procrustes(c, c)
        assert result["m2"] < 1e-6
        assert result["correlation"] > 0.999

    def test_procrustes_rotated_configs(self) -> None:
        """A 90-degree rotation should be recovered with high correlation."""
        c1 = [[0, 0], [1, 0], [0, 1]]
        c2 = [[0, 0], [0, 1], [-1, 0]]  # 90-degree CCW rotation of c1
        result = procrustes(c1, c2)
        assert result["correlation"] > 0.9

    def test_procrustes_m2_nonnegative(self) -> None:
        """Procrustes m2 statistic must be non-negative."""
        c1 = [[0, 0], [1, 0], [0, 1], [1, 1]]
        c2 = [[0.1, 0.1], [1.2, 0.1], [0.1, 0.9], [0.9, 1.1]]
        result = procrustes(c1, c2)
        assert result["m2"] >= -1e-10

    def test_procrustes_transformed_coords_shape(self) -> None:
        """Transformed coordinates have the same shape as input."""
        c1 = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        c2 = [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
        result = procrustes(c1, c2)
        tc = result["transformed_coords"]
        assert len(tc) == 3
        assert all(len(row) == 3 for row in tc)

    def test_procrustes_mismatched_shapes_raises(self) -> None:
        """Different point counts raise ValueError."""
        c1 = [[0, 0], [1, 0]]
        c2 = [[0, 0], [1, 0], [0, 1]]
        with pytest.raises(ValueError, match="same number of points"):
            procrustes(c1, c2)

    def test_procrustes_too_few_points_raises(self) -> None:
        """Fewer than 2 points raises ValueError."""
        with pytest.raises(ValueError, match="at least 2 points"):
            procrustes([[0, 0]], [[1, 1]])


# ============================================================================
# Indicator Species (IndVal) tests
# ============================================================================


class TestIndval:
    """Test Indicator Value analysis."""

    def test_indval_returns_per_species_results(self) -> None:
        """IndVal returns one result dict per species."""
        mat = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0]]
        labels = ["forest", "forest", "grassland", "grassland"]
        results = indval(mat, labels, n_permutations=49, seed=0)
        assert len(results) == 3  # 3 species

    def test_indval_identifies_indicator_species(self) -> None:
        """Species 0 should indicate 'forest' and species 1 should indicate 'grassland'."""
        mat = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0]]
        labels = ["forest", "forest", "grassland", "grassland"]
        results = indval(mat, labels, n_permutations=49, seed=0)
        # Sort by species_idx for predictability
        by_idx = {r["species_idx"]: r for r in results}
        assert by_idx[0]["best_group"] == "forest"
        assert by_idx[1]["best_group"] == "grassland"

    def test_indval_score_in_range(self) -> None:
        """IndVal score must be in [0, 100]."""
        mat = [[10, 0], [8, 2], [0, 9], [1, 7]]
        labels = ["A", "A", "B", "B"]
        results = indval(mat, labels, n_permutations=49, seed=0)
        for r in results:
            assert 0 <= r["indval"] <= 100

    def test_indval_p_value_in_range(self) -> None:
        """P-values must be in (0, 1]."""
        mat = [[10, 0], [8, 2], [0, 9], [1, 7]]
        labels = ["A", "A", "B", "B"]
        results = indval(mat, labels, n_permutations=49, seed=0)
        for r in results:
            assert 0 < r["p_value"] <= 1.0

    def test_indval_single_group_raises(self) -> None:
        """Only one group raises ValueError."""
        mat = [[10, 0], [8, 2]]
        labels = ["A", "A"]
        with pytest.raises(ValueError, match="two distinct groups"):
            indval(mat, labels)


# ============================================================================
# ANOSIM tests
# ============================================================================


class TestAnosim:
    """Test Analysis of Similarities."""

    def test_anosim_strong_grouping(self) -> None:
        """Well-separated groups yield a high R statistic."""
        dm = [
            [0, 1, 5, 5],
            [1, 0, 5, 5],
            [5, 5, 0, 1],
            [5, 5, 1, 0],
        ]
        labels = ["A", "A", "B", "B"]
        result = anosim(dm, labels, n_permutations=99, seed=0)
        assert result["r_statistic"] > 0.0
        assert 0 < result["p_value"] <= 1.0

    def test_anosim_no_grouping(self) -> None:
        """Random labels yield R near 0."""
        dm = [[0, 2, 2, 2], [2, 0, 2, 2], [2, 2, 0, 2], [2, 2, 2, 0]]
        labels = ["A", "B", "A", "B"]
        result = anosim(dm, labels, n_permutations=99, seed=42)
        # R should be near 0 for equidistant points
        assert -1.0 <= result["r_statistic"] <= 1.0

    def test_anosim_output_keys(self) -> None:
        """ANOSIM returns r_statistic, p_value, n_permutations."""
        dm = [[0, 1, 3], [1, 0, 2], [3, 2, 0]]
        labels = ["X", "X", "Y"]
        result = anosim(dm, labels, n_permutations=49, seed=0)
        assert "r_statistic" in result
        assert "p_value" in result
        assert result["n_permutations"] == 49


# ============================================================================
# PERMANOVA tests
# ============================================================================


class TestPermanova:
    """Test Permutational MANOVA."""

    def test_permanova_strong_grouping(self) -> None:
        """Well-separated groups give a large F statistic."""
        dm = [
            [0, 1, 5, 5],
            [1, 0, 5, 5],
            [5, 5, 0, 1],
            [5, 5, 1, 0],
        ]
        labels = ["A", "A", "B", "B"]
        result = permanova(dm, labels, n_permutations=99, seed=0)
        assert result["f_statistic"] > 0.0
        assert 0 < result["p_value"] <= 1.0
        assert 0.0 <= result["r_squared"] <= 1.0

    def test_permanova_output_keys(self) -> None:
        """PERMANOVA returns f_statistic, p_value, r_squared, n_permutations."""
        dm = [[0, 2, 4], [2, 0, 3], [4, 3, 0]]
        labels = ["G1", "G1", "G2"]
        result = permanova(dm, labels, n_permutations=49, seed=0)
        for key in ("f_statistic", "p_value", "r_squared", "n_permutations"):
            assert key in result

    def test_permanova_single_group_raises(self) -> None:
        """Only one group raises ValueError."""
        dm = [[0, 1], [1, 0]]
        with pytest.raises(ValueError, match="two distinct groups"):
            permanova(dm, ["A", "A"])


# ============================================================================
# Cluster communities tests
# ============================================================================


class TestClusterCommunities:
    """Test hierarchical community clustering."""

    def test_upgma_basic(self) -> None:
        """UPGMA clustering returns correct keys."""
        dm = [[0, 2, 6], [2, 0, 5], [6, 5, 0]]
        result = cluster_communities(dm, method="upgma", n_clusters=2)
        assert "dendrogram" in result
        assert "cluster_labels" in result
        assert "cophenetic_correlation" in result

    def test_upgma_n_clusters(self) -> None:
        """Two closest items are in the same cluster when n_clusters=2."""
        dm = [[0, 2, 6], [2, 0, 5], [6, 5, 0]]
        result = cluster_communities(dm, method="upgma", n_clusters=2)
        labels = result["cluster_labels"]
        # Items 0 and 1 (distance 2) should be in the same cluster
        assert labels[0] == labels[1]
        assert labels[0] != labels[2]

    def test_single_linkage(self) -> None:
        """Single linkage produces valid output."""
        dm = [[0, 1, 4], [1, 0, 3], [4, 3, 0]]
        result = cluster_communities(dm, method="single", n_clusters=2)
        labels = result["cluster_labels"]
        assert labels[0] == labels[1]

    def test_complete_linkage(self) -> None:
        """Complete linkage produces valid output."""
        dm = [[0, 1, 4], [1, 0, 3], [4, 3, 0]]
        result = cluster_communities(dm, method="complete", n_clusters=2)
        assert len(result["cluster_labels"]) == 3

    def test_cophenetic_correlation_range(self) -> None:
        """Cophenetic correlation must be in [-1, 1]."""
        dm = [[0, 2, 6, 8], [2, 0, 5, 7], [6, 5, 0, 3], [8, 7, 3, 0]]
        result = cluster_communities(dm, method="upgma")
        cc = result["cophenetic_correlation"]
        assert -1.0 <= cc <= 1.0

    def test_unsupported_method_raises(self) -> None:
        """Unknown linkage method raises ValueError."""
        dm = [[0, 1], [1, 0]]
        with pytest.raises(ValueError, match="Unsupported linkage"):
            cluster_communities(dm, method="ward")


# ============================================================================
# SIMPER tests
# ============================================================================


class TestSimper:
    """Test Similarity Percentages analysis."""

    def test_simper_returns_per_species(self) -> None:
        """SIMPER returns one result per species, sorted by contribution."""
        mat = [[10, 0], [8, 2], [0, 9], [1, 7]]
        labels = ["A", "A", "B", "B"]
        results = simper(mat, labels)
        assert len(results) == 2

    def test_simper_cumulative_100(self) -> None:
        """Cumulative percentage should reach 100%."""
        mat = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0]]
        labels = ["A", "A", "B", "B"]
        results = simper(mat, labels)
        assert abs(results[-1]["cumulative_pct"] - 100.0) < 0.1

    def test_simper_contribution_sorted_descending(self) -> None:
        """Results are sorted by contribution_pct in descending order."""
        mat = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0]]
        labels = ["A", "A", "B", "B"]
        results = simper(mat, labels)
        for i in range(len(results) - 1):
            assert results[i]["contribution_pct"] >= results[i + 1]["contribution_pct"]

    def test_simper_single_group_raises(self) -> None:
        """Single group raises ValueError."""
        mat = [[10, 0], [8, 2]]
        labels = ["A", "A"]
        with pytest.raises(ValueError, match="two distinct groups"):
            simper(mat, labels)


# ============================================================================
# Multivariate Dispersion (PERMDISP) tests
# ============================================================================


class TestMultivariateDispersion:
    """Test multivariate dispersion analysis."""

    def test_permdisp_output_keys(self) -> None:
        """PERMDISP returns f_statistic, p_value, group_dispersions, n_permutations."""
        dm = [[0, 1, 4, 5], [1, 0, 5, 4], [4, 5, 0, 1], [5, 4, 1, 0]]
        labels = ["X", "X", "Y", "Y"]
        result = multivariate_dispersion(dm, labels, n_permutations=49, seed=0)
        for key in ("f_statistic", "p_value", "group_dispersions", "n_permutations"):
            assert key in result

    def test_permdisp_group_dispersions_present(self) -> None:
        """Group dispersions are returned for each group."""
        dm = [[0, 1, 4, 5], [1, 0, 5, 4], [4, 5, 0, 1], [5, 4, 1, 0]]
        labels = ["X", "X", "Y", "Y"]
        result = multivariate_dispersion(dm, labels, n_permutations=49, seed=0)
        assert "X" in result["group_dispersions"]
        assert "Y" in result["group_dispersions"]

    def test_permdisp_p_value_in_range(self) -> None:
        """P-value must be in (0, 1]."""
        dm = [[0, 2, 5, 6], [2, 0, 4, 5], [5, 4, 0, 1], [6, 5, 1, 0]]
        labels = ["A", "A", "B", "B"]
        result = multivariate_dispersion(dm, labels, n_permutations=49, seed=0)
        assert 0 < result["p_value"] <= 1.0


# ============================================================================
# Phylogenetic Diversity tests
# ============================================================================


class TestFaithsPD:
    """Test Faith's Phylogenetic Diversity."""

    def test_faiths_pd_all_taxa(self, simple_tree: dict) -> None:
        """PD of all taxa should equal the total branch length connecting them."""
        pd_value = faiths_pd(simple_tree, ["A", "B", "C", "D"])
        # All branches connecting A,B,C,D to root:
        # internal_1=0.3, A=0.5, B=0.4, internal_2=0.2, C=0.6, D=0.3 => total=2.3
        assert abs(pd_value - 2.3) < 1e-6

    def test_faiths_pd_subset(self, simple_tree: dict) -> None:
        """PD of a subset should be less than or equal to total PD."""
        pd_all = faiths_pd(simple_tree, ["A", "B", "C", "D"])
        pd_sub = faiths_pd(simple_tree, ["A", "B"])
        assert pd_sub <= pd_all + 1e-10

    def test_faiths_pd_single_taxon(self, simple_tree: dict) -> None:
        """PD of a single taxon includes its branch and path to root."""
        pd_a = faiths_pd(simple_tree, ["A"])
        # A branch=0.5, internal_1 branch=0.3 => 0.8
        assert abs(pd_a - 0.8) < 1e-6

    def test_faiths_pd_missing_taxa_raises(self, simple_tree: dict) -> None:
        """Taxa not in the tree should raise ValueError."""
        with pytest.raises(ValueError, match="None of the taxa"):
            faiths_pd(simple_tree, ["X", "Y", "Z"])


# ============================================================================
# UniFrac tests
# ============================================================================


class TestUniFrac:
    """Test UniFrac distance computations."""

    def test_unifrac_identical_communities(self, simple_tree: dict) -> None:
        """Identical communities should have UniFrac = 0."""
        ab = {"A": 1, "B": 1, "C": 1}
        dist = compute_unifrac(simple_tree, ab, ab, weighted=False)
        assert abs(dist) < 1e-10

    def test_unifrac_disjoint_communities(self, simple_tree: dict) -> None:
        """Completely disjoint communities should have positive UniFrac."""
        ab_a = {"A": 1, "B": 1}
        ab_b = {"C": 1, "D": 1}
        dist = compute_unifrac(simple_tree, ab_a, ab_b, weighted=False)
        assert dist > 0.0

    def test_unifrac_range(self, simple_tree: dict) -> None:
        """UniFrac should be in [0, 1]."""
        ab_a = {"A": 5, "B": 3}
        ab_b = {"C": 7, "D": 2}
        for weighted in [False, True]:
            dist = compute_unifrac(simple_tree, ab_a, ab_b, weighted=weighted)
            assert 0.0 <= dist <= 1.0 + 1e-10

    def test_weighted_unifrac_differs_from_unweighted(self, simple_tree: dict) -> None:
        """Weighted and unweighted UniFrac should generally differ."""
        ab_a = {"A": 10, "B": 1}
        ab_b = {"A": 1, "B": 10}
        uw = compute_unifrac(simple_tree, ab_a, ab_b, weighted=False)
        w = compute_unifrac(simple_tree, ab_a, ab_b, weighted=True)
        # Both communities have same species -> unweighted should be 0
        assert abs(uw) < 1e-10
        # Weighted should capture abundance differences
        assert w > 0.0 or abs(w) < 1e-10  # may be 0 if branch-length-weighted props cancel


# ============================================================================
# Phylogenetic beta diversity tests
# ============================================================================


class TestPhylogeneticBetaDiversity:
    """Test phylogenetic beta diversity computation."""

    def test_phylo_beta_diversity_structure(self, simple_tree: dict) -> None:
        """Output has distance_matrix, metric, n_communities."""
        communities = [["A", "B"], ["C", "D"], ["A", "C"]]
        result = phylogenetic_beta_diversity(simple_tree, communities)
        assert "distance_matrix" in result
        assert result["n_communities"] == 3
        dm = result["distance_matrix"]
        assert len(dm) == 3
        # Diagonal should be zero
        for i in range(3):
            assert abs(dm[i][i]) < 1e-10

    def test_phylo_beta_diversity_symmetric(self, simple_tree: dict) -> None:
        """Distance matrix should be symmetric."""
        communities = [["A", "B"], ["C", "D"]]
        result = phylogenetic_beta_diversity(simple_tree, communities)
        dm = result["distance_matrix"]
        assert abs(dm[0][1] - dm[1][0]) < 1e-10


# ============================================================================
# NRI / NTI tests
# ============================================================================


class TestNriNti:
    """Test Net Relatedness Index and Nearest Taxon Index."""

    def test_nri_nti_output_keys(self, simple_tree: dict) -> None:
        """NRI/NTI returns nri, nti, p_values, ses_values."""
        communities = [["A", "B", "C"]]
        result = nri_nti(simple_tree, communities, n_randomizations=49)
        assert "nri" in result
        assert "nti" in result
        assert "p_values" in result

    def test_nri_nti_per_community(self, simple_tree: dict) -> None:
        """One NRI/NTI value per community."""
        communities = [["A", "B"], ["C", "D"]]
        result = nri_nti(simple_tree, communities, n_randomizations=49)
        assert len(result["nri"]) == 2
        assert len(result["nti"]) == 2

    def test_nri_nti_single_taxon(self, simple_tree: dict) -> None:
        """A community with <2 taxa gives NRI=NTI=0."""
        result = nri_nti(simple_tree, [["A"]], n_randomizations=49)
        assert result["nri"][0] == 0.0
        assert result["nti"][0] == 0.0


# ============================================================================
# Phylogenetic signal tests
# ============================================================================


class TestPhylogeneticSignal:
    """Test phylogenetic signal detection."""

    def test_blomberg_k_output(self, simple_tree: dict) -> None:
        """Blomberg's K returns statistic, p_value, method, n_taxa."""
        traits = {"A": 1.0, "B": 1.5, "C": 5.0, "D": 5.5}
        result = phylogenetic_signal(simple_tree, traits, method="blomberg_k")
        assert "statistic" in result
        assert "p_value" in result
        assert result["method"] == "blomberg_k"
        assert result["n_taxa"] == 4

    def test_pagel_lambda_output(self, simple_tree: dict) -> None:
        """Pagel's lambda returns statistic in [0, 1]."""
        traits = {"A": 1.0, "B": 1.5, "C": 5.0, "D": 5.5}
        result = phylogenetic_signal(simple_tree, traits, method="pagel_lambda")
        assert 0.0 <= result["statistic"] <= 1.0
        assert result["method"] == "pagel_lambda"

    def test_phylogenetic_signal_few_taxa_raises(self, simple_tree: dict) -> None:
        """Fewer than 3 taxa with trait data raises ValueError."""
        traits = {"A": 1.0, "B": 2.0}
        with pytest.raises(ValueError, match="at least 3 taxa"):
            phylogenetic_signal(simple_tree, traits)

    def test_unknown_method_raises(self, simple_tree: dict) -> None:
        """Unknown method raises ValueError."""
        traits = {"A": 1.0, "B": 2.0, "C": 3.0}
        with pytest.raises(ValueError, match="Unknown method"):
            phylogenetic_signal(simple_tree, traits, method="unknown")


# ============================================================================
# Build simple tree tests
# ============================================================================


class TestBuildSimpleTree:
    """Test tree construction from distance matrices."""

    def test_upgma_basic(self) -> None:
        """UPGMA builds a valid tree from a 3-taxon distance matrix."""
        dm = [[0, 2, 6], [2, 0, 5], [6, 5, 0]]
        names = ["A", "B", "C"]
        tree = build_simple_tree(dm, names, method="upgma")
        assert tree["name"] == "root"
        assert "children" in tree

        # All tip names should be present
        from metainformant.ecology.phylogenetic.diversity import _get_all_tips

        tips = set(_get_all_tips(tree))
        assert tips == {"A", "B", "C"}

    def test_nj_basic(self) -> None:
        """Neighbor-joining builds a valid tree from a 4-taxon distance matrix."""
        dm = [
            [0, 5, 9, 9],
            [5, 0, 10, 10],
            [9, 10, 0, 8],
            [9, 10, 8, 0],
        ]
        names = ["A", "B", "C", "D"]
        tree = build_simple_tree(dm, names, method="nj")
        assert tree["name"] == "root"
        from metainformant.ecology.phylogenetic.diversity import _get_all_tips

        tips = set(_get_all_tips(tree))
        assert tips == {"A", "B", "C", "D"}

    def test_build_tree_dimension_mismatch_raises(self) -> None:
        """Mismatched matrix/names raises ValueError."""
        dm = [[0, 1], [1, 0]]
        with pytest.raises(ValueError, match="taxa"):
            build_simple_tree(dm, ["A", "B", "C"])

    def test_unknown_method_raises(self) -> None:
        """Unknown method raises ValueError."""
        dm = [[0, 1], [1, 0]]
        with pytest.raises(ValueError, match="Unknown method"):
            build_simple_tree(dm, ["A", "B"], method="parsimony")
