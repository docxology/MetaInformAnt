"""Round-4 regression tests for the ecology package.

Every test pins a result that was wrong before the fix, against an oracle that does not share the code under
test: SciPy (``isotonic_regression``, ``norm``, ``brentq``), closed forms (tetrahedron and octahedron volumes,
the star-tree value of Blomberg's K, hand-summed UniFrac branch lengths) or a brute-force restatement of the
statistic. The defects:

- ``ordination``: the pool-adjacent-violators routine behind NMDS disparities was wrong for runs of three or more
  violators (it gave [2.2, 2.2, 2.2] for [3, 2, 1]); PCoA ``variance_explained`` was normalised by the returned
  axes, so it summed to 1 for any number of axes; ``procrustes`` rotated the second configuration the wrong way
  (m2 = 4 for an exact rotation) and returned a perfect fit for a degenerate one.
- ``macroecology``: ``_normal_cdf`` was off by up to 0.037; ``fit_logseries`` reported alpha = 0 for one species
  and a runaway number for all-singleton communities; ``chi_squared_gof`` did not reject negative frequencies as
  documented; degenerate regressions returned flat fits, which also pulled bootstrap intervals towards 0.
- ``community``: ``calculate_biodiversity_indices`` turned any failure into a list of zeros;
  ``alpha_beta_gamma_diversity`` returned beta = 0 for communities with nothing in common and could report
  gamma < alpha.
- ``indicators``: PERMANOVA / PERMDISP returned F = 0, p = 1 for perfectly separated groups; ANOSIM accepted
  designs without within- or between-group pairs; permutation p-values were rounded to 4 places (0.0).
- ``functional``: FRic in 3+ dimensions was a bounding-box product; FDiv centred on the wrong point.
- ``phylogenetic.diversity``: ``"unweighted_unifrac"`` computed weighted UniFrac (it contains "weighted");
  unweighted UniFrac divided by the whole tree's length; NRI/NTI reported 0 and p = 1 for undefined cases, p = 0
  for finite randomizations, and used the global ``random`` state; "Blomberg's K" and "Pagel's lambda" were not
  those statistics.
- ``compare_sad_models`` ranked models by AICs built from chi-squared statistics on different bins (abundance
  classes, octaves, ranks), which are not comparable; ``rarefaction_curve`` truncated fractional abundances;
  ``cluster_communities`` reported an undefined cophenetic correlation as 0.
- plots: ``output_path`` saved pyplot's current figure instead of the figure of the axes passed in.
"""

from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy import optimize, stats

from metainformant.ecology.analysis import community, functional, indicators, macroecology, ordination
from metainformant.ecology.phylogenetic import diversity

matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Shared trees
# ---------------------------------------------------------------------------


@pytest.fixture
def four_tip_tree() -> Dict[str, Any]:
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


def _balanced_tree(depth: int, prefix: str = "t") -> Dict[str, Any]:
    """Complete binary tree with unit branch lengths and ``2 ** depth`` tips named ``<prefix>00``, ``<prefix>01``..."""
    counter = [0]

    def build(level: int, tag: str) -> Dict[str, Any]:
        if level == 0:
            counter[0] += 1
            return {"name": f"{prefix}{counter[0] - 1:02d}", "branch_length": 1.0, "children": []}
        return {
            "name": f"n{tag}",
            "branch_length": 1.0,
            "children": [build(level - 1, tag + "0"), build(level - 1, tag + "1")],
        }

    tree = build(depth, "r")
    tree["name"] = "root"
    tree["branch_length"] = 0.0
    return tree


# ---------------------------------------------------------------------------
# ordination
# ---------------------------------------------------------------------------


class TestIsotonicRegression:
    def test_three_violators_pool_to_their_mean(self) -> None:
        assert ordination._isotonic_regression(np.array([3.0, 2.0, 1.0])).tolist() == pytest.approx([2.0, 2.0, 2.0])

    @pytest.mark.parametrize("seed", range(8))
    def test_matches_scipy_on_random_weighted_vectors(self, seed: int) -> None:
        rng = np.random.default_rng(seed)
        n = int(rng.integers(3, 40))
        values = rng.normal(size=n)
        weights = rng.uniform(0.2, 3.0, size=n)
        fitted = ordination._isotonic_regression(values, weights)
        expected = optimize.isotonic_regression(values, weights=weights).x
        np.testing.assert_allclose(fitted, expected, atol=1e-10)

    def test_result_is_monotone_and_preserves_the_weighted_mean(self) -> None:
        values = np.array([5.0, 1.0, 4.0, 2.0, 3.0, 9.0, 0.5])
        fitted = ordination._isotonic_regression(values)
        assert np.all(np.diff(fitted) >= -1e-12)
        assert fitted.mean() == pytest.approx(values.mean())


class TestPcoaVarianceExplained:
    @staticmethod
    def _distances() -> List[List[float]]:
        points = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3], [1, 1, 1]], dtype=float)
        return np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(-1)).tolist()

    def test_share_of_total_variance_not_of_the_returned_axes(self) -> None:
        distances = np.array(self._distances())
        centring = np.eye(5) - np.ones((5, 5)) / 5
        eigenvalues = np.sort(np.linalg.eigvalsh(-0.5 * centring @ (distances**2) @ centring))[::-1]
        total = eigenvalues[eigenvalues > 0].sum()

        one_axis = ordination.pcoa(self._distances(), n_components=1)["variance_explained"]

        assert one_axis == pytest.approx([eigenvalues[0] / total], rel=1e-9)
        assert one_axis[0] < 0.7  # the first axis is not "100% of the variance"

    def test_all_positive_axes_add_up_to_one_and_fewer_axes_to_less(self) -> None:
        three_axes = ordination.pcoa(self._distances(), n_components=3)["variance_explained"]
        two_axes = ordination.pcoa(self._distances(), n_components=2)["variance_explained"]
        assert sum(three_axes) == pytest.approx(1.0)
        assert sum(two_axes) < 0.95
        assert two_axes == pytest.approx(three_axes[:2])


class TestProcrustes:
    def test_exact_rotation_is_recovered(self) -> None:
        reference = [[0, 0], [1, 0], [0, 1]]
        rotated = [[0, 0], [0, 1], [-1, 0]]  # the reference turned by 90 degrees

        result = ordination.procrustes(reference, rotated)

        assert result["m2"] == pytest.approx(0.0, abs=1e-12)
        assert result["correlation"] == pytest.approx(1.0)
        np.testing.assert_allclose(np.array(result["transformed_coords"]), np.array(reference, float), atol=1e-9)

    def test_scaled_translated_rotation_in_three_dimensions(self) -> None:
        rng = np.random.default_rng(3)
        reference = rng.normal(size=(7, 3))
        angle = 0.9
        rotation = np.array(
            [[math.cos(angle), -math.sin(angle), 0.0], [math.sin(angle), math.cos(angle), 0.0], [0.0, 0.0, 1.0]]
        )
        target = 2.5 * reference @ rotation + np.array([4.0, -1.0, 2.0])

        result = ordination.procrustes(reference.tolist(), target.tolist())

        assert result["m2"] == pytest.approx(0.0, abs=1e-10)
        np.testing.assert_allclose(np.array(result["transformed_coords"]), reference, atol=1e-8)

    def test_mirror_image_cannot_be_matched_by_a_proper_rotation(self) -> None:
        reference = np.array([[0, 0], [2, 0], [0, 1], [1.5, 2.5]], dtype=float)
        mirrored = reference * np.array([-1.0, 1.0])

        result = ordination.procrustes(reference.tolist(), mirrored.tolist())

        def unit(config: np.ndarray) -> np.ndarray:
            centred = config - config.mean(axis=0)
            return centred / np.linalg.norm(centred)

        singular_values = np.linalg.svd(unit(mirrored).T @ unit(reference), compute_uv=False)
        # excluding the reflection costs the smallest singular value twice: tr = s1 - s2
        expected = singular_values[0] - singular_values[1]
        assert result["correlation"] == pytest.approx(expected, rel=1e-9)
        assert result["m2"] == pytest.approx(2.0 - 2.0 * expected, rel=1e-9)
        assert result["correlation"] < 0.5

    def test_degenerate_configuration_is_rejected_not_reported_as_a_perfect_fit(self) -> None:
        with pytest.raises(ValueError, match="do not all coincide"):
            ordination.procrustes([[0, 0], [1, 0], [0, 1]], [[2, 2], [2, 2], [2, 2]])


# ---------------------------------------------------------------------------
# macroecology
# ---------------------------------------------------------------------------


class TestNormalCdf:
    @pytest.mark.parametrize("z", [-9.0, -6.0, -3.0, -1.5, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0, 8.5])
    def test_matches_scipy(self, z: float) -> None:
        assert macroecology._normal_cdf(z) == pytest.approx(stats.norm.cdf(z), rel=1e-12, abs=1e-300)

    def test_lognormal_expectations_use_the_exact_cdf(self) -> None:
        abundances = [100.0, 60.0, 35.0, 20.0, 12.0, 6.0, 3.0, 2.0, 1.0, 1.0]
        fit = macroecology.fit_lognormal(abundances)
        mu, sigma, s = fit["mu"], fit["sigma"], len(abundances)

        for octave, expected in enumerate(fit["expected_frequencies"]):
            lower, upper = math.log(2**octave), math.log(2 ** (octave + 1))
            reference = s * (stats.norm.cdf((upper - mu) / sigma) - stats.norm.cdf((lower - mu) / sigma))
            assert expected == pytest.approx(reference, rel=1e-9, abs=1e-12)
        veil = stats.norm.sf((math.log(0.5) - mu) / sigma)
        assert fit["s_star"] == pytest.approx(s / veil, rel=1e-9)


class TestChiSquaredGofValidation:
    @pytest.mark.parametrize(
        ("observed", "expected"),
        [([-1.0, 2.0], [1.0, 2.0]), ([1.0, 2.0], [1.0, -2.0]), ([1.0, float("nan")], [1.0, 2.0])],
    )
    def test_negative_or_nan_frequencies_raise(self, observed: List[float], expected: List[float]) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            macroecology.chi_squared_gof(observed, expected)


class TestFisherLogSeries:
    @pytest.mark.parametrize(
        "abundances",
        [[100], [5], [2, 1, 1], [1000, 1], [20, 15, 12, 8, 5, 3, 2, 1, 1, 1], [2] + [1] * 99],
    )
    def test_alpha_solves_the_defining_equation(self, abundances: List[float]) -> None:
        species, individuals = len(abundances), sum(abundances)
        reference = optimize.brentq(lambda a: a * math.log1p(individuals / a) - species, 1e-12, 1e15, xtol=1e-14)

        fit = macroecology.fit_logseries(abundances)

        assert fit["alpha"] == pytest.approx(reference, rel=1e-9)
        assert fit["alpha"] > 0
        assert fit["x"] == pytest.approx(individuals / (individuals + fit["alpha"]))

    def test_one_species_has_positive_alpha(self) -> None:
        # S = 1 used to return alpha = 0 (and x = 0), which does not satisfy S = alpha ln(1 + N / alpha)
        fit = macroecology.fit_logseries([100])
        assert fit["alpha"] == pytest.approx(0.1544496866793, rel=1e-9)
        assert fit["x"] == pytest.approx(100 / (100 + fit["alpha"]))

    def test_all_singletons_have_infinite_alpha(self) -> None:
        # S = N has no finite root; the old Newton iteration stopped at a runaway number (~1.3e8 for four species)
        fit = macroecology.fit_logseries([1, 1, 1, 1])
        assert math.isinf(fit["alpha"]) and fit["alpha"] > 0
        assert fit["x"] == 0.0
        assert fit["expected_frequencies"] == [4.0]

    def test_fractional_abundances_below_the_species_count_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="counts"):
            macroecology.fit_logseries([0.4, 0.3, 0.2])


class TestDegenerateRegressions:
    @pytest.mark.parametrize(
        "call",
        [
            lambda: macroecology.species_area_power([10, 10, 10], [5, 6, 7]),
            lambda: macroecology.species_area_logarithmic([10, 10, 10], [5, 6, 7]),
            lambda: macroecology.distance_decay([3.0, 3.0, 3.0], [0.9, 0.5, 0.2]),
            lambda: macroecology.metabolic_scaling([2.0, 2.0, 2.0], [1.0, 2.0, 3.0]),
            lambda: macroecology.taylors_power_law([4.0, 4.0, 4.0], [1.0, 2.0, 3.0]),
        ],
    )
    def test_constant_predictor_is_an_error_not_a_flat_fit(self, call: Any) -> None:
        with pytest.raises(ValueError, match="must not all be equal"):
            call()

    def test_bootstrap_interval_does_not_count_slopes_of_unfit_resamples_as_zero(self) -> None:
        # With three areas, 1 resample in 9 draws a single area; its "slope" used to be counted as 0, which at
        # 1000 resamples pulled the lower end of the interval to exactly 0 for a clearly positive exponent.
        fit = macroecology.species_area_power([1, 10, 100], [5, 14, 50], n_bootstrap=1000)
        lower, upper = fit["ci_z"]
        assert fit["z"] > 0.4
        assert lower > 0.3
        assert lower <= fit["z"] <= upper

    def test_distance_decay_reports_how_many_pairs_were_fitted(self) -> None:
        result = macroecology.distance_decay([1.0, 2.0, 3.0, 4.0, 5.0], [0.9, 0.0, 0.5, 0.3, 0.1])
        assert result["n_used"] == 4  # the zero similarity has no logarithm and is left out


# ---------------------------------------------------------------------------
# community
# ---------------------------------------------------------------------------


class TestBiodiversityIndices:
    def test_a_failing_index_is_an_error_not_a_list_of_zeros(self) -> None:
        with pytest.raises(TypeError):
            community.calculate_biodiversity_indices([[1.0, "many"], [2.0, 3.0]])

    def test_mapping_input_is_rejected_instead_of_zeroed(self) -> None:
        with pytest.raises(TypeError, match="list of abundance vectors"):
            community.calculate_biodiversity_indices({"a": [1.0, 2.0]})  # type: ignore[arg-type]


class TestAlphaBetaGamma:
    def test_disjoint_single_species_communities_have_beta_equal_to_their_number(self) -> None:
        result = community.alpha_beta_gamma_diversity([[5, 0, 0], [0, 5, 0], [0, 0, 5]])
        assert result["alpha"] == pytest.approx(0.0)
        assert result["gamma"] == pytest.approx(math.log(3))
        assert result["beta"] == pytest.approx(3.0)  # was 0.0: alpha was 0, so the ratio fell back to 0

    def test_identical_communities_have_beta_one(self) -> None:
        result = community.alpha_beta_gamma_diversity([[5, 3, 1], [5, 3, 1]])
        assert result["beta"] == pytest.approx(1.0)
        assert result["gamma"] == pytest.approx(result["alpha"])

    def test_gamma_never_falls_below_alpha_for_unequal_sample_sizes(self) -> None:
        # a 1001-individual community and a 20-individual one: an unweighted mean of the entropies exceeded the
        # entropy of the pool
        result = community.alpha_beta_gamma_diversity([[1000, 1], [5, 5, 5, 5]])
        assert result["gamma"] >= result["alpha"]
        assert result["beta"] >= 1.0
        totals = np.array([1001.0, 20.0])
        within = [stats.entropy([1000, 1]), stats.entropy([5, 5, 5, 5])]
        pooled = stats.entropy([1005, 6, 5, 5])
        assert result["alpha"] == pytest.approx(float(totals @ within / totals.sum()))
        assert result["gamma"] == pytest.approx(pooled)

    def test_no_individuals_gives_zeros(self) -> None:
        assert community.alpha_beta_gamma_diversity([[0, 0], [0, 0]]) == {"alpha": 0.0, "beta": 0.0, "gamma": 0.0}


# ---------------------------------------------------------------------------
# indicators
# ---------------------------------------------------------------------------


def _separated_groups_matrix() -> List[List[float]]:
    """Six observations in two tight groups of three (within distance 0) that are 5 apart."""
    return [[0.0 if (i < 3) == (j < 3) else 5.0 for j in range(6)] for i in range(6)]


class TestPermanovaDegenerateDesigns:
    def test_perfect_separation_has_infinite_f_and_unit_r_squared(self) -> None:
        labels = ["a", "a", "a", "b", "b", "b"]
        result = indicators.permanova(_separated_groups_matrix(), labels, n_permutations=499, seed=1)
        assert math.isinf(result["f_statistic"]) and result["f_statistic"] > 0
        assert result["r_squared"] == pytest.approx(1.0)
        # only the 2 labellings that reproduce the grouping (of 20) match it: p is near 0.1, not 1
        assert result["p_value"] < 0.3

    def test_more_observations_than_groups_are_required(self) -> None:
        with pytest.raises(ValueError, match="more observations than groups"):
            indicators.permanova([[0, 1], [1, 0]], ["a", "b"])

    def test_all_zero_distances_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="All distances are zero"):
            indicators.permanova([[0, 0, 0], [0, 0, 0], [0, 0, 0]], ["a", "a", "b"])


class TestPermdispDegenerateDesigns:
    @staticmethod
    def _unequal_dispersion_matrix() -> List[List[float]]:
        """Two groups of three points, each an equilateral triangle: side 1 in one group, side 3 in the other."""
        size = 6
        matrix = [[0.0] * size for _ in range(size)]
        for i in range(size):
            for j in range(size):
                if i == j:
                    continue
                if i < 3 and j < 3:
                    matrix[i][j] = 1.0
                elif i >= 3 and j >= 3:
                    matrix[i][j] = 3.0
                else:
                    matrix[i][j] = 5.0
        return matrix

    def test_groups_of_unequal_dispersion_but_no_spread_inside_have_infinite_f(self) -> None:
        labels = ["a", "a", "a", "b", "b", "b"]
        result = indicators.multivariate_dispersion(
            self._unequal_dispersion_matrix(), labels, n_permutations=199, seed=2
        )
        assert math.isinf(result["f_statistic"]) and result["f_statistic"] > 0
        assert result["group_dispersions"]["a"] < result["group_dispersions"]["b"]
        assert result["p_value"] < 0.5

    def test_more_observations_than_groups_are_required(self) -> None:
        with pytest.raises(ValueError, match="more observations than groups"):
            indicators.multivariate_dispersion([[0, 1], [1, 0]], ["a", "b"])


class TestAnosimDegenerateDesigns:
    def test_a_single_group_has_no_between_group_pairs(self) -> None:
        with pytest.raises(ValueError, match="two distinct groups"):
            indicators.anosim([[0, 1, 2], [1, 0, 3], [2, 3, 0]], ["a", "a", "a"])

    def test_all_groups_of_one_have_no_within_group_pairs(self) -> None:
        with pytest.raises(ValueError, match="two observations"):
            indicators.anosim([[0, 1, 2], [1, 0, 3], [2, 3, 0]], ["a", "b", "c"])


class TestDistanceMatrixValidation:
    def test_asymmetric_matrix_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="not symmetric"):
            indicators.anosim([[0, 1, 2], [9, 0, 3], [2, 3, 0]], ["a", "a", "b"])

    def test_non_finite_distance_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="non-finite"):
            indicators.permanova([[0, float("nan"), 2], [float("nan"), 0, 3], [2, 3, 0]], ["a", "a", "b"])


class TestPermutationPValuesAreNotRoundedToZero:
    def test_indval_p_value_keeps_its_resolution(self) -> None:
        sites = [[10.0]] * 10 + [[0.0]] * 10
        labels = ["a"] * 10 + ["b"] * 10

        result = indicators.indval(sites, labels, n_permutations=30000, seed=0)[0]

        # (hits + 1) / 30001 is at least 3.3e-5: a p-value rounded to 4 places printed it as 0.0
        assert result["p_value"] >= 1.0 / 30001 - 1e-12
        assert 0.0 < result["p_value"] < 1e-3
        assert result["p_value"] * 30001 == pytest.approx(round(result["p_value"] * 30001), abs=1e-6)


# ---------------------------------------------------------------------------
# functional
# ---------------------------------------------------------------------------


class TestFunctionalRichnessHull:
    def test_unit_tetrahedron_has_volume_one_sixth(self) -> None:
        # a bounding box would give 1.0
        value = functional.functional_richness([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        assert value == pytest.approx(1.0 / 6.0)

    def test_octahedron_has_volume_four_thirds(self) -> None:
        points = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        assert functional.functional_richness(points) == pytest.approx(4.0 / 3.0)  # box: 8

    def test_interior_points_do_not_change_the_volume(self) -> None:
        cube = [[x, y, z] for x in (0, 2) for y in (0, 2) for z in (0, 2)]
        assert functional.functional_richness(cube + [[1, 1, 1], [0.5, 1.5, 1.0]]) == pytest.approx(8.0)

    def test_coplanar_species_enclose_no_volume(self) -> None:
        # four species on one plane: the old product of the two non-constant ranges (an area) was returned
        assert functional.functional_richness([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]) == 0.0

    def test_too_few_species_for_the_dimension_enclose_no_volume(self) -> None:
        assert functional.functional_richness([[0, 0, 0], [1, 2, 3], [4, 1, 0]]) == 0.0

    def test_one_dimension_returns_a_float(self) -> None:
        value = functional.functional_richness([[0], [5], [10]])
        assert isinstance(value, float) and value == 10.0


class TestFunctionalDivergenceCentre:
    @staticmethod
    def _fdiv(abundances: List[float], traits: List[float]) -> float:
        """Villeger et al. 2008 for one trait: G is the midpoint of the extreme species, not the weighted mean."""
        centre = (min(traits) + max(traits)) / 2.0
        distances = [abs(t - centre) for t in traits]
        mean_distance = sum(distances) / len(distances)
        total = sum(abundances)
        weights = [a / total for a in abundances]
        delta = sum(w * (d - mean_distance) for w, d in zip(weights, distances))
        delta_abs = sum(w * abs(d - mean_distance) for w, d in zip(weights, distances))
        return (delta + mean_distance) / (delta_abs + mean_distance)

    def test_asymmetric_abundances_use_the_hull_vertex_centre(self) -> None:
        traits = [0.0, 5.0, 10.0]
        abundances = [1.0, 1.0, 10.0]
        value = functional.functional_divergence([[t] for t in traits], abundances)
        assert value == pytest.approx(self._fdiv(abundances, traits))
        assert value == pytest.approx(0.8918918918918919)

    def test_interior_species_do_not_move_the_centre(self) -> None:
        # the centre depends on the hull only, so abundance on interior species changes the weights, not G
        traits = [0.0, 2.0, 3.0, 7.0, 10.0]
        abundances = [3.0, 1.0, 6.0, 2.0, 1.0]
        value = functional.functional_divergence([[t] for t in traits], abundances)
        assert value == pytest.approx(self._fdiv(abundances, traits))


# ---------------------------------------------------------------------------
# phylogenetic diversity
# ---------------------------------------------------------------------------


class TestUniFracMetricNames:
    def test_unweighted_unifrac_is_unweighted(self, four_tip_tree: Dict[str, Any]) -> None:
        communities = [["A", "B"], ["A", "C", "D"]]
        presence = [{t: 1.0 for t in c} for c in communities]
        unweighted = diversity.compute_unifrac(four_tip_tree, presence[0], presence[1], weighted=False)
        weighted = diversity.compute_unifrac(four_tip_tree, presence[0], presence[1], weighted=True)
        assert unweighted != pytest.approx(weighted, abs=1e-3)  # the metrics must be distinguishable here

        result = diversity.phylogenetic_beta_diversity(four_tip_tree, communities, metric="unweighted_unifrac")

        assert result["distance_matrix"][0][1] == pytest.approx(unweighted)
        assert result["distance_matrix"][1][0] == pytest.approx(unweighted)

    def test_weighted_and_default_metrics(self, four_tip_tree: Dict[str, Any]) -> None:
        communities = [["A", "B"], ["A", "C", "D"]]
        presence = [{t: 1.0 for t in c} for c in communities]
        weighted = diversity.phylogenetic_beta_diversity(four_tip_tree, communities, metric="weighted_unifrac")
        default = diversity.phylogenetic_beta_diversity(four_tip_tree, communities)
        assert weighted["distance_matrix"][0][1] == pytest.approx(
            diversity.compute_unifrac(four_tip_tree, presence[0], presence[1], weighted=True)
        )
        assert default["distance_matrix"][0][1] == pytest.approx(
            diversity.compute_unifrac(four_tip_tree, presence[0], presence[1], weighted=False)
        )

    def test_unknown_metric_raises_instead_of_computing_something_else(self, four_tip_tree: Dict[str, Any]) -> None:
        with pytest.raises(ValueError, match="Unsupported metric"):
            diversity.phylogenetic_beta_diversity(four_tip_tree, [["A"], ["B"]], metric="bray_curtis")


class TestUnweightedUniFracDenominator:
    def test_branches_neither_community_reaches_are_outside_the_denominator(
        self, four_tip_tree: Dict[str, Any]
    ) -> None:
        # A = {A}, B = {B}: branch A (0.5) and B (0.4) are unique, internal_1 (0.3) is shared, and the C/D side
        # (0.2 + 0.6 + 0.3) leads to neither. unique / observed = 0.9 / 1.2; dividing by the whole tree
        # (2.3) gave 0.391.
        value = diversity.compute_unifrac(four_tip_tree, {"A": 1}, {"B": 1}, weighted=False)
        assert value == pytest.approx(0.9 / 1.2)

    def test_partial_overlap_by_hand(self, four_tip_tree: Dict[str, Any]) -> None:
        # {A, B} against {A, C, D}: unique branches B 0.4, C 0.6, D 0.3, internal_2 0.2 = 1.5 of 2.3 observed
        value = diversity.compute_unifrac(four_tip_tree, {"A": 1, "B": 1}, {"A": 1, "C": 1, "D": 1}, weighted=False)
        assert value == pytest.approx(1.5 / 2.3)

    def test_adding_unrelated_tips_to_the_tree_does_not_change_the_distance(
        self, four_tip_tree: Dict[str, Any]
    ) -> None:
        bigger = {
            "name": "outer",
            "branch_length": 0.0,
            "children": [
                four_tip_tree,
                {
                    "name": "extra",
                    "branch_length": 1.0,
                    "children": [{"name": f"X{i}", "branch_length": 1.0, "children": []} for i in range(5)],
                },
            ],
        }
        small = diversity.compute_unifrac(four_tip_tree, {"A": 1}, {"B": 1}, weighted=False)
        large = diversity.compute_unifrac(bigger, {"A": 1}, {"B": 1}, weighted=False)
        assert large == pytest.approx(small)


class TestNriNtiHonestStatistics:
    def test_p_value_of_a_finite_randomization_is_never_zero(self) -> None:
        tree = _balanced_tree(8)  # 256 tips
        tips = diversity._get_all_tips(tree)
        result = diversity.nri_nti(tree, [[tips[0], tips[1]]], n_randomizations=19, seed=0)
        # the sister pair is at least as clustered as every one of 19 random pairs: (1 + 0) / (1 + 19)
        assert result["p_values"]["nri"][0] == pytest.approx(1.0 / 20.0)
        assert result["p_values"]["nti"][0] == pytest.approx(1.0 / 20.0)
        assert result["nri"][0] > 0

    def test_undefined_statistics_are_nan(self, four_tip_tree: Dict[str, Any]) -> None:
        result = diversity.nri_nti(four_tip_tree, [["A"], ["Z", "Q"], ["A", "B", "C"]], n_randomizations=29, seed=1)
        for key in ("nri", "nti"):
            assert math.isnan(result[key][0]) and math.isnan(result[key][1])
            assert not math.isnan(result[key][2])
        assert math.isnan(result["p_values"]["nri"][0]) and math.isnan(result["p_values"]["nti"][1])
        assert 0.0 < result["p_values"]["nri"][2] <= 1.0

    def test_zero_spread_null_makes_the_standardised_index_undefined(self) -> None:
        star = {
            "name": "root",
            "branch_length": 0.0,
            "children": [{"name": n, "branch_length": 1.0, "children": []} for n in "ABCDE"],
        }
        result = diversity.nri_nti(star, [["A", "B", "C"]], n_randomizations=29, seed=3)
        assert math.isnan(result["nri"][0]) and math.isnan(result["nti"][0])
        assert result["p_values"]["nri"][0] == pytest.approx(1.0)  # every null value ties the observed one

    def test_global_random_state_is_neither_used_nor_changed(self, four_tip_tree: Dict[str, Any]) -> None:
        random.seed(12345)
        state = random.getstate()
        diversity.nri_nti(four_tip_tree, [["A", "B", "C"]], n_randomizations=25)  # seed=None
        diversity.phylogenetic_signal(four_tip_tree, {"A": 1.0, "B": 1.5, "C": 5.0, "D": 5.5})
        assert random.getstate() == state


class TestBlombergK:
    @staticmethod
    def _reference_k(traits: np.ndarray, vcv: np.ndarray) -> float:
        """Blomberg et al. 2003 written out with explicit matrices."""
        n = len(traits)
        inverse = np.linalg.inv(vcv)
        ones = np.ones(n)
        a_hat = (ones @ inverse @ traits) / (ones @ inverse @ ones)
        mse = ((traits - a_hat) @ inverse @ (traits - a_hat)) / (n - 1)
        mse0 = ((traits - traits.mean()) @ (traits - traits.mean())) / (n - 1)
        expected = (np.trace(vcv) - n / (ones @ inverse @ ones)) / (n - 1)
        return float((mse0 / mse) / expected)

    def test_star_tree_gives_k_of_exactly_one_for_any_trait(self) -> None:
        # With independent tips of equal variance the phylogenetic mean is the arithmetic mean and the expected
        # ratio is the root-to-tip length: K = 1 whatever the data
        star = {
            "name": "root",
            "branch_length": 0.0,
            "children": [{"name": n, "branch_length": 2.0, "children": []} for n in "ABCDEF"],
        }
        traits = {"A": 1.0, "B": 7.0, "C": -2.0, "D": 3.5, "E": 3.0, "F": 11.0}
        result = diversity.phylogenetic_signal(star, traits, "blomberg_k", seed=0)
        assert result["statistic"] == pytest.approx(1.0, rel=1e-12)

    def test_matches_the_matrix_formula(self, four_tip_tree: Dict[str, Any]) -> None:
        traits = {"A": 1.0, "B": 1.5, "C": 5.0, "D": 5.5}
        vcv = np.array([[0.8, 0.3, 0.0, 0.0], [0.3, 0.7, 0.0, 0.0], [0.0, 0.0, 0.8, 0.2], [0.0, 0.0, 0.2, 0.5]])
        result = diversity.phylogenetic_signal(four_tip_tree, traits, "blomberg_k", seed=1)
        assert result["statistic"] == pytest.approx(self._reference_k(np.array([1.0, 1.5, 5.0, 5.5]), vcv), rel=1e-9)

    def test_detects_brownian_signal_and_not_noise(self) -> None:
        tree = _balanced_tree(5)
        tips = diversity._get_all_tips(tree)
        covariance = diversity._phylogenetic_vcv(tree, tips)
        rng = np.random.default_rng(0)
        brownian = np.linalg.cholesky(covariance) @ rng.normal(size=len(tips))
        noise = rng.normal(size=len(tips))

        with_signal = diversity.phylogenetic_signal(tree, dict(zip(tips, brownian)), "blomberg_k", seed=1)
        without = diversity.phylogenetic_signal(tree, dict(zip(tips, noise)), "blomberg_k", seed=1)

        assert with_signal["statistic"] == pytest.approx(self._reference_k(brownian, covariance), rel=1e-9)
        assert with_signal["p_value"] < 0.01
        assert without["p_value"] > 0.05
        assert without["statistic"] < with_signal["statistic"]

    def test_a_constant_trait_is_an_error_not_k_equals_zero_with_p_equals_one(
        self, four_tip_tree: Dict[str, Any]
    ) -> None:
        with pytest.raises(ValueError, match="no variation"):
            diversity.phylogenetic_signal(four_tip_tree, {"A": 2.0, "B": 2.0, "C": 2.0, "D": 2.0}, "blomberg_k")

    def test_seeded_p_value_is_reproducible(self, four_tip_tree: Dict[str, Any]) -> None:
        traits = {"A": 1.0, "B": 1.5, "C": 5.0, "D": 5.5}
        first = diversity.phylogenetic_signal(four_tip_tree, traits, "blomberg_k", seed=9)
        second = diversity.phylogenetic_signal(four_tip_tree, traits, "blomberg_k", seed=9)
        assert first == second


class TestPagelLambda:
    @staticmethod
    def _log_likelihood(traits: np.ndarray, vcv: np.ndarray, lam: float) -> float:
        n = len(traits)
        scaled = np.diag(np.diag(vcv)) + lam * (vcv - np.diag(np.diag(vcv)))
        inverse = np.linalg.inv(scaled)
        ones = np.ones(n)
        mean = (ones @ inverse @ traits) / (ones @ inverse @ ones)
        residual = traits - mean
        sigma2 = (residual @ inverse @ residual) / n
        return float(-0.5 * (n * math.log(2 * math.pi * sigma2) + np.log(np.linalg.eigvalsh(scaled)).sum() + n))

    def test_estimate_is_the_maximum_of_the_brownian_likelihood(self) -> None:
        tree = _balanced_tree(4)
        tips = diversity._get_all_tips(tree)
        covariance = diversity._phylogenetic_vcv(tree, tips)
        rng = np.random.default_rng(5)
        # a trait between pure noise and Brownian motion, so that the optimum is interior
        trait = 0.5 * (np.linalg.cholesky(covariance) @ rng.normal(size=len(tips))) + 1.5 * rng.normal(size=len(tips))

        result = diversity.phylogenetic_signal(tree, dict(zip(tips, trait)), "pagel_lambda")

        grid = np.linspace(0.0, 1.0, 1001)
        likelihoods = [self._log_likelihood(trait, covariance, float(lam)) for lam in grid]
        assert result["statistic"] == pytest.approx(float(grid[int(np.argmax(likelihoods))]), abs=2e-3)
        assert result["log_likelihood"] == pytest.approx(max(likelihoods), abs=1e-6)
        assert result["log_likelihood_lambda0"] == pytest.approx(likelihoods[0], abs=1e-9)
        statistic = 2.0 * (max(likelihoods) - likelihoods[0])
        assert result["p_value"] == pytest.approx(float(stats.chi2.sf(statistic, df=1)), rel=1e-3)

    def test_detects_brownian_signal_and_not_noise(self) -> None:
        tree = _balanced_tree(5)
        tips = diversity._get_all_tips(tree)
        covariance = diversity._phylogenetic_vcv(tree, tips)
        rng = np.random.default_rng(0)
        brownian = np.linalg.cholesky(covariance) @ rng.normal(size=len(tips))
        noise = rng.normal(size=len(tips))

        with_signal = diversity.phylogenetic_signal(tree, dict(zip(tips, brownian)), "pagel_lambda")
        without = diversity.phylogenetic_signal(tree, dict(zip(tips, noise)), "pagel_lambda")

        assert with_signal["statistic"] == pytest.approx(1.0, abs=1e-6)
        assert with_signal["p_value"] < 1e-6
        assert without["statistic"] == pytest.approx(0.0, abs=1e-6)
        assert without["p_value"] == pytest.approx(1.0)

    def test_is_deterministic(self, four_tip_tree: Dict[str, Any]) -> None:
        traits = {"A": 1.0, "B": 1.5, "C": 5.0, "D": 5.5}
        assert diversity.phylogenetic_signal(four_tip_tree, traits, "pagel_lambda") == diversity.phylogenetic_signal(
            four_tip_tree, traits, "pagel_lambda"
        )

    def test_star_tree_has_no_identifiable_lambda(self) -> None:
        star = {
            "name": "root",
            "branch_length": 0.0,
            "children": [{"name": n, "branch_length": 1.0, "children": []} for n in "ABCD"],
        }
        with pytest.raises(ValueError, match="not identifiable"):
            diversity.phylogenetic_signal(star, {"A": 1.0, "B": 2.0, "C": 4.0, "D": 8.0}, "pagel_lambda")


# ---------------------------------------------------------------------------
# SAD model comparison, rarefaction, clustering
# ---------------------------------------------------------------------------


class TestCompareSadModels:
    COMMUNITY = [100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 2.0, 1.0, 1.0]

    @staticmethod
    def _aicc(observed: List[float], predicted: List[float], k: int) -> float:
        n = len(observed)
        rss = sum((math.log(o) - math.log(e)) ** 2 for o, e in zip(observed, predicted))
        value = n * math.log(rss / n) + 2 * k
        if n > k + 1:
            value += 2 * k * (k + 1) / (n - k - 1)
        return value

    def test_every_model_is_scored_the_same_way_on_its_rank_abundance_curve(self) -> None:
        result = macroecology.compare_sad_models(self.COMMUNITY)
        observed = sorted(self.COMMUNITY, reverse=True)
        broken_stick = macroecology.fit_broken_stick(self.COMMUNITY)["expected_abundances"]
        geometric = macroecology.fit_geometric_series(self.COMMUNITY)
        assert result["broken_stick"]["aic"] == pytest.approx(self._aicc(observed, broken_stick, 0), rel=1e-9)
        assert result["geometric_series"]["aic"] == pytest.approx(
            self._aicc(observed, geometric["expected_abundances"], 1)
        )
        for name, k in (("logseries", 1), ("lognormal", 2)):
            assert result[name]["aic"] == pytest.approx(self._aicc(observed, result[name]["expected_abundances"], k))

    def test_lognormal_prediction_is_the_quantile_curve(self) -> None:
        result = macroecology.compare_sad_models(self.COMMUNITY)
        fit = macroecology.fit_lognormal(self.COMMUNITY)
        s = len(self.COMMUNITY)
        expected = [
            math.exp(fit["mu"] + fit["sigma"] * stats.norm.ppf(1.0 - (rank - 0.5) / s)) for rank in range(1, s + 1)
        ]
        np.testing.assert_allclose(result["lognormal"]["expected_abundances"], expected, rtol=1e-9)

    def test_logseries_prediction_follows_scipys_log_series(self) -> None:
        result = macroecology.compare_sad_models(self.COMMUNITY)
        fit = macroecology.fit_logseries(self.COMMUNITY)
        s = len(self.COMMUNITY)
        for rank, predicted in enumerate(result["logseries"]["expected_abundances"], start=1):
            # largest a with (expected species with >= a individuals) >= rank - 1/2
            candidates = [a for a in range(1, 4000) if s * stats.logser.sf(a - 1, fit["x"]) >= rank - 0.5]
            assert predicted == max(candidates)

    def test_data_that_follow_a_model_exactly_rank_it_first(self) -> None:
        n, s = 100.0, 6
        broken_stick = [(n / s) * sum(1.0 / i for i in range(k, s + 1)) for k in range(1, s + 1)]
        assert macroecology.compare_sad_models(broken_stick)["broken_stick"]["rank"] == 1

        k = 0.4
        raw = [k * (1 - k) ** r for r in range(8)]
        geometric = [100.0 * value / sum(raw) for value in raw]
        assert macroecology.compare_sad_models(geometric)["geometric_series"]["rank"] == 1

    def test_the_chi_squared_of_different_bins_no_longer_decides_the_ranking(self) -> None:
        # The log-series' chi-squared is taken over ~100 abundance classes (242.7 here against 4.8 for the
        # lognormal's octaves). Turned into AICs those were 245 and 10.5: a 235-point gap that only measured the
        # number of bins, and put the log-series last. On a common footing the four models are within a few points.
        result = macroecology.compare_sad_models(self.COMMUNITY)
        assert result["logseries"]["gof"]["chi_squared"] > 100 > result["lognormal"]["gof"]["chi_squared"]
        spread = max(entry["aic"] for entry in result.values()) - min(entry["aic"] for entry in result.values())
        assert spread < 40
        assert result["logseries"]["rank"] < 4  # no longer last
        assert sorted(entry["rank"] for entry in result.values()) == [1, 2, 3, 4]

    def test_all_singleton_community_does_not_break_the_comparison(self) -> None:
        result = macroecology.compare_sad_models([1.0, 1.0, 1.0, 1.0, 1.0])
        assert result["logseries"]["expected_abundances"] == [1.0] * 5
        assert math.isinf(result["logseries"]["parameters"]["alpha"])


class TestRarefactionNeedsCounts:
    def test_fractional_abundances_are_rejected_not_truncated(self) -> None:
        with pytest.raises(ValueError, match="whole-number"):
            community.rarefaction_curve([10.0, 2.7, 1.0])

    def test_whole_number_floats_are_fine(self) -> None:
        curve = community.rarefaction_curve([10.0, 5.0, 1.0], max_samples=3)
        assert [n for n, _ in curve] == [1, 2, 3]
        assert curve[0][1] == pytest.approx(1.0)


class TestCopheneticCorrelationUndefined:
    def test_two_observations_have_no_correlation(self) -> None:
        result = indicators.cluster_communities([[0.0, 1.0], [1.0, 0.0]])
        assert math.isnan(result["cophenetic_correlation"])  # was 0.0

    def test_equal_distances_have_no_correlation(self) -> None:
        result = indicators.cluster_communities([[0.0, 2.0, 2.0], [2.0, 0.0, 2.0], [2.0, 2.0, 0.0]])
        assert math.isnan(result["cophenetic_correlation"])

    def test_defined_correlations_are_unchanged(self) -> None:
        distances = [[0, 2, 6], [2, 0, 5], [6, 5, 0]]
        result = indicators.cluster_communities(distances, method="upgma")
        # UPGMA merges (0, 1) at 2 then joins 2 at the mean 5.5: cophenetic distances 2, 5.5, 5.5
        reference = stats.pearsonr([2, 6, 5], [2, 5.5, 5.5]).statistic
        assert result["cophenetic_correlation"] == pytest.approx(reference, abs=1e-5)


# ---------------------------------------------------------------------------
# plot export
# ---------------------------------------------------------------------------


class TestPlotExportSavesTheAxesFigure:
    @staticmethod
    def _inked_fraction(path: Path) -> float:
        image = plt.imread(path)
        return float((image[..., :3] < 0.99).any(axis=-1).mean())

    def test_ecology_plot_saves_the_figure_of_the_axes_passed_in(self, tmp_path: Path) -> None:
        from metainformant.ecology.visualization.visualization import plot_species_abundance_distribution

        figure, axes = plt.subplots()
        other = plt.figure()  # a later, empty figure is now pyplot's current one
        target = tmp_path / "abundance.png"
        try:
            plot_species_abundance_distribution(np.array([10.0, 5.0, 3.0, 1.0]), ax=axes, output_path=target)
        finally:
            plt.close(figure)
            plt.close(other)
        assert self._inked_fraction(target) > 0.01

    def test_phenotype_plot_saves_the_figure_of_the_axes_passed_in(self, tmp_path: Path) -> None:
        pd = pytest.importorskip("pandas")
        from metainformant.phenotype.visualization.visualization import plot_trait_correlation_matrix

        rng = np.random.default_rng(0)
        traits = pd.DataFrame(rng.normal(size=(30, 3)), columns=["length", "width", "mass"])
        figure, axes = plt.subplots()
        other = plt.figure()  # a later, empty figure is now pyplot's current one
        target = tmp_path / "correlations.png"
        try:
            plot_trait_correlation_matrix(traits, ax=axes, output_path=target)
        finally:
            plt.close(figure)
            plt.close(other)
        assert self._inked_fraction(target) > 0.01
