"""Metabolomics analysis functions do not turn degenerate input into numbers.

Regression tests for ``metainformant.metabolomics.analysis.identification``:

* ``differential_abundance`` accepted groups of one sample. Welch's test needs each group's variance, which one
  observation cannot give, so the variance was taken to be 0.0 and two single observations that differed got
  ``t = inf`` and ``p = 0.0`` (a one-versus-one "highly significant" difference);
* ``identify_metabolites`` divided by the observed m/z, so a negative m/z had a negative error that is "below" every
  tolerance: it matched the whole database, with scores above 1. Zero divided by zero and NaN matched nothing
  without saying why;
* ``missing_value_imputation(method="knn")`` took the "most correlated" rows even when they had no correlation at
  all (too little overlap, correlation -inf), and averaged their unrelated values into the gap.

Oracles: scipy's Welch t-test and a direct implementation of the k-nearest-row rule. Real computation; no test
doubles (real-implementation policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy import stats

from metainformant.metabolomics.analysis.identification import (
    differential_abundance,
    identify_metabolites,
    identify_with_adducts,
    missing_value_imputation,
)
from metainformant.metabolomics.pathways.enrichment import pathway_activity_scoring

DATABASE = {"glucose": 180.0634, "fructose": 180.0634, "lactate": 90.0317}


class TestDifferentialAbundanceNeedsTwoSamplesPerGroup:
    DATA = np.array([[1.0, 2.0, 3.0, 4.0, 5.0], [5.0, 5.0, 5.0, 9.0, 9.0]])

    @pytest.mark.parametrize(
        ("group_a", "group_b"),
        [([0], [1]), ([0], [1, 2, 3]), ([0, 1, 2], [3]), ([], [1, 2]), ([0, 1], [])],
    )
    def test_a_group_with_fewer_than_two_samples_raises(self, group_a: list[int], group_b: list[int]) -> None:
        with pytest.raises(ValueError, match=r"at least 2 samples in each group"):
            differential_abundance(self.DATA, group_a, group_b)

    def test_one_versus_one_no_longer_reports_p_equal_zero(self) -> None:
        # two single observations that differ: used to give t = -inf and p = 0.0
        with pytest.raises(ValueError, match=r"got 1 in group A and 1 in group B"):
            differential_abundance(np.array([[1.0, 2.0]]), [0], [1])

    def test_two_per_group_is_the_smallest_welch_test_and_matches_scipy(self) -> None:
        rng = np.random.default_rng(1)
        data = rng.normal(size=(6, 4))

        t_stats, p_values = differential_abundance(data, [0, 1], [2, 3])

        reference = stats.ttest_ind(data[:, :2], data[:, 2:], axis=1, equal_var=False)
        np.testing.assert_allclose(t_stats, reference.statistic, rtol=1e-10)
        np.testing.assert_allclose(p_values, reference.pvalue, rtol=1e-9)

    def test_unequal_groups_still_match_scipy(self) -> None:
        rng = np.random.default_rng(2)
        data = rng.normal(size=(5, 9))

        t_stats, p_values = differential_abundance(data, [0, 1], [2, 3, 4, 5, 6, 7, 8])

        reference = stats.ttest_ind(data[:, :2], data[:, 2:], axis=1, equal_var=False)
        np.testing.assert_allclose(t_stats, reference.statistic, rtol=1e-10)
        np.testing.assert_allclose(p_values, reference.pvalue, rtol=1e-9)

    def test_constant_groups_keep_their_documented_degenerate_values(self) -> None:
        t_stats, p_values = differential_abundance(self.DATA[:, :4], [0, 1], [2, 3])

        # metabolite 1 is constant (5, 5) in group A and (5, 9) -> not constant in B; metabolite 0: 1,2 vs 3,4
        assert np.isfinite(t_stats).all() and ((0.0 <= p_values) & (p_values <= 1.0)).all()
        both_constant, p_constant = differential_abundance(np.array([[2.0, 2.0, 7.0, 7.0]]), [0, 1], [2, 3])
        assert both_constant[0] == -np.inf and p_constant[0] == 0.0  # constants that differ
        same_constant, p_same = differential_abundance(np.array([[2.0, 2.0, 2.0, 2.0]]), [0, 1], [2, 3])
        assert same_constant[0] == 0.0 and p_same[0] == 1.0  # constants that agree

    def test_a_missing_value_gives_nan_for_that_metabolite_only(self) -> None:
        data = np.array([[1.0, 2.0, 5.0, 6.0, 7.0], [1.0, np.nan, 5.0, 6.0, 7.0]])

        t_stats, p_values = differential_abundance(data, [0, 1], [2, 3, 4])

        assert np.isfinite(t_stats[0]) and np.isnan(t_stats[1])
        assert np.isfinite(p_values[0]) and np.isnan(p_values[1])


class TestObservedMzMustBePositiveAndFinite:
    @pytest.mark.parametrize("bad", [-180.0634, -1e-9, 0.0, float("nan"), float("inf"), float("-inf")])
    def test_identify_metabolites_rejects_a_bad_m_z_and_names_it(self, bad: float) -> None:
        with pytest.raises(ValueError, match=r"observed_mz must contain only positive finite m/z values; element 1 is"):
            identify_metabolites(np.array([180.0634, bad, 90.0317]), DATABASE)

    def test_a_negative_m_z_no_longer_matches_the_whole_database(self) -> None:
        with pytest.raises(ValueError):
            identify_metabolites(np.array([-180.0634]), DATABASE, ppm_tolerance=10.0)

    @pytest.mark.parametrize("bad", [-5.0, 0.0, float("nan"), float("inf")])
    def test_identify_with_adducts_rejects_a_bad_m_z(self, bad: float) -> None:
        with pytest.raises(ValueError, match=r"observed_mz must contain only positive finite m/z values"):
            identify_with_adducts(np.array([bad]), DATABASE)

    def test_the_count_of_bad_values_is_reported(self) -> None:
        with pytest.raises(ValueError, match=r"\(3 invalid value\(s\)\)"):
            identify_metabolites(np.array([1.0, -1.0, 0.0, float("nan")]), DATABASE)

    def test_valid_input_in_any_array_like_is_unchanged(self) -> None:
        reference = identify_metabolites(np.array([180.0634, 90.0317]), DATABASE)

        assert identify_metabolites([180.0634, 90.0317], DATABASE) == reference
        assert [m.matched_name for m in reference[1]] == ["lactate"]

    def test_an_empty_array_gives_no_results(self) -> None:
        assert identify_metabolites(np.array([]), DATABASE) == []
        assert identify_with_adducts(np.array([]), DATABASE) == []


def _knn_oracle(data: np.ndarray, k: int = 5) -> np.ndarray:
    """The documented rule written out directly: fill a gap with the mean of its k most correlated metabolites."""
    filled = data.copy()
    for i in range(data.shape[0]):
        missing = np.isnan(data[i])
        if not missing.any():
            continue
        scored: list[tuple[float, int]] = []
        for j in range(data.shape[0]):
            both = ~(np.isnan(data[i]) | np.isnan(data[j]))
            if j == i or both.sum() < 2:
                continue  # no usable overlap: not a neighbor
            scored.append((float(np.corrcoef(data[i, both], data[j, both])[0, 1]), j))
        neighbors = [j for _, j in sorted(scored)[-min(k, data.shape[0] - 1) :]]
        for col in np.flatnonzero(missing):
            values = data[neighbors, col] if neighbors else np.array([])
            values = values[~np.isnan(values)]
            filled[i, col] = values.mean() if values.size else 0.0
    return filled


class TestKnnImputationUsesOnlyRealNeighbors:
    def test_a_metabolite_with_one_observed_value_has_no_neighbors_and_is_not_filled_from_strangers(self) -> None:
        rng = np.random.default_rng(3)
        data = rng.uniform(10.0, 50.0, size=(7, 6))
        data[0, 1:] = np.nan  # a single observed value: no row shares two samples with it

        result = missing_value_imputation(data, method="knn")

        assert result[0, 0] == data[0, 0]
        np.testing.assert_array_equal(
            result[0, 1:], 0.0
        )  # nothing to impute from; the old code used other rows' values
        np.testing.assert_array_equal(result[1:], data[1:])

    def test_with_enough_overlap_the_result_follows_the_k_nearest_row_rule(self) -> None:
        # only one metabolite has gaps, so no earlier imputed value can feed a later one
        rng = np.random.default_rng(4)
        data = rng.uniform(10.0, 50.0, size=(9, 8))
        data[3, [1, 4, 6]] = np.nan

        result = missing_value_imputation(data, method="knn")

        np.testing.assert_allclose(result, _knn_oracle(data), rtol=1e-12)
        assert not np.allclose(result[3], data[3], equal_nan=True)

    def test_fewer_valid_neighbors_than_k_means_only_the_valid_ones_are_averaged(self) -> None:
        rng = np.random.default_rng(5)
        data = rng.uniform(10.0, 50.0, size=(7, 6))
        data[0, 3:] = np.nan  # row 0 is observed in samples 0-2
        for stranger in (1, 4, 5, 6):  # these rows overlap row 0 in one sample only: no correlation, no neighbor
            data[stranger, :2] = np.nan
        valid_neighbors = [2, 3]  # complete rows: three shared samples with row 0

        result = missing_value_imputation(data, method="knn")

        expected = data[valid_neighbors][:, 3:].mean(axis=0)  # the old code averaged the strangers in as well
        np.testing.assert_allclose(result[0, 3:], expected, rtol=1e-12)
        np.testing.assert_array_equal(result[0, :3], data[0, :3])

    def test_unknown_method_still_raises(self) -> None:
        with pytest.raises(ValueError, match=r"Unknown imputation method"):
            missing_value_imputation(np.ones((2, 2)), method="typo")

    def test_a_metabolite_with_no_observed_value_stays_at_the_missing_marker(self) -> None:
        data = np.array([[np.nan, np.nan, np.nan], [4.0, 5.0, 6.0]])

        for method in ("min_half", "median", "knn"):
            result: Any = missing_value_imputation(data, method=method)
            np.testing.assert_array_equal(result[0], 0.0)
            np.testing.assert_array_equal(result[1], data[1])


class TestPathwayActivityNeedsMeasuredMetabolites:
    PATHWAYS = {"measured": ["a", "b"], "unmeasured": ["x", "y"], "empty": []}

    def test_a_pathway_with_no_measured_metabolite_is_not_given_a_made_up_score(self) -> None:
        # with min_coverage = 0 the old code reported the unmeasured pathways with activity_score 0.0
        results = pathway_activity_scoring({"a": 2.0}, self.PATHWAYS, min_coverage=0.0)

        assert [r.pathway_name for r in results] == ["measured"]
        assert results[0].activity_score == 2.0 and results[0].n_measured == 1 and results[0].n_total == 2

    def test_the_default_coverage_filter_is_unchanged(self) -> None:
        results = pathway_activity_scoring({"a": 2.0, "b": 4.0}, self.PATHWAYS)

        assert [r.pathway_name for r in results] == ["measured"]
        assert results[0].activity_score == 3.0

    def test_nothing_measured_gives_nothing(self) -> None:
        assert pathway_activity_scoring({}, self.PATHWAYS, min_coverage=0.0) == []
