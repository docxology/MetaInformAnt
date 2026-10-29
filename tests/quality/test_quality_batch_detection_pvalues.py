"""``detect_batch_effects`` counts features with exact F-test p-values and never turns a failure into p = 1.

Regression tests for ``metainformant.quality.batch.detection.detect_batch_effects``:

* the per-feature p-value was ``1.0 - f_dist.cdf(F, ...)``. The cdf rounds to 1.0 for every upper tail below about
  1e-16, so every p-value from 1e-17 to 1e-300 came out as exactly 0.0 and passed any ``alpha`` below it, however
  small (a Bonferroni-style ``alpha`` for a genome-wide feature count is easily 1e-20). ``sf`` is exact there;
* a feature that does not vary inside any batch but differs between batches (the strongest possible batch effect)
  got ``F = 0`` and so ``p = 1.0``, "not significant", because the within-batch mean square was exactly zero;
* a NaN made its feature's F-statistic 0 (p = 1.0) and, through a silent fallback after the SVD failed, turned the
  PCA-space silhouette into a silhouette of the raw features, so a dataset with one missing value was reported as
  having no batch effect;
* a label list of the wrong length, or labels that were sample identifiers, raised an obscure ``IndexError`` or
  produced a report built on no within-batch variation.

The oracles are scipy's own ``f_oneway`` and ``f.sf``. Real computation; no test doubles (real-implementation
policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy import stats

from metainformant.quality.batch.detection import detect_batch_effects

N_PER_BATCH = 20
LABELS = ["A"] * N_PER_BATCH + ["B"] * N_PER_BATCH
# alternating -1/+1 inside every batch: mean 0, so a shift of d between the batches gives F = 9.5 d^2 on (1, 38) df
PATTERN = np.tile([-1.0, 1.0], N_PER_BATCH // 2)


def _shifted_features(shifts: list[float]) -> np.ndarray:
    return np.column_stack([np.concatenate([PATTERN, PATTERN + d]) for d in shifts])


def _reference_p_values(data: np.ndarray) -> np.ndarray:
    return np.array([stats.f_oneway(data[:N_PER_BATCH, j], data[N_PER_BATCH:, j]).pvalue for j in range(data.shape[1])])


class TestExactTailProbabilities:
    SHIFTS = [6.0, 10.0, 15.0, 20.0]  # p = 1e-20, 2e-28, 5e-35, 1e-39

    def test_the_scenario_is_one_where_one_minus_cdf_cannot_tell_the_features_apart(self) -> None:
        data = _shifted_features(self.SHIFTS)
        reference = _reference_p_values(data)
        f_stats = 9.5 * np.array(self.SHIFTS) ** 2

        assert np.all(reference < 1e-19)
        assert np.all(1.0 - stats.f.cdf(f_stats, 1, 2 * N_PER_BATCH - 2) == 0.0)  # what the old code compared
        assert reference == pytest.approx(stats.f.sf(f_stats, 1, 2 * N_PER_BATCH - 2), rel=1e-9)

    def test_a_tiny_alpha_separates_features_by_their_true_p_values(self) -> None:
        data = _shifted_features(self.SHIFTS)
        alpha = 1e-30
        expected = int((_reference_p_values(data) < alpha).sum())

        report = detect_batch_effects(data, LABELS, alpha=alpha)

        assert expected == 2  # the shifts 15 and 20 are below 1e-30, 6 and 10 are not
        assert report.n_significant_features == expected  # 1 - cdf counted all four

    @pytest.mark.parametrize("alpha", [1e-18, 1e-25, 1e-33, 1e-37, 1e-45])
    def test_the_count_follows_the_exact_p_values_at_every_alpha(self, alpha: float) -> None:
        data = _shifted_features(self.SHIFTS)

        report = detect_batch_effects(data, LABELS, alpha=alpha)

        assert report.n_significant_features == int((_reference_p_values(data) < alpha).sum())

    @pytest.mark.parametrize("alpha", [0.5, 0.05, 0.01, 1e-4, 1e-9])
    def test_ordinary_alphas_agree_with_scipy_on_random_data(self, alpha: float) -> None:
        rng = np.random.default_rng(12)
        data = rng.normal(size=(2 * N_PER_BATCH, 60))
        data[N_PER_BATCH:, :25] += np.linspace(0.2, 3.0, 25)  # a graded batch effect on the first features
        expected = int((_reference_p_values(data) < alpha).sum())

        assert detect_batch_effects(data, LABELS, alpha=alpha).n_significant_features == expected

    def test_three_batches_use_the_right_degrees_of_freedom(self) -> None:
        rng = np.random.default_rng(5)
        data = rng.normal(size=(30, 40))
        data[10:20, :10] += 1.5
        data[20:, 5:15] -= 1.0
        labels = ["x"] * 10 + ["y"] * 10 + ["z"] * 10
        p_ref = np.array(
            [stats.f_oneway(data[:10, j], data[10:20, j], data[20:, j]).pvalue for j in range(data.shape[1])]
        )

        for alpha in (0.05, 0.001):
            assert detect_batch_effects(data, labels, alpha=alpha).n_significant_features == int((p_ref < alpha).sum())

    def test_unequal_batch_sizes_use_the_right_degrees_of_freedom(self) -> None:
        rng = np.random.default_rng(8)
        data = rng.normal(size=(23, 30))
        data[:5, :12] += 2.0
        labels = ["small"] * 5 + ["large"] * 18
        p_ref = np.array([stats.f_oneway(data[:5, j], data[5:, j]).pvalue for j in range(data.shape[1])])

        for alpha in (0.05, 0.001):
            assert detect_batch_effects(data, labels, alpha=alpha).n_significant_features == int((p_ref < alpha).sum())


class TestDegenerateFeatures:
    def test_a_feature_constant_inside_each_batch_but_different_between_them_is_significant(self) -> None:
        noise = np.random.default_rng(3).normal(size=(2 * N_PER_BATCH, 1))
        perfect = np.concatenate([np.zeros(N_PER_BATCH), np.full(N_PER_BATCH, 5.0)])[:, None]
        data = np.hstack([perfect, noise])

        report = detect_batch_effects(data, LABELS, alpha=0.001)

        assert report.n_significant_features == 1  # the perfectly separated feature; its F used to be set to 0

    def test_a_feature_that_is_constant_everywhere_is_not_associated_with_the_batch(self) -> None:
        data = np.hstack([np.full((2 * N_PER_BATCH, 1), 3.0), _shifted_features([8.0])])

        report = detect_batch_effects(data, LABELS, alpha=0.05)

        assert report.n_significant_features == 1

    def test_no_numpy_warning_is_raised_for_zero_within_batch_variance(self) -> None:
        # pytest turns RuntimeWarning (divide by zero, invalid value) into an error, so this fails on a bare 0/0
        data = np.zeros((10, 3))

        report = detect_batch_effects(data, ["a"] * 5 + ["b"] * 5)

        assert report.n_significant_features == 0
        assert report.pvca_variance == {"batch": 0.0, "residual": 1.0}

    def test_a_single_batch_has_nothing_to_detect(self) -> None:
        data = np.random.default_rng(1).normal(size=(10, 4))

        report = detect_batch_effects(data, ["only"] * 10)

        assert (report.n_batches, report.n_significant_features, report.silhouette_score) == (1, 0, 0.0)
        assert report.severity == "low"

    def test_a_batch_of_one_sample_is_allowed_when_another_batch_has_replicates(self) -> None:
        data = np.random.default_rng(2).normal(size=(9, 5))

        report = detect_batch_effects(data, ["a"] * 4 + ["b"] * 4 + ["c"])

        assert report.n_batches == 3


class TestInputsThatCannotBeAnalysed:
    @pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
    def test_a_non_finite_value_raises_instead_of_looking_like_no_batch_effect(self, bad: float) -> None:
        data = _shifted_features([8.0, 8.0])
        data[3, 1] = bad

        with pytest.raises(ValueError, match=r"NaN or infinite"):
            detect_batch_effects(data, LABELS)

    def test_the_number_of_labels_must_match_the_number_of_samples(self) -> None:
        data = _shifted_features([8.0])

        with pytest.raises(ValueError, match=r"batch_labels has 39 entries for 40 samples"):
            detect_batch_effects(data, LABELS[:-1])
        with pytest.raises(ValueError, match=r"batch_labels has 41 entries for 40 samples"):
            detect_batch_effects(data, LABELS + ["B"])

    def test_every_sample_in_its_own_batch_is_refused(self) -> None:
        data = np.random.default_rng(4).normal(size=(6, 3))

        with pytest.raises(ValueError, match=r"batches of their own"):
            detect_batch_effects(data, [f"sample{i}" for i in range(6)])

    def test_a_one_dimensional_array_raises(self) -> None:
        with pytest.raises(ValueError, match=r"2D"):
            detect_batch_effects(np.arange(10.0), ["a"] * 5 + ["b"] * 5)

    def test_an_empty_array_raises(self) -> None:
        with pytest.raises(ValueError, match=r"no samples"):
            detect_batch_effects(np.empty((0, 3)), [])

    @pytest.mark.parametrize("bad", [0.0, 1.0, -0.05, 1.5, float("nan")])
    def test_alpha_must_lie_strictly_between_zero_and_one(self, bad: float) -> None:
        data = _shifted_features([8.0])

        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            detect_batch_effects(data, LABELS, alpha=bad)

    def test_a_list_of_lists_is_accepted_like_an_array(self) -> None:
        data = _shifted_features([8.0, 2.0]).tolist()

        report: Any = detect_batch_effects(data, LABELS)  # type: ignore[arg-type]

        assert report.n_samples == 2 * N_PER_BATCH
