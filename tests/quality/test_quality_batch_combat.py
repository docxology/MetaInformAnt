"""Parametric empirical Bayes ComBat (Johnson, Li & Rabinovic 2007) in ``correct_batch_combat``.

Regression tests for ``metainformant.quality.batch.detection.correct_batch_combat``. The old code

* recomputed the grand standard deviation inside the loop over batches, from data that earlier iterations had
  already corrected, so the result depended on the order in which batches were processed (it processed them in
  sorted-label order, so merely renaming the batches changed the values), and
* "shrunk" every scale by a constant 0.5 instead of estimating hyperparameters from the data and iterating to the
  empirical Bayes posterior means.

Oracle. Neither ``pycombat`` nor ``inmoose`` is installed and there is no network, so the reference is
``_hand_combat`` below: a plain-Python transcription of the paper's equations (loops over features, batches and
samples; a scalar fixed-point iteration per feature; no code shared with the module). Because a transcription error
could be made twice, the tests also check properties of the method that such an error would break: the two
posterior-mean equations hold at the returned values, the posterior means lie between the raw estimate and the prior
mean, and on data simulated from the ComBat model the empirical Bayes estimates are closer to the true parameters
than the raw ones.

Real computation; no test doubles.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.quality.batch import detection
from metainformant.quality.batch.detection import correct_batch_combat


def _hand_combat(data: np.ndarray, labels: list[str]) -> np.ndarray:
    """ComBat, parametric, no covariates, written one scalar at a time straight from the paper."""
    n_samples, n_features = data.shape
    batches: dict[str, list[int]] = {}
    for sample, label in enumerate(labels):
        batches.setdefault(label, []).append(sample)

    # standardise once: alpha_g = grand mean, sigma_g^2 = pooled within-batch variance divided by N
    alpha = [sum(data[j][g] for j in range(n_samples)) / n_samples for g in range(n_features)]
    sigma_sq = []
    for g in range(n_features):
        residual_ss = 0.0
        for members in batches.values():
            mean = sum(data[j][g] for j in members) / len(members)
            residual_ss += sum((data[j][g] - mean) ** 2 for j in members)
        sigma_sq.append(residual_ss / n_samples)
    z = [[(data[j][g] - alpha[g]) / math.sqrt(sigma_sq[g]) for g in range(n_features)] for j in range(n_samples)]

    out = [[0.0] * n_features for _ in range(n_samples)]
    for members in batches.values():
        n = len(members)
        gamma_hat = [sum(z[j][g] for j in members) / n for g in range(n_features)]
        delta_hat_sq = [sum((z[j][g] - gamma_hat[g]) ** 2 for j in members) / (n - 1) for g in range(n_features)]
        # hyperparameters by the method of moments over features
        gamma_bar = sum(gamma_hat) / n_features
        tau_sq = sum((v - gamma_bar) ** 2 for v in gamma_hat) / (n_features - 1)
        mean_d = sum(delta_hat_sq) / n_features
        var_d = sum((v - mean_d) ** 2 for v in delta_hat_sq) / (n_features - 1)
        lam = (2 * var_d + mean_d**2) / var_d
        theta = (mean_d * var_d + mean_d**3) / var_d
        for g in range(n_features):
            gamma_star, delta_star_sq = gamma_hat[g], delta_hat_sq[g]
            for _ in range(10_000):
                new_gamma = (n * tau_sq * gamma_hat[g] + delta_star_sq * gamma_bar) / (n * tau_sq + delta_star_sq)
                new_delta = (theta + 0.5 * sum((z[j][g] - new_gamma) ** 2 for j in members)) / (n / 2 + lam - 1)
                converged = abs(new_gamma - gamma_star) < 1e-14 and abs(new_delta - delta_star_sq) < 1e-14
                gamma_star, delta_star_sq = new_gamma, new_delta
                if converged:
                    break
            for j in members:
                out[j][g] = math.sqrt(sigma_sq[g]) / math.sqrt(delta_star_sq) * (z[j][g] - gamma_star) + alpha[g]
    return np.array(out)


def _batched_data(
    sizes: tuple[int, ...],
    shifts: tuple[float, ...],
    scales: tuple[float, ...],
    n_features: int,
    seed: int,
) -> tuple[np.ndarray, list[str]]:
    """Batches with their own location and scale on top of per-feature baselines."""
    rng = np.random.default_rng(seed)
    baseline = rng.normal(0.0, 3.0, size=n_features)
    blocks, labels = [], []
    for index, (size, shift, scale) in enumerate(zip(sizes, shifts, scales)):
        blocks.append(baseline + rng.normal(shift, scale, size=(size, n_features)))
        labels += [f"batch{index}"] * size
    return np.vstack(blocks), labels


def _median_feature_sd(block: np.ndarray) -> float:
    return float(np.median(block.std(axis=0, ddof=1)))


class TestAgainstHandImplementation:
    def test_three_unequal_batches(self) -> None:
        data, labels = _batched_data((5, 4, 6), (0.0, 2.0, -1.5), (1.0, 2.0, 0.5), n_features=25, seed=3)

        np.testing.assert_allclose(correct_batch_combat(data, labels), _hand_combat(data, labels), rtol=1e-8, atol=1e-9)

    def test_two_batches(self) -> None:
        data, labels = _batched_data((7, 9), (0.0, 3.0), (1.0, 1.5), n_features=18, seed=11)

        np.testing.assert_allclose(correct_batch_combat(data, labels), _hand_combat(data, labels), rtol=1e-8, atol=1e-9)

    def test_many_tiny_batches(self) -> None:
        data, labels = _batched_data((2, 3, 2, 4), (0.0, 1.0, -2.0, 0.5), (1.0, 0.7, 1.8, 1.0), n_features=30, seed=17)

        np.testing.assert_allclose(correct_batch_combat(data, labels), _hand_combat(data, labels), rtol=1e-8, atol=1e-9)


class TestBatchOrderInvariance:
    DATA, LABELS = _batched_data((6, 8, 5), (0.0, 2.0, -1.5), (1.0, 2.5, 0.6), n_features=40, seed=11)

    def test_renaming_batches_so_that_they_sort_in_another_order_changes_nothing(self) -> None:
        # the old code processed batches in sorted-label order: this reverses that order for the same grouping
        renamed = {"batch0": "z", "batch1": "y", "batch2": "x"}
        relabelled = [renamed[label] for label in self.LABELS]

        np.testing.assert_allclose(
            correct_batch_combat(self.DATA, relabelled),
            correct_batch_combat(self.DATA, self.LABELS),
            rtol=0,
            atol=1e-12,
        )

    def test_shuffling_the_samples_only_shuffles_the_result(self) -> None:
        permutation = np.random.default_rng(0).permutation(len(self.LABELS))
        shuffled_labels = [self.LABELS[i] for i in permutation]

        np.testing.assert_allclose(
            correct_batch_combat(self.DATA[permutation], shuffled_labels),
            correct_batch_combat(self.DATA, self.LABELS)[permutation],
            rtol=0,
            atol=1e-12,
        )

    def test_stacking_the_batches_in_another_order_gives_each_sample_the_same_value(self) -> None:
        by_batch = {label: np.flatnonzero(np.array(self.LABELS) == label) for label in ("batch0", "batch1", "batch2")}
        order = np.concatenate([by_batch["batch2"], by_batch["batch0"], by_batch["batch1"]])

        np.testing.assert_allclose(
            correct_batch_combat(self.DATA[order], [self.LABELS[i] for i in order]),
            correct_batch_combat(self.DATA, self.LABELS)[order],
            rtol=0,
            atol=1e-12,
        )


class TestEmpiricalBayes:
    @staticmethod
    def _standardised_batch(n: int = 6, n_features: int = 60, seed: int = 4) -> np.ndarray:
        rng = np.random.default_rng(seed)
        gamma = rng.normal(0.4, 0.5, n_features)
        delta = np.sqrt(2.0 / rng.gamma(5.0, 1.0, n_features) * 4.0)
        return gamma + delta * rng.standard_normal((n, n_features))

    def test_the_returned_posterior_means_satisfy_both_posterior_equations(self) -> None:
        z = self._standardised_batch()
        n = z.shape[0]
        gamma_star, delta_star_sq = detection._combat_posterior(z)

        # hyperparameters from their definitions (method of moments over features)
        gamma_hat = np.array([np.mean(z[:, g]) for g in range(z.shape[1])])
        delta_hat_sq = np.array([np.var(z[:, g], ddof=1) for g in range(z.shape[1])])
        gamma_bar, tau_sq = np.mean(gamma_hat), np.var(gamma_hat, ddof=1)
        mean_d, var_d = np.mean(delta_hat_sq), np.var(delta_hat_sq, ddof=1)
        lam = (2 * var_d + mean_d**2) / var_d
        theta = (mean_d * var_d + mean_d**3) / var_d

        expected_gamma = (n * tau_sq * gamma_hat + delta_star_sq * gamma_bar) / (n * tau_sq + delta_star_sq)
        expected_delta = (theta + 0.5 * ((z - gamma_star) ** 2).sum(axis=0)) / (n / 2 + lam - 1)
        np.testing.assert_allclose(gamma_star, expected_gamma, rtol=0, atol=1e-8)
        np.testing.assert_allclose(delta_star_sq, expected_delta, rtol=1e-8, atol=0)

    def test_the_location_posterior_lies_between_the_raw_estimate_and_the_prior_mean(self) -> None:
        z = self._standardised_batch()
        gamma_star, _ = detection._combat_posterior(z)
        gamma_hat = z.mean(axis=0)
        gamma_bar = gamma_hat.mean()

        # a convex combination of the two, so it moves towards the prior mean and never past the raw estimate
        assert np.all((gamma_star - gamma_bar) * (gamma_hat - gamma_bar) >= -1e-12)
        assert np.all(np.abs(gamma_star - gamma_bar) <= np.abs(gamma_hat - gamma_bar) + 1e-12)
        assert np.any(np.abs(gamma_star - gamma_hat) > 1e-3)  # and it really does shrink

    def test_the_estimates_are_closer_to_the_truth_than_the_raw_batch_statistics(self) -> None:
        # Data drawn from the ComBat model with known parameters: gamma_g ~ N(0.5, 0.4^2) and
        # delta_g^2 ~ InverseGamma(6, 8) (mean 1.6), four samples per feature. The raw batch mean has mean squared
        # error delta^2 / n = 0.4 and the raw variance estimate (3 degrees of freedom) is very noisy. ComBat's moment
        # estimate of the prior spread is the spread of the raw estimates, which still contains their sampling noise,
        # so the shrinkage is milder than with the true prior; measured on this seed the empirical Bayes error is 0.42
        # (location) and 0.27 (scale) of the raw error. The assertions only require it to be at least halved.
        rng = np.random.default_rng(2024)
        n, n_features = 4, 4000
        gamma = rng.normal(0.5, 0.4, n_features)
        delta_sq = 8.0 / rng.gamma(6.0, 1.0, n_features)
        z = gamma + np.sqrt(delta_sq) * rng.standard_normal((n, n_features))

        gamma_star, delta_star_sq = detection._combat_posterior(z)
        gamma_hat, delta_hat_sq = z.mean(axis=0), z.var(axis=0, ddof=1)

        assert np.mean((gamma_star - gamma) ** 2) < 0.5 * np.mean((gamma_hat - gamma) ** 2)
        assert np.mean((delta_star_sq - delta_sq) ** 2) < 0.5 * np.mean((delta_hat_sq - delta_sq) ** 2)

    def test_the_scale_adjustment_is_estimated_from_the_data_not_a_fixed_fraction(self) -> None:
        # one batch is three times as noisy as the other two; location is identical everywhere
        data, labels = _batched_data((8, 8, 8), (0.0, 0.0, 0.0), (1.0, 1.0, 3.0), n_features=200, seed=0)
        blocks = {name: slice(8 * i, 8 * i + 8) for i, name in enumerate(("batch0", "batch1", "batch2"))}
        corrected = correct_batch_combat(data, labels)

        before = [_median_feature_sd(data[s]) for s in blocks.values()]
        after = [_median_feature_sd(corrected[s]) for s in blocks.values()]

        assert max(before) / min(before) > 2.5  # the batch effect is there ...
        assert max(after) / min(after) < 1.25  # ... and the per-batch spreads are now alike (a constant 0.5 left ~1.7)

    def test_a_location_effect_is_removed_without_touching_the_scale(self) -> None:
        data, labels = _batched_data((8, 8, 8), (0.0, 4.0, -3.0), (1.0, 1.0, 1.0), n_features=200, seed=1)
        corrected = correct_batch_combat(data, labels)
        blocks = [slice(0, 8), slice(8, 16), slice(16, 24)]

        def mean_offset(matrix: np.ndarray) -> float:
            """Average over features of how far the batch means are from the grand mean."""
            grand = matrix.mean(axis=0)
            return float(np.mean([np.abs(matrix[s].mean(axis=0) - grand).mean() for s in blocks]))

        assert mean_offset(corrected) < 0.1 * mean_offset(data)
        for s in blocks:
            assert _median_feature_sd(corrected[s]) == pytest.approx(_median_feature_sd(data[s]), rel=0.2)

    def test_features_that_all_share_one_variance_and_one_location_degenerate_to_standardisation(self) -> None:
        # identical features: the priors are point masses (tau^2 = 0 and no spread of delta^2), so the posterior
        # means equal the raw batch statistics and each batch is just centred and scaled to the pooled spread
        rng = np.random.default_rng(8)
        column = np.concatenate([rng.normal(0.0, 1.0, 5), rng.normal(3.0, 2.0, 6)])
        data = np.tile(column[:, None], (1, 4))
        labels = ["a"] * 5 + ["b"] * 6
        groups = [slice(0, 5), slice(5, 11)]

        corrected = correct_batch_combat(data, labels)

        grand = column.mean()
        pooled_sd = math.sqrt(sum(((column[s] - column[s].mean()) ** 2).sum() for s in groups) / len(column))
        expected = np.empty_like(column)
        for s in groups:
            expected[s] = (column[s] - column[s].mean()) / column[s].std(ddof=1) * pooled_sd + grand
        np.testing.assert_allclose(corrected, np.tile(expected[:, None], (1, 4)), rtol=1e-9, atol=1e-10)


class TestEdgeCases:
    def test_a_single_batch_is_returned_unchanged_as_a_copy(self) -> None:
        data = np.random.default_rng(5).normal(size=(8, 4))

        corrected = correct_batch_combat(data, ["A"] * 8)

        np.testing.assert_array_equal(corrected, data)
        assert corrected is not data

    def test_the_input_is_not_modified(self) -> None:
        data, labels = _batched_data((5, 5), (0.0, 2.0), (1.0, 1.0), n_features=10, seed=2)
        snapshot = data.copy()

        correct_batch_combat(data, labels)

        np.testing.assert_array_equal(data, snapshot)

    def test_features_without_variation_inside_a_batch_are_left_unchanged(self) -> None:
        data, labels = _batched_data((5, 5, 5), (0.0, 2.0, -1.0), (1.0, 1.0, 1.0), n_features=8, seed=6)
        data[5:10, 0] = 7.5  # constant inside batch1 only: its variance there cannot be estimated
        data[:, 1] = 3.0  # constant everywhere

        corrected = correct_batch_combat(data, labels)

        np.testing.assert_array_equal(corrected[:, 0], data[:, 0])
        np.testing.assert_array_equal(corrected[:, 1], data[:, 1])
        assert np.all(np.isfinite(corrected))
        assert not np.allclose(corrected[:, 2:], data[:, 2:])  # the estimable features were corrected

    def test_all_constant_data_comes_back_unchanged_and_finite(self) -> None:
        data = np.ones((10, 3))

        corrected = correct_batch_combat(data, ["A"] * 5 + ["B"] * 5)

        np.testing.assert_array_equal(corrected, data)

    @pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
    def test_non_finite_values_are_rejected_not_spread_over_every_feature(self, bad: float) -> None:
        data, labels = _batched_data((4, 4), (0.0, 1.0), (1.0, 1.0), n_features=6, seed=1)
        data[2, 3] = bad

        with pytest.raises(ValueError, match=r"NaN or infinite"):
            correct_batch_combat(data, labels)

    def test_a_label_list_of_the_wrong_length_is_rejected(self) -> None:
        with pytest.raises(ValueError, match=r"batch_labels has 3 entries for 4 samples"):
            correct_batch_combat(np.ones((4, 3)), ["A", "A", "B"])

    def test_one_dimensional_data_is_rejected(self) -> None:
        with pytest.raises(ValueError, match=r"2D samples x features"):
            correct_batch_combat(np.ones(6), ["A"] * 3 + ["B"] * 3)

    def test_a_batch_with_one_sample_is_rejected_because_its_variance_is_not_estimable(self) -> None:
        data, labels = _batched_data((4, 1, 4), (0.0, 1.0, 2.0), (1.0, 1.0, 1.0), n_features=6, seed=3)

        with pytest.raises(ValueError, match=r"at least 2 samples.*\['batch1'\]"):
            correct_batch_combat(data, labels)

    def test_a_single_estimable_feature_cannot_fit_the_priors(self) -> None:
        data, labels = _batched_data((4, 4), (0.0, 1.0), (1.0, 1.0), n_features=5, seed=3)
        data[:, 1:] = 2.0  # only feature 0 varies

        with pytest.raises(ValueError, match=r"at least 2 features"):
            correct_batch_combat(data, labels)

    def test_two_estimable_features_are_enough(self) -> None:
        data, labels = _batched_data((4, 4), (0.0, 1.0), (1.0, 1.0), n_features=5, seed=3)
        data[:, 2:] = 2.0

        corrected = correct_batch_combat(data, labels)

        assert np.all(np.isfinite(corrected))
        np.testing.assert_array_equal(corrected[:, 2:], data[:, 2:])

    def test_batches_of_two_samples_are_corrected_without_numerical_warnings(self) -> None:
        # pytest turns every warning into an error, so a 0/0 or an overflow inside the estimation would fail here
        data, labels = _batched_data((2, 2, 2), (0.0, 1.0, -1.0), (1.0, 1.0, 1.0), n_features=20, seed=9)

        assert np.all(np.isfinite(correct_batch_combat(data, labels)))

    def test_non_string_hashable_labels_work_like_strings(self) -> None:
        data, labels = _batched_data((5, 5), (0.0, 2.0), (1.0, 1.5), n_features=12, seed=4)
        numeric = [0 if label == "batch0" else 1 for label in labels]

        np.testing.assert_array_equal(correct_batch_combat(data, numeric), correct_batch_combat(data, labels))
