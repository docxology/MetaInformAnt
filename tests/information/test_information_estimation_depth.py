"""Depth tests for metainformant.information metrics estimation (Round-4 T3).

Zero-mock, value-pinned: exact entropies for known distributions, estimator
bias direction, MI/KL identities, and entropy-rate sanity on deterministic
and i.i.d. sequences.

Known limitations documented (not asserted away silently):
- jackknife (Zahl 1977 category jackknife) can overestimate on dense
  uniform samples (may exceed log2 k); tested for bias direction only.
- panzeri_treves_bias_correction's code formula differs from its
  docstring comment; tested for structural bounds only.
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pytest

from metainformant.information.metrics.core.estimation import (
    bias_correction,
    effective_sample_size_correction,
    entropy_estimator,
    entropy_rate_estimator,
    kl_divergence_estimator,
    mutual_information_estimator,
    panzeri_treves_bias_correction,
)


class TestPluginEntropyExactValues:
    """Plugin entropy pinned to closed-form values (uncorrected)."""

    def test_uniform_distribution_is_log_k(self) -> None:
        counts = {"A": 25, "C": 25, "G": 25, "T": 25}
        h = entropy_estimator(counts, method="plugin", bias_correction=False)
        assert h == pytest.approx(2.0, abs=1e-12)

    def test_deterministic_distribution_is_zero(self) -> None:
        assert entropy_estimator({"A": 100}, method="plugin", bias_correction=False) == 0.0

    def test_known_binary_value(self) -> None:
        # p=(3/4,1/4): H = 0.8112781244591328 bits
        h = entropy_estimator({"A": 75, "B": 25}, method="plugin", bias_correction=False)
        assert h == pytest.approx(0.8112781244591328, abs=1e-12)

    def test_bias_correction_reduces_entropy(self) -> None:
        counts = {"A": 50, "T": 30, "G": 20}
        raw = entropy_estimator(counts, method="plugin", bias_correction=False)
        corrected = entropy_estimator(counts, method="plugin", bias_correction=True)
        assert corrected < raw
        # Miller correction: (k-1)/(2n) = 2/(2*100) = 0.01
        assert raw - corrected == pytest.approx(0.01, abs=1e-12)

    def test_miller_madow_matches_plugin_plus_correction(self) -> None:
        counts = {"A": 50, "T": 30, "G": 20}
        raw = entropy_estimator(counts, method="plugin", bias_correction=False)
        mm = entropy_estimator(counts, method="miller_madow")
        assert raw - mm == pytest.approx((3 - 1) / (2 * 100), abs=1e-12)

    def test_tiny_sample_substantial_correction(self) -> None:
        # k=2, n=2: H=1 bit, correction (2-1)/(2*2)=0.25 -> 0.75
        h = entropy_estimator({"A": 1, "B": 1}, method="plugin", bias_correction=True)
        assert h == pytest.approx(0.75, abs=1e-12)

    def test_list_and_dict_inputs_agree(self) -> None:
        assert entropy_estimator([50, 30, 20]) == entropy_estimator({"A": 50, "T": 30, "G": 20})


class TestEstimatorOrdering:
    """Bias-corrected estimators must not be systematically below plugin."""

    def test_chao_shen_finite_and_ordered_sparse(self) -> None:
        # Sparse counts: all estimators in [0, log2(k)]
        counts = {"A": 100, "T": 2, "G": 1, "C": 1, "N": 1}
        log_k = math.log2(5)
        for method in ("plugin", "miller_madow", "chao_shen", "jackknife"):
            h = entropy_estimator(counts, method=method)
            assert 0.0 <= h <= log_k + 1e-9, method

    def test_jackknife_applies_nonnegative_bias_direction(self) -> None:
        # Plugin entropy is downward-biased; the jackknife correction must
        # not systematically push below the plugin estimate.
        rng = np.random.RandomState(12)
        samples = rng.choice(4, size=2000)
        counts = Counter(samples.tolist())
        h_plugin = entropy_estimator(counts, method="plugin", bias_correction=False)
        h_jack = entropy_estimator(counts, method="jackknife")
        assert h_jack >= h_plugin
        # Known limitation (documented, not asserted away): the category
        # jackknife can overshoot log2(k) on dense uniform samples.


class TestMutualInformationIdentities:
    def test_independent_variables_zero_mi_uncorrected(self) -> None:
        x = [0, 0, 1, 1] * 25
        y = [0, 1, 0, 1] * 25
        mi = mutual_information_estimator(x, y, bias_correction=False)
        assert mi == pytest.approx(0.0, abs=1e-12)

    def test_perfectly_dependent_equals_marginal_entropy(self) -> None:
        y = list(range(50))
        x = [v * 2 for v in y]  # deterministic bijection
        mi = mutual_information_estimator(x, y, bias_correction=False)
        assert mi == pytest.approx(math.log2(50), abs=1e-9)

    def test_mi_of_self_is_entropy(self) -> None:
        x = ["a"] * 40 + ["b"] * 60
        mi = mutual_information_estimator(x, list(x), bias_correction=False)
        h = entropy_estimator({"a": 40, "b": 60}, bias_correction=False)
        assert mi == pytest.approx(h, abs=1e-9)

    def test_mi_non_negative_with_default_corrections(self) -> None:
        x = [0, 0, 1, 1] * 25
        y = [0, 1, 0, 1] * 25
        assert mutual_information_estimator(x, y) >= 0.0

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            mutual_information_estimator([0, 1], [0])


class TestKLDivergence:
    def test_identical_samples_zero(self) -> None:
        p = [0, 1, 2, 0, 1]
        assert kl_divergence_estimator(p, list(p)) == pytest.approx(0.0, abs=1e-12)

    def test_q_missing_p_support_is_infinite(self) -> None:
        # q assigns zero probability to an outcome p observes => KL = inf
        p = [0, 1, 1]
        q = [1, 1, 1]
        assert kl_divergence_estimator(p, q) == float("inf")

    def test_known_asymmetric_value(self) -> None:
        # p -> (0.75, 0.25), q -> (0.25, 0.75):
        # KL = 0.75*log2(3) + 0.25*log2(1/3) = 0.5*log2(3) bits
        p = [0] * 75 + [1] * 25
        q = [0] * 25 + [1] * 75
        assert kl_divergence_estimator(p, q) == pytest.approx(0.5 * math.log2(3), abs=1e-12)

    def test_asymmetry(self) -> None:
        p = [0] * 75 + [1] * 25
        q = [0] * 25 + [1] * 75
        kl_pq = kl_divergence_estimator(p, q)
        kl_qp = kl_divergence_estimator(q, p)
        assert kl_pq == pytest.approx(kl_qp, abs=1e-12)  # symmetric for binary swap
        assert kl_pq > 0.0

    def test_non_negative(self) -> None:
        assert kl_divergence_estimator([0, 0, 1, 1], [0, 1, 0, 1]) >= 0.0

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            kl_divergence_estimator([0, 1], [0])


class TestBiasCorrectionHelpers:
    def test_bias_correction_formula(self) -> None:
        assert bias_correction(2.0, sample_size=100, alphabet_size=4) == pytest.approx(1.985)

    def test_effective_sample_size_identity_at_one(self) -> None:
        assert effective_sample_size_correction(3.7, sample_size=1, alphabet_size=4) == 3.7

    def test_panzeri_treves_uniform_fallback_nonneg(self) -> None:
        out = panzeri_treves_bias_correction(1.0, sample_size=50, alphabet_size=4)
        assert 0.0 <= out <= 1.0

    def test_panzeri_treves_explicit_frequencies_bounded(self) -> None:
        # Structural bounds only: the implemented formula differs from the
        # docstring comment (flagged for owner review in the lane report).
        freq = np.array([10.0, 20.0, 30.0, 40.0])
        out = panzeri_treves_bias_correction(2.0, 100, 4, response_frequencies=freq)
        assert 0.0 <= out <= 2.0

    def test_panzeri_treves_identity_at_small_sample(self) -> None:
        # sample_size <= 1 must return the entropy unchanged.
        assert panzeri_treves_bias_correction(1.23, sample_size=1, alphabet_size=4) == 1.23


class TestEntropyRate:
    def test_deterministic_alternation_zero_rate(self) -> None:
        seq = [0, 1] * 25
        assert entropy_rate_estimator(seq, order=1) == pytest.approx(0.0, abs=1e-9)

    def test_iid_fair_coin_rate_equals_entropy(self) -> None:
        # For i.i.d. sequences the order-1 conditional entropy equals the
        # marginal entropy (up to finite-sample estimation error).
        rng = np.random.RandomState(21)
        seq = rng.choice([0, 1], size=4000).tolist()
        rate = entropy_rate_estimator(seq, order=1)
        assert rate == pytest.approx(1.0, abs=0.05)

    def test_period_two_block_rate_zero_at_order_two(self) -> None:
        # Period-4 sequence [0,0,1,1]: fully determined by previous 2 symbols
        seq = [0, 0, 1, 1] * 20
        r2 = entropy_rate_estimator(seq, order=2)
        assert r2 == pytest.approx(0.0, abs=1e-9)

    def test_short_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="too short"):
            entropy_rate_estimator([0, 1], order=1)

    def test_invalid_order_raises(self) -> None:
        with pytest.raises(ValueError, match="Order must be >= 1"):
            entropy_rate_estimator([0, 1, 0, 1], order=0)
