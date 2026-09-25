"""Regression tests for the Chao-Shen entropy estimator.

The current implementation is the published Chao-Shen (2003) estimator --
identical to ``entropy.ChaoShen`` in the R `entropy` package (Strimmer):
coverage C = 1 - f1/n, adjusted probabilities pa = C*p, Horvitz-Thompson
inclusion probabilities la = 1 - (1-pa)^n, H = -sum pa*log2(pa)/la.
These tests pin it against hand-computed values and structural invariants.
(An earlier implementation used an ad-hoc singleton/non-singleton surrogate;
its pinned values were replaced by the published-estimator values below.)
real implementations: all values computed from real arithmetic.
"""

from __future__ import annotations


import numpy as np

from metainformant.information.metrics.core.estimation import (
    _chao_shen_entropy_estimator,
    entropy_bootstrap_confidence,
    entropy_estimator,
)


def _reference_chao_shen(counts: list[int]) -> float:
    """Independent hand-rolled reference implementation for cross-checking."""
    arr = np.array(counts, dtype=float)
    arr = arr[arr > 0]
    total = int(arr.sum())
    if total == 0:
        return 0.0
    f1 = int(np.sum(arr == 1))
    coverage = 1.0 - f1 / total
    if coverage <= 0.0:
        return 0.0
    adjusted = coverage * arr / total
    inclusion = 1.0 - (1.0 - adjusted) ** total
    return max(0.0, float(-np.sum(adjusted * np.log2(adjusted) / inclusion)))


class TestChaoShenEstimator:
    """Value-pinned regression tests for the coverage-adjusted estimator."""

    def test_matches_hand_computed_value_sparse(self) -> None:
        # counts {A:50, T:2, G:1, C:1}: n=54, f1=2, C=52/54; value of the
        # published estimator (was 0.5243226408289577 under the removed
        # ad-hoc surrogate).
        counts = np.array([50, 2, 1, 1])
        expected = 0.6805093130489859
        got = _chao_shen_entropy_estimator(counts, total=54)
        assert abs(got - expected) < 1e-12

    def test_matches_independent_reference_implementation(self) -> None:
        for counts in ([50, 2, 1, 1], [10, 10, 10, 1], [5, 3, 2, 1, 1, 1], [100, 1]):
            got = _chao_shen_entropy_estimator(np.array(counts), total=sum(counts))
            ref = _reference_chao_shen(counts)
            assert abs(got - ref) < 1e-12

    def test_no_singletons_equals_uncorrected_plugin(self) -> None:
        # With no singletons, coverage C = 1 and adjusted probs are the
        # empirical probs, so Chao-Shen must equal the *uncorrected* plugin
        # (Miller-Madow correction is a separate method concern).
        counts = {"A": 50, "T": 30}
        cs = entropy_estimator(counts, method="chao_shen")
        p = np.array([50, 30]) / 80.0
        plugin_uncorrected = -float(np.sum(p * np.log2(p)))
        assert abs(cs - plugin_uncorrected) < 1e-12

    def test_all_singletons_returns_zero(self) -> None:
        # Every category observed once => coverage estimate C = 0.
        assert _chao_shen_entropy_estimator(np.array([1, 1, 1]), total=3) == 0.0

    def test_exposed_via_public_dispatcher(self) -> None:
        counts = {"A": 50, "T": 2, "G": 1, "C": 1}
        h = entropy_estimator(counts, method="chao_shen")
        assert h >= 0.0
        assert abs(h - 0.6805093130489859) < 1e-12

    def test_zero_and_empty_inputs(self) -> None:
        assert _chao_shen_entropy_estimator(np.array([0, 0]), total=0) == 0.0
        assert _chao_shen_entropy_estimator(np.array([], dtype=int), total=0) == 0.0

    def test_bootstrap_confidence_reproducible(self) -> None:
        counts = {"A": 50, "T": 30, "G": 20}
        r1 = entropy_bootstrap_confidence(counts, n_bootstraps=50, random_state=7)
        r2 = entropy_bootstrap_confidence(counts, n_bootstraps=50, random_state=7)
        assert r1 == r2
        assert r1["ci_lower"] <= r1["entropy"] <= r1["ci_upper"]
