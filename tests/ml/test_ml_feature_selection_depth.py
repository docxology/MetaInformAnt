"""Depth tests for metainformant.ml interpretability feature selection (Round-4 T3).

Value-pinned, zero-mock tests: real numeric optimization and information
computations on synthetic data with planted structure. Regression tests for
the recursive_elimination model-contract fix (fail-fast on non-estimators).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from sklearn.linear_model import LinearRegression, LogisticRegression

from metainformant.ml.interpretability.feature_selection import (
    _binomial_test_pvalue,
    _fit_lasso_pure,
    _mi_continuous,
    _mi_discrete_continuous,
    mutual_information_selection,
    recursive_elimination,
    stability_selection,
)


class TestLassoCoordinateDescent:
    """The pure-Python lasso must recover planted sparse structure."""

    def test_recovers_planted_sparse_coefficients(self) -> None:
        rng = np.random.RandomState(0)
        n, p = 200, 6
        X = rng.randn(n, p)
        beta = np.array([2.0, 0.0, -1.5, 0.0, 0.0, 0.0])
        y = X @ beta + rng.randn(n) * 0.1

        coefs = _fit_lasso_pure(X.tolist(), y.tolist(), alpha=0.01, max_iter=2000, tol=1e-8)
        assert len(coefs) == p
        # Signal coefficients recovered approximately; noise coefficients ~0
        assert abs(coefs[0] - 2.0) < 0.15
        assert abs(coefs[2] + 1.5) < 0.15
        for noise_idx in (1, 3, 4, 5):
            assert abs(coefs[noise_idx]) < 0.15

    def test_strong_alpha_zeroes_all_coefficients(self) -> None:
        rng = np.random.RandomState(1)
        X = rng.randn(100, 3)
        y = X @ np.array([1.0, 1.0, 1.0])
        coefs = _fit_lasso_pure(X.tolist(), y.tolist(), alpha=1e6, max_iter=100)
        assert all(c == 0.0 for c in coefs)

    def test_centering_invariant_to_constant_shift(self) -> None:
        # Lasso is fit on centered data: shifting y by a constant must not
        # change coefficients.
        rng = np.random.RandomState(2)
        X = rng.randn(150, 4)
        y = X @ np.array([1.0, -0.5, 0.0, 0.7]) + 3.0 + rng.randn(150) * 0.05
        c1 = _fit_lasso_pure(X.tolist(), y.tolist(), alpha=0.05)
        c2 = _fit_lasso_pure(X.tolist(), (y + 100.0).tolist(), alpha=0.05)
        assert all(abs(a - b) < 1e-6 for a, b in zip(c1, c2))

    def test_empty_input_returns_empty(self) -> None:
        assert _fit_lasso_pure([], [], alpha=0.1) == []
        assert _fit_lasso_pure([[]], [], alpha=0.1) == []


class TestMutualInformation:
    """MI estimators must rank informative features above noise features."""

    def test_informative_features_rank_first(self) -> None:
        rng = np.random.RandomState(3)
        n = 400
        informative = rng.randn(n, 3)
        noise = rng.randn(n, 10)
        y = (informative.sum(axis=1) > 0).astype(float)
        X = np.hstack([informative, noise])

        result = mutual_information_selection(X, y, n_features=3, discrete_target=True)
        top3 = set(result["selected_features"])
        assert top3 == {0, 1, 2}

    def test_mi_discrete_continuous_perfect_separation(self) -> None:
        # Two classes with disjoint value ranges => MI == H(class) in nats.
        continuous = [0.0] * 50 + [10.0] * 50
        discrete = [0.0] * 50 + [1.0] * 50
        mi = _mi_discrete_continuous(discrete, continuous, n_bins=10)
        h_class = math.log(2)  # balanced classes
        assert abs(mi - h_class) < 0.01

    def test_mi_independent_variables_near_zero(self) -> None:
        rng = np.random.RandomState(4)
        x = rng.randn(5000).tolist()
        y = rng.randn(5000).tolist()
        mi = _mi_continuous(x, y, n_bins=8)
        assert mi < 0.02

    def test_mi_non_negative(self) -> None:
        rng = np.random.RandomState(5)
        x = rng.randn(200).tolist()
        y = (np.array(x) > 0).astype(float).tolist()
        assert _mi_continuous(x, y) >= 0.0


class TestBinomialTestPvalue:
    """Upper-tail binomial p-value approximations pinned at known points."""

    def test_null_hypothesis_gives_large_pvalue(self) -> None:
        # k at the null mean should not reject.
        p = _binomial_test_pvalue(k=50, n=100, p=0.5)
        assert 0.4 < p <= 1.0

    def test_extreme_success_gives_small_pvalue(self) -> None:
        p = _binomial_test_pvalue(k=100, n=100, p=0.5)
        assert p < 1e-10

    def test_zero_trials_returns_one(self) -> None:
        assert _binomial_test_pvalue(k=0, n=0, p=0.5) == 1.0

    def test_degenerate_variance_zero_success(self) -> None:
        # variance ~ 0 and k == mean => not significant
        assert _binomial_test_pvalue(k=0, n=10, p=0.0) == 1.0

    def test_degenerate_variance_success(self) -> None:
        assert _binomial_test_pvalue(k=1, n=10, p=0.0) == 0.0

    def test_monotone_in_k(self) -> None:
        ps = [_binomial_test_pvalue(k, 100, 0.3) for k in range(30, 101, 10)]
        assert all(b <= a for a, b in zip(ps, ps[1:]))


class TestRecursiveEliminationContract:
    """RFE contract, including the fail-fast model-contract regression."""

    def test_rejects_non_estimator_with_typeerror(self) -> None:
        # Regression: a string model previously degraded silently into
        # correlation-based importances with only a log warning.
        with pytest.raises(TypeError, match="fit.*predict|estimator"):
            recursive_elimination(model="ridge", X=[[1.0, 2.0], [3.0, 4.0]], y=[0, 1], n_features=1)

    def test_rejects_object_without_predict(self) -> None:
        class FitOnly:
            def fit(self, X, y):
                return self

        with pytest.raises(TypeError, match="estimator"):
            recursive_elimination(model=FitOnly(), X=[[1.0], [2.0]], y=[0, 1], n_features=1)

    def test_selects_signal_feature_with_real_estimator(self) -> None:
        rng = np.random.RandomState(6)
        n, p = 80, 8
        X = rng.randn(n, p)
        y = X[:, 0] * 2 + rng.randn(n) * 0.1  # only feature 0 matters

        result = recursive_elimination(model=LinearRegression(), X=X, y=y, n_features=3, step=1, cv=3)
        assert len(result["selected_features"]) == 3
        assert 0 in result["selected_features"]
        assert len(result["ranking"]) == p
        # Documented contract: selected features carry rank 1
        selected_ranks = sorted(result["ranking"][i] for i in result["selected_features"])
        assert selected_ranks == [1, 1, 1]
        # cv_scores recorded for every step plus the final fit
        assert all("n_features" in s and "cv_score" in s for s in result["cv_scores"])

    def test_classification_path_selects_signal(self) -> None:
        rng = np.random.RandomState(9)
        n, p = 120, 6
        X = rng.randn(n, p)
        y = (X[:, 2] > 0).astype(int)

        result = recursive_elimination(model=LogisticRegression(max_iter=1000), X=X, y=y, n_features=2, step=1, cv=3)
        assert 2 in result["selected_features"]

    def test_n_features_exceeding_total_raises(self) -> None:
        with pytest.raises(ValueError, match="must be <="):
            recursive_elimination(model=LinearRegression(), X=[[1, 2]], y=[0], n_features=5)


class TestStabilitySelectionContract:
    """Stability selection on strongly separable data keeps signal features."""

    def test_signal_features_most_frequently_selected(self) -> None:
        rng = np.random.RandomState(7)
        n = 200
        X_signal = rng.randn(n, 2)
        X_noise = rng.randn(n, 6)
        X = np.hstack([X_signal, X_noise])
        y = (X_signal[:, 0] + X_signal[:, 1] > 0).astype(int)

        result = stability_selection(X, y, n_bootstrap=20, random_state=42)
        probs = result["selection_probabilities"]
        assert len(probs) == 8
        # Each signal feature selected at least as often as any noise feature
        assert min(probs[:2]) >= max(probs[2:])

    def test_reproducible_with_random_state(self) -> None:
        rng = np.random.RandomState(8)
        X = rng.randn(120, 5)
        y = (X[:, 0] > 0).astype(int)
        r1 = stability_selection(X, y, n_bootstrap=10, random_state=11)
        r2 = stability_selection(X, y, n_bootstrap=10, random_state=11)
        assert r1["selection_probabilities"] == r2["selection_probabilities"]
        assert r1["selected_features"] == r2["selected_features"]
