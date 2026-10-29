"""Boruta's binomial test is exact, and recursive elimination does not score a failed fold as 0.0.

* ``_binomial_test_pvalue`` is documented as the one-sided binomial test and computed the normal approximation without
  continuity correction. That halves the tail for few trials (8 hits in 10 iterations: 0.029 for the exact 0.0547), so
  Boruta confirmed features at ``alpha = 0.05`` that the exact test does not confirm.
* ``_cross_validate_simple`` (the scorer of ``recursive_elimination``) caught every exception in a fold and scored the
  fold 0.0, so a model that could not be fitted looked like an inaccurate one, and returned 0.0 when ``cv`` left no
  fold at all.

Oracles: ``scipy.stats.binom.sf`` and fold arithmetic done by hand. No test doubles.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

stats = pytest.importorskip("scipy.stats")

from metainformant.ml.interpretability.feature_selection import (  # noqa: E402
    _binomial_test_pvalue,
    _cross_validate_simple,
    recursive_elimination,
)


@pytest.mark.parametrize("n", [1, 2, 5, 10, 20, 50, 100])
def test_the_binomial_tail_is_exact_for_the_boruta_null(n: int) -> None:
    for k in range(0, n + 2):
        assert _binomial_test_pvalue(k, n, 0.5) == pytest.approx(
            float(stats.binom.sf(k - 1, n, 0.5)), rel=1e-10, abs=1e-300
        )


@pytest.mark.parametrize("p", [0.05, 0.3, 0.9])
def test_the_binomial_tail_is_exact_for_other_probabilities(p: float) -> None:
    for n in (7, 40, 150):
        for k in (0, 1, n // 4, n // 2, n - 1, n):
            assert _binomial_test_pvalue(k, n, p) == pytest.approx(
                float(stats.binom.sf(k - 1, n, p)), rel=1e-10, abs=1e-300
            )


def test_eight_hits_in_ten_iterations_is_not_significant_at_five_percent() -> None:
    # exact: P(X >= 8) = 56 / 1024 = 0.0547 > 0.05; the normal approximation said 0.029 and Boruta confirmed the feature
    assert _binomial_test_pvalue(8, 10, 0.5) == pytest.approx(56 / 1024, rel=1e-12, abs=1e-300)
    assert _binomial_test_pvalue(8, 10, 0.5) > 0.05
    assert _binomial_test_pvalue(9, 10, 0.5) == pytest.approx(11 / 1024, rel=1e-12, abs=1e-300)
    assert _binomial_test_pvalue(9, 10, 0.5) < 0.05


def test_degenerate_arguments() -> None:
    assert _binomial_test_pvalue(0, 0, 0.5) == 1.0  # no trials, nothing observed
    assert _binomial_test_pvalue(3, 0, 0.5) == 0.0  # more successes than trials cannot happen
    assert _binomial_test_pvalue(0, 12, 0.5) == 1.0
    assert _binomial_test_pvalue(1, 12, 0.0) == 0.0 and _binomial_test_pvalue(12, 12, 1.0) == 1.0


class MajorityClassifier:
    def fit(self, X, y):
        values = list(np.asarray(y).tolist())
        self.label = max(set(values), key=values.count)
        return self

    def predict(self, X):
        return [self.label] * len(X)


class Exploding:
    def fit(self, X, y):
        raise RuntimeError("cannot fit")

    def predict(self, X):
        return []


def test_a_working_model_gets_its_fold_accuracies_averaged() -> None:
    X = [[float(i)] for i in range(20)]
    y = [0] * 14 + [1] * 6

    # contiguous folds of four: the majority label of the training part is 0 in every fold
    expected = [1.0, 1.0, 1.0, 0.5, 0.0]
    assert _cross_validate_simple(MajorityClassifier(), X, y, cv=5) == pytest.approx(sum(expected) / 5)


def test_a_model_that_cannot_be_fitted_raises_instead_of_scoring_zero() -> None:
    X = [[float(i)] for i in range(20)]
    y = [0] * 20

    with pytest.raises(RuntimeError, match="cannot fit"):
        _cross_validate_simple(Exploding(), X, y, cv=5)


def test_a_cv_that_leaves_no_fold_raises() -> None:
    with pytest.raises(ValueError, match="no fold"):
        _cross_validate_simple(MajorityClassifier(), [[0.0], [1.0]], [0, 1], cv=5)  # fold size 0


def test_recursive_elimination_reports_a_model_that_cannot_be_fitted() -> None:
    rng = np.random.default_rng(0)
    X = rng.normal(size=(30, 4)).tolist()
    y = [int(row[0] > 0) for row in X]

    with pytest.raises(RuntimeError, match="cannot fit"):
        recursive_elimination(Exploding(), X, y, n_features=2, cv=3)


def test_recursive_elimination_still_works_with_a_real_model() -> None:
    sklearn_linear = pytest.importorskip("sklearn.linear_model")
    rng = np.random.default_rng(1)
    X = rng.normal(size=(60, 4))
    y = (X[:, 0] + 0.1 * rng.normal(size=60) > 0).astype(int)
    order = rng.permutation(60)  # contiguous folds need shuffled rows

    result = recursive_elimination(
        sklearn_linear.LogisticRegression(max_iter=500), X[order].tolist(), y[order].tolist(), n_features=2, cv=3
    )

    final_score = result["cv_scores"][-1]["cv_score"]
    assert 0 in result["selected_features"]
    assert math.isfinite(final_score) and final_score > 0.7
