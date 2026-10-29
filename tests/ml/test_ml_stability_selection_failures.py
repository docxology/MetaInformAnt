"""``select_features_stability`` validates its method and counts only the bootstrap iterations that ran.

* An unknown ``method`` raised inside the bootstrap loop, was caught, logged and skipped every time; the all-zero
  selection frequencies then made ``np.argsort(...)[::-1]`` hand back the *last* features as "top-ranked" with only a
  warning, so ``method="bogus"`` returned an arbitrary subset. It now raises ``ValueError`` before any work.
* A bootstrap iteration that failed (for example an L1 logistic regression on a resample that holds a single class)
  still counted in the denominator of the selection frequency, as a vote against every feature. The frequency now
  uses the iterations that completed, and when none did the call raises ``RuntimeError``.

The oracle repeats the documented procedure with scikit-learn directly (same bootstrap draws from the same generator).
No test doubles.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

sklearn = pytest.importorskip("sklearn")
from sklearn.feature_selection import SelectFromModel  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402

from metainformant.ml.features.features import select_features_stability  # noqa: E402

pytestmark = pytest.mark.filterwarnings("ignore")


def _rare_class_data() -> tuple[np.ndarray, np.ndarray]:
    """40 samples, 3 positives: a few bootstrap resamples contain no positive at all."""
    rng = np.random.default_rng(8)
    X = rng.normal(size=(40, 6))
    y = np.zeros(40, dtype=int)
    y[rng.choice(40, size=3, replace=False)] = 1
    X[y == 1, 0] += 2.5
    return X, y


def test_an_unknown_method_raises_instead_of_returning_arbitrary_features() -> None:
    X, y = _rare_class_data()

    with pytest.raises(ValueError, match="Unknown stability method"):
        select_features_stability(X, y, method="bogus", n_bootstraps=5, random_state=1)


def test_a_non_positive_number_of_bootstraps_raises() -> None:
    X, y = _rare_class_data()

    with pytest.raises(ValueError, match="n_bootstraps"):
        select_features_stability(X, y, n_bootstraps=0, random_state=1)


def test_failed_iterations_are_not_votes_against_every_feature() -> None:
    X, y = _rare_class_data()
    n_bootstraps, threshold, seed = 30, 0.5, 7

    generator = np.random.default_rng(seed)
    counts, completed = np.zeros(X.shape[1]), 0
    for _ in range(n_bootstraps):
        rows = generator.choice(len(y), size=len(y), replace=True)
        selector = SelectFromModel(LogisticRegression(penalty="l1", solver="liblinear", random_state=seed))
        try:
            selector.fit(X[rows], y[rows])
        except ValueError:
            continue  # a resample with a single class cannot be fitted
        counts += selector.get_support().astype(int)
        completed += 1
    assert 0 < completed < n_bootstraps  # the data do make some iterations fail
    expected = np.where(counts / completed >= threshold)[0].tolist()
    assert expected != np.where(counts / n_bootstraps >= threshold)[0].tolist()  # and it matters

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, selected = select_features_stability(
            X, y, method="lasso", n_bootstraps=n_bootstraps, threshold=threshold, random_state=seed
        )

    assert selected == expected


def test_when_every_iteration_fails_the_call_raises() -> None:
    X = np.random.default_rng(0).normal(size=(30, 4))
    y = np.zeros(30, dtype=int)  # one class only: the L1 logistic regression cannot be fitted on any resample

    with pytest.raises(RuntimeError, match="bootstrap iterations failed"):
        select_features_stability(X, y, method="lasso", n_bootstraps=6, random_state=1)


def test_a_normal_run_is_unchanged() -> None:
    rng = np.random.default_rng(3)
    X = rng.normal(size=(120, 8))
    y = (X[:, 2] + 0.2 * rng.normal(size=120) > 0).astype(int)

    _, selected = select_features_stability(X, y, method="rf", n_bootstraps=10, threshold=0.7, random_state=1)

    assert 2 in selected
