"""AutoML cross-validation scores models for real, and a model that cannot be scored is reported as such.

``_cross_validate`` had two routes. For scikit-learn estimators it called ``cross_val_score(scoring=metric)``, which
rejects the documented metric name ``"mse"``, and any exception (that one, or a genuine fit failure) was swallowed
with ``except Exception: pass``. The "pure Python fallback" it then ran appended to a list called ``scores`` that only
existed on the scikit-learn route, so every fold raised ``NameError``, was caught by ``except Exception`` and scored
0.0: every model that is not a scikit-learn estimator, and every ``metric="mse"`` search, scored exactly 0.0.
A failing fit also scored 0.0, which is the best possible negated MSE.

``model_selection`` ranked a model that failed with ``score = 0.0`` next to the working ones (above every model with a
negative R^2) and let a NaN score through as "success". Scores are now computed or the error is raised, and failed
models carry ``score = None`` and are left out of the ranking.

Oracles: the fold arithmetic is redone by hand in the test, and ``cross_val_score`` with ``neg_mean_squared_error`` for
the scikit-learn route. No test doubles: the plain classes below are real ``fit`` / ``predict`` objects.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.ml.automl import optimization
from metainformant.ml.automl.optimization import _cross_validate, model_selection, random_search

sklearn = pytest.importorskip("sklearn")
from sklearn.linear_model import LinearRegression  # noqa: E402
from sklearn.model_selection import cross_val_score  # noqa: E402


class MajorityClassifier:
    """A plain ``fit`` / ``predict`` object, not a scikit-learn estimator."""

    def fit(self, X, y):
        values = list(np.asarray(y).tolist())
        self.label = max(set(values), key=values.count)
        return self

    def predict(self, X):
        return [self.label] * len(X)


class MeanRegressor:
    def fit(self, X, y):
        self.mean = float(np.mean(np.asarray(y, dtype=float)))
        return self

    def predict(self, X):
        return [self.mean] * len(X)


class Exploding:
    def fit(self, X, y):
        raise RuntimeError("cannot fit")

    def predict(self, X):
        return []


def _contiguous_folds(n: int, cv: int) -> list[tuple[list[int], list[int]]]:
    """(train indices, test indices) of the contiguous folds the pure-Python route uses."""
    size = n // cv
    folds = []
    for fold in range(cv):
        start = fold * size
        end = start + size if fold < cv - 1 else n
        test = list(range(start, end))
        folds.append(([i for i in range(n) if i not in set(test)], test))
    return folds


X = [[float(i)] for i in range(20)]
CLASSES = [0] * 14 + [1] * 6
TARGET = [float(i % 5) for i in range(20)]


class TestPurePythonRoute:
    def test_accuracy_is_computed_fold_by_fold(self) -> None:
        expected = []
        for train, test in _contiguous_folds(20, 5):
            labels = [CLASSES[i] for i in train]
            majority = max(set(labels), key=labels.count)
            expected.append(sum(CLASSES[i] == majority for i in test) / len(test))

        assert expected == [1.0, 1.0, 1.0, 0.5, 0.0]
        assert _cross_validate(MajorityClassifier(), X, CLASSES, cv=5, metric="accuracy") == pytest.approx(
            sum(expected) / 5
        )  # 0.7; the old route returned 0.0

    def test_mse_is_negated_and_r2_is_one_minus_the_residual_share(self) -> None:
        mse_folds, r2_folds = [], []
        for train, test in _contiguous_folds(20, 5):
            mean = sum(TARGET[i] for i in train) / len(train)
            truth = [TARGET[i] for i in test]
            residual = sum((t - mean) ** 2 for t in truth)
            total = sum((t - sum(truth) / len(truth)) ** 2 for t in truth)
            mse_folds.append(-residual / len(truth))
            r2_folds.append(1.0 - residual / total if total > 1e-15 else 0.0)

        assert _cross_validate(MeanRegressor(), X, TARGET, cv=5, metric="mse") == pytest.approx(sum(mse_folds) / 5)
        assert _cross_validate(MeanRegressor(), X, TARGET, cv=5, metric="r2") == pytest.approx(sum(r2_folds) / 5)
        assert sum(mse_folds) / 5 < -1.0  # a real, clearly non-zero score

    def test_a_model_that_cannot_be_fitted_raises_instead_of_scoring_zero(self) -> None:
        with pytest.raises(RuntimeError, match="cannot fit"):
            _cross_validate(Exploding(), X, CLASSES, cv=5, metric="accuracy")

    @pytest.mark.parametrize("estimator", [MajorityClassifier(), LinearRegression()])
    def test_an_unknown_metric_raises_for_both_routes(self, estimator) -> None:
        with pytest.raises(ValueError, match="Unsupported metric"):
            _cross_validate(estimator, X, TARGET, cv=3, metric="f-score")


class TestScikitLearnRoute:
    def test_mse_means_negated_mean_squared_error(self) -> None:
        rng = np.random.default_rng(0)
        features = rng.normal(size=(30, 2))
        target = (features @ np.array([1.5, -2.0]) + rng.normal(scale=0.3, size=30)).tolist()
        rows = features.tolist()

        got = _cross_validate(LinearRegression(), rows, target, cv=3, metric="mse")

        oracle = cross_val_score(LinearRegression(), features, np.array(target), cv=3, scoring="neg_mean_squared_error")
        assert got == pytest.approx(float(oracle.mean()), rel=1e-12)
        assert -1.0 < got < 0.0  # the old code answered 0.0 (a perfect fit) for the documented name "mse"

    def test_accuracy_and_r2_are_passed_through(self) -> None:
        rng = np.random.default_rng(1)
        features = rng.normal(size=(30, 2))
        target = (features @ np.array([1.0, 1.0])).tolist()

        r2 = _cross_validate(LinearRegression(), features.tolist(), target, cv=3, metric="r2")

        assert r2 == pytest.approx(float(cross_val_score(LinearRegression(), features, np.array(target), cv=3).mean()))
        assert r2 > 0.99

    def test_random_search_scores_real_mse_values(self) -> None:
        from sklearn.linear_model import Ridge

        rng = np.random.default_rng(2)
        features = rng.normal(size=(40, 3))
        target = (features @ np.array([2.0, 0.0, -1.0]) + rng.normal(scale=0.2, size=40)).tolist()

        result = random_search(
            lambda alpha: Ridge(alpha=alpha),
            {"alpha": [0.01, 1.0, 100.0, 10000.0]},
            features.tolist(),
            target,
            n_iter=8,
            cv=3,
            metric="mse",
            random_state=3,
        )

        scores = [r["score"] for r in result["all_results"]]
        assert all(math.isfinite(s) and s < 0.0 for s in scores)  # the old code scored every iteration 0.0
        by_alpha = {r["params"]["alpha"]: r["score"] for r in result["all_results"]}
        for alpha, score in by_alpha.items():
            oracle = cross_val_score(
                Ridge(alpha=alpha), features, np.array(target), cv=3, scoring="neg_mean_squared_error"
            ).mean()
            assert score == pytest.approx(float(oracle), abs=1e-6)
        assert result["best_params"]["alpha"] == max(by_alpha, key=by_alpha.get)


class TestModelSelection:
    def test_a_model_that_cannot_be_scored_is_reported_and_not_ranked(self) -> None:
        # four samples in two folds: the training folds have two samples, fewer than KNN's 3 neighbours
        features = [[0.0], [1.0], [2.0], [3.0]]
        target = [3.0, -1.0, 2.5, -2.0]

        result = model_selection(features, target, task="regression", cv=2)

        knn = result["cv_results_per_model"]["knn"]
        assert knn["score"] is None and knn["status"].startswith("failed")
        ranked = [name for name, _ in result["rankings"]]
        assert "knn" not in ranked
        assert ranked and all(math.isfinite(score) for _, score in result["rankings"])
        assert result["best_model_type"] == ranked[0]
        # the working models are scored, and the best of them is not the failed one
        assert all(
            info["status"] == "success" and math.isfinite(info["score"])
            for name, info in result["cv_results_per_model"].items()
            if name != "knn"
        )

    def test_when_every_model_fails_there_is_no_best_model(self) -> None:
        result = model_selection([[0.0], [1.0]], [0.0, 1.0], task="regression", cv=2)

        # R^2 is undefined on a single test sample, so no model can be scored
        assert result["rankings"] == []
        assert result["best_model_type"] == "none"
        assert all(
            info["score"] is None and info["status"].startswith("failed")
            for info in result["cv_results_per_model"].values()
        )

    def test_working_models_keep_their_scores(self) -> None:
        rng = np.random.default_rng(5)
        features = rng.normal(size=(60, 3))
        labels = (features[:, 0] > 0).astype(int).tolist()

        result = model_selection(features.tolist(), labels, task="classification", cv=3)

        assert result["best_model_type"] in {name for name, _ in result["rankings"]}
        assert all(0.0 <= score <= 1.0 for _, score in result["rankings"])
        assert max(score for _, score in result["rankings"]) > 0.8
        assert optimization.HAS_SKLEARN


def test_pure_python_singleton_r2_folds_are_rejected() -> None:
    with pytest.raises(ValueError, match="undefined"):
        _cross_validate(MeanRegressor(), [[float(i)] for i in range(4)], [0.0, 1.0, 2.0, 3.0], cv=4, metric="r2")


@pytest.mark.parametrize("prediction, expected", [([2.0, 2.0], 1.0), ([1.0, 1.0], 0.0)])
def test_pure_python_constant_r2_uses_finite_convention(prediction, expected) -> None:
    assert optimization._score_predictions("r2", [2.0, 2.0], prediction) == expected


def test_pure_python_r2_retains_tiny_variance() -> None:
    assert optimization._score_predictions("r2", [1e-9, 2e-9], [1e-9, 2e-9]) == 1.0
