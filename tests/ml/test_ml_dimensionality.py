"""Tests for the consolidated dimensionality reduction family.

``metainformant.ml.features.dimensionality`` keeps a single standardization
path: ``pca_reduction`` and ``ica_reduction`` standardize through the shared
``_scale_features`` helper, and the ``reduce_dimensions_*`` functions are thin
delegating wrappers (``standardize`` is only a deprecated alias of
``scale_data``). Tests exercise real scikit-learn models on deterministic
data; no test doubles.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.ml.features.dimensionality import (
    compare_dimensionality_methods,
    ica_reduction,
    pca_reduction,
    reduce_dimensions_pca,
)

sklearn_decomposition = pytest.importorskip("sklearn.decomposition")
sklearn_preprocessing = pytest.importorskip("sklearn.preprocessing")


def _correlated_features(n_samples: int = 60, seed: int = 7) -> np.ndarray:
    """Deterministic 5-feature matrix with offset, correlated columns."""
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(n_samples, 2))
    X = np.empty((n_samples, 5))
    X[:, 0] = base[:, 0] * 10.0 + 100.0
    X[:, 1] = base[:, 0] * 9.0 + base[:, 1] + 50.0
    X[:, 2] = base[:, 1] * 3.0 + 7.0
    X[:, 3] = base[:, 0] + base[:, 1] - 5.0
    X[:, 4] = rng.normal(scale=0.1, size=n_samples)
    return X


class TestPcaReduction:
    def test_standardized_scores_match_single_standardization_path(self) -> None:
        """Scores must equal the projection of StandardScaler-scaled inputs."""
        X = _correlated_features()
        X_pca, model = pca_reduction(X, n_components=4, random_state=0)

        assert X_pca.shape == (60, 4)
        scaler = sklearn_preprocessing.StandardScaler()
        X_scaled = scaler.fit_transform(X)

        assert np.allclose(X_scaled @ model.components_.T, X_pca)

    def test_pca_scores_mean_and_variance(self) -> None:
        """Standardized PCA scores are zero-mean with variance == explained_variance_."""
        X = _correlated_features()
        X_pca, model = pca_reduction(X, n_components=4, random_state=0)

        assert np.allclose(X_pca.mean(axis=0), 0.0, atol=1e-10)
        assert np.allclose(X_pca.var(axis=0, ddof=1), model.explained_variance_)

        ratios = model.explained_variance_ratio_
        assert np.all(np.diff(ratios) <= 0), (
            "explained variance must be sorted descending"
        )
        assert 0.0 < ratios.sum() <= 1.0 + 1e-9

    def test_unscaled_pca_uses_raw_mean_and_variance(self) -> None:
        """scale_data=False must center (never rescale) the raw matrix."""
        X = _correlated_features()
        X_raw, model_raw = pca_reduction(
            X, n_components=3, scale_data=False, random_state=0
        )

        X_centered = X - X.mean(axis=0)
        assert np.allclose(X_centered @ model_raw.components_.T, X_raw)
        assert np.allclose(X_raw.var(axis=0, ddof=1), model_raw.explained_variance_)

        _, model_scaled = pca_reduction(
            X, n_components=3, scale_data=True, random_state=0
        )
        assert not np.allclose(
            model_raw.explained_variance_, model_scaled.explained_variance_
        )


class TestIcaReduction:
    def test_ica_returns_finite_components_with_seed(self) -> None:
        X = _correlated_features()
        X_ica, ica = ica_reduction(X, n_components=2, random_state=0, max_iter=1000)

        assert X_ica.shape == (60, 2)
        assert np.all(np.isfinite(X_ica))
        assert ica.mixing_.shape == (5, 2)

    def test_ica_is_reproducible_with_random_state(self) -> None:
        X = _correlated_features()
        X_a, _ = ica_reduction(X, n_components=2, random_state=0, max_iter=1000)
        X_b, _ = ica_reduction(X, n_components=2, random_state=0, max_iter=1000)

        assert np.allclose(X_a, X_b)


class TestReduceDimensionsWrappers:
    def test_wrapper_delegates_to_pca_reduction(self) -> None:
        """reduce_dimensions_pca repackages pca_reduction, standardization included."""
        X = _correlated_features()
        X_reduced, components, explained_var = reduce_dimensions_pca(
            X, n_components=3, random_state=0
        )

        X_pca, model = pca_reduction(X, n_components=3, random_state=0)

        assert np.allclose(X_reduced, X_pca)
        assert np.allclose(components, model.components_.T)
        assert np.allclose(explained_var, model.explained_variance_ratio_)

    def test_standardize_alias_matches_scale_data(self) -> None:
        X = _correlated_features()
        X_scaled, _, _ = reduce_dimensions_pca(
            X, n_components=3, scale_data=False, random_state=0
        )
        X_aliased, _, _ = reduce_dimensions_pca(
            X, n_components=3, standardize=False, random_state=0
        )

        assert np.allclose(X_scaled, X_aliased)

    def test_default_component_count_is_capped(self) -> None:
        X = _correlated_features()
        X_reduced, components, explained_var = reduce_dimensions_pca(X, random_state=0)

        assert X_reduced.shape == (60, 5)
        assert components.shape == (5, 5)
        assert explained_var.shape == (5,)


class TestCompareDimensionalityMethods:
    def test_compare_pca_and_ica(self) -> None:
        X = _correlated_features()

        report = compare_dimensionality_methods(
            X, methods=["pca", "ica"], n_components=2, random_state=0
        )

        assert report["input_shape"] == (60, 5)
        assert report["n_components"] == 2
        assert set(report["embeddings"]) == {"pca", "ica"}
        for method in ("pca", "ica"):
            assert report["comparison"][method]["success"] is True
            assert report["comparison"][method]["shape"] == (60, 2)
            assert report["embeddings"][method].shape == (60, 2)
        pca_variance = report["comparison"]["pca"]["explained_variance"]
        assert pca_variance is not None and 0.0 < pca_variance <= 1.0 + 1e-9
        assert report["comparison"]["ica"]["explained_variance"] is None

    def test_unknown_method_is_recorded_as_failure(self) -> None:
        X = _correlated_features()

        report = compare_dimensionality_methods(X, methods=["nope"], n_components=2)

        assert report["comparison"]["nope"]["success"] is False
        assert "Unknown method" in report["comparison"]["nope"]["error"]
        assert "nope" not in report["embeddings"]
