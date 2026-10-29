"""Unregularized CCA is refused where every canonical correlation would be 1, and both CCAs share one core.

Regression tests for ``metainformant.multiomics.analysis.integration.canonical_correlation`` and
``metainformant.multiomics.methods.factorization.canonical_correlation``.

Defects: with ``regularization=0`` and a layer of ``n_samples - 1`` or more informative features, both functions
returned canonical correlations of exactly 1 whatever the data were (the samples, once centred, span only
``n - 1`` dimensions, so such a layer spans all of them and contains the other layer). The integration version only
logged a warning; the factorization version whitened with a ``1e-10`` jitter and returned near-1 numbers silently.
The two were also separate implementations (a covariance-matrix SVD with arbitrary signs against a thin-SVD dual
form with a fixed sign convention), and the integration version filled a missing value with the mean of the *whole*
layer, an unrelated feature's scale.

Now both run the same numerical core and both raise ``ValueError`` asking for ``regularization > 0``.

Oracles: a QR/SVD principal-angle computation of classical CCA, a dense ``p x p`` ridge CCA, and ``np.corrcoef``.
Real computation; no test doubles (real-implementation policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import linalg

from metainformant.multiomics.analysis.integration import canonical_correlation as integration_cca
from metainformant.multiomics.methods.factorization import canonical_correlation as factorization_cca

DEGENERATE = r"Unregularized CCA is degenerate.*regularization > 0"


def _layers(x: np.ndarray, y: np.ndarray) -> dict[str, pd.DataFrame]:
    index = [f"S{i}" for i in range(x.shape[0])]
    return {"rna": pd.DataFrame(x, index=index), "protein": pd.DataFrame(y, index=index)}


def _standardize(a: np.ndarray) -> np.ndarray:
    """Columns centred and scaled to unit sample variance, the scale the integration version works on."""
    return (a - a.mean(axis=0)) / a.std(axis=0, ddof=1)


def _random_layers(n: int, p: int, q: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    latent = rng.standard_normal((n, 2))
    x = latent @ rng.standard_normal((2, p)) + rng.standard_normal((n, p))
    y = latent @ rng.standard_normal((2, q)) + rng.standard_normal((n, q))
    return x, y


def _classical_correlations(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Classical CCA as the cosines of the principal angles between the column spaces of the centred layers.

    ``scipy.linalg.orth`` gives a rank-revealing orthonormal basis, so layers whose columns are linearly dependent
    are handled (a plain QR would keep directions that are not in the column space).
    """
    qx = linalg.orth(x - x.mean(axis=0))
    qy = linalg.orth(y - y.mean(axis=0))
    return np.linalg.svd(qx.T @ qy, compute_uv=False)


def _run_both(x: np.ndarray, y: np.ndarray, **kwargs: Any) -> None:
    integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=2, **kwargs)
    factorization_cca(x, y, n_components=2, **kwargs)


class TestUnregularizedCcaIsRefusedWhereItIsDegenerate:
    @pytest.mark.parametrize("n_features", [29, 30, 31, 60, 400])
    def test_a_layer_with_n_minus_one_or_more_features_raises_in_both_functions(self, n_features: int) -> None:
        # 30 samples: after centring they span 29 dimensions, so 29 informative features already span all of them
        x, y = _random_layers(30, n_features, 6)

        with pytest.raises(ValueError, match=DEGENERATE):
            integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=3, regularization=0.0)
        with pytest.raises(ValueError, match=DEGENERATE):
            factorization_cca(x, y, n_components=3, regularization=0.0)

    def test_it_does_not_matter_which_layer_is_the_wide_one(self) -> None:
        x, y = _random_layers(30, 6, 45)

        with pytest.raises(ValueError, match=DEGENERATE):
            integration_cca(_layers(x, y), layers=("rna", "protein"), regularization=0.0)
        with pytest.raises(ValueError, match=DEGENERATE):
            factorization_cca(x, y, regularization=0.0)

    def test_the_default_of_the_integration_function_is_unregularized_and_so_is_refused_for_wide_data(self) -> None:
        x, y = _random_layers(25, 80, 70)

        with pytest.raises(ValueError, match=DEGENERATE):
            integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=4)

    def test_the_message_explains_the_degeneracy_and_names_the_layer(self) -> None:
        x, y = _random_layers(20, 5, 30)

        with pytest.raises(ValueError) as from_integration:
            integration_cca(_layers(x, y), layers=("rna", "protein"), regularization=0.0)
        with pytest.raises(ValueError) as from_factorization:
            factorization_cca(x, y, regularization=0.0)

        message = str(from_integration.value)
        assert "layer 'protein' has 30 features of numerical rank 19" in message
        assert "20 centred samples span only 19 dimensions" in message
        assert "every canonical correlation is exactly 1" in message
        assert "regularization > 0" in message
        assert "layer 'Y' has 30 features of numerical rank 19" in str(from_factorization.value)

    @pytest.mark.parametrize("n_samples", [5, 12, 30])
    def test_the_boundary_is_n_minus_one_features(self, n_samples: int) -> None:
        narrow, other = _random_layers(n_samples, n_samples - 2, 1)  # rank n - 2: room for a one-feature partner
        edge, other_2 = _random_layers(n_samples, n_samples - 1, 1)  # rank n - 1: spans everything

        _run_both(narrow, other, regularization=0.0)  # no error
        with pytest.raises(ValueError, match=DEGENERATE):
            _run_both(edge, other_2, regularization=0.0)

    def test_two_samples_always_correlate_perfectly_so_they_need_a_ridge(self) -> None:
        x, y = _random_layers(2, 3, 2)

        with pytest.raises(ValueError, match=DEGENERATE):
            _run_both(x, y, regularization=0.0)
        _run_both(x, y, regularization=0.1)  # a ridge gives a well-defined (if uninformative) answer

    def test_the_refusal_is_justified_every_correlation_really_is_one(self) -> None:
        # the claim in the message, checked by classical CCA computed from orthonormal bases
        x, y = _random_layers(30, 40, 6)

        np.testing.assert_allclose(_classical_correlations(x, y), np.ones(6), atol=1e-10)

    def test_any_positive_ridge_is_accepted_for_the_same_data(self) -> None:
        x, y = _random_layers(30, 400, 60)

        for lam in (1e-8, 0.01, 0.5, 100.0):
            _, _, _, _, corr = integration_cca(
                _layers(x, y), layers=("rna", "protein"), n_components=3, regularization=lam
            )
            result = factorization_cca(x, y, n_components=3, regularization=lam)

            assert np.all((0.0 <= corr) & (corr <= 1.0))
            assert np.all((0.0 <= np.array(result["correlations"])) & (np.array(result["correlations"]) <= 1.0))

    def test_the_default_ridge_of_the_factorization_function_never_raises_for_wide_data(self) -> None:
        x, y = _random_layers(20, 200, 150)

        result = factorization_cca(x, y, n_components=3)  # regularization=0.1

        assert len(result["correlations"]) == 3 and max(result["correlations"]) < 1.0


class TestRankNotFeatureCountDecides:
    def test_many_features_of_low_rank_are_not_degenerate(self) -> None:
        # 60 features that are linear combinations of 3 signals: rank 3 on 30 samples, nowhere near 29
        rng = np.random.default_rng(3)
        signals = rng.standard_normal((30, 3))
        x = signals @ rng.standard_normal((3, 60))
        y = signals[:, :2] @ rng.standard_normal((2, 5)) + 0.3 * rng.standard_normal((30, 5))

        _, _, _, _, corr = integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=2, regularization=0.0)

        np.testing.assert_allclose(corr, _classical_correlations(x, y)[:2], atol=1e-9)
        assert corr[0] > 0.9 and corr[1] < 1.0 - 1e-6  # informative, not the trivial all-ones answer

    def test_constant_features_do_not_count_towards_the_rank(self) -> None:
        x, y = _random_layers(30, 29, 3)  # 29 informative features: rank 29 = n - 1
        x_padded = np.hstack([x, np.full((30, 20), 4.0)])

        with pytest.raises(ValueError, match=r"layer 'rna' has 49 features of numerical rank 29"):
            integration_cca(_layers(x_padded, y), layers=("rna", "protein"), regularization=0.0)

        x_fewer, _ = _random_layers(30, 20, 3)
        x_fewer_padded = np.hstack([x_fewer, np.full((30, 20), 4.0)])  # 40 features, rank 20: fine
        integration_cca(_layers(x_fewer_padded, y), layers=("rna", "protein"), n_components=2, regularization=0.0)

    def test_partly_degenerate_data_is_refused_and_ridge_remains_available(self) -> None:
        # ranks 20 + 15 = 35 > 29 force six correlations to one.
        x, y = _random_layers(30, 20, 15, seed=4)
        classical = _classical_correlations(x, y)[:12]
        np.testing.assert_allclose(classical[:6], np.ones(6), atol=1e-9)
        assert np.all(classical[6:] < 1.0 - 1e-6)

        with pytest.raises(ValueError, match=DEGENERATE):
            integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=12, regularization=0.0)
        with pytest.raises(ValueError, match=DEGENERATE):
            factorization_cca(x, y, n_components=12, regularization=0.0)

        _, _, _, _, corr = integration_cca(
            _layers(x, y), layers=("rna", "protein"), n_components=12, regularization=0.1
        )
        assert np.all((corr >= 0.0) & (corr < 1.0))


class TestBothFunctionsRunOneCore:
    SHAPES = [(120, 8, 6), (40, 60, 25), (25, 10, 200), (30, 29, 29)]

    @pytest.mark.parametrize(("n", "p", "q"), SHAPES)
    @pytest.mark.parametrize("lam", [0.01, 0.1, 1.0, 25.0])
    def test_same_correlations_variates_and_weights_for_the_same_standardized_data(
        self, n: int, p: int, q: int, lam: float
    ) -> None:
        x, y = _random_layers(n, p, q, seed=n + p)
        x_c, y_c, x_w, y_w, corr = integration_cca(
            _layers(x, y), layers=("rna", "protein"), n_components=4, regularization=lam
        )

        other = factorization_cca(_standardize(x), _standardize(y), n_components=4, regularization=lam)

        np.testing.assert_allclose(corr, other["correlations"], rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(x_c, np.array(other["x_scores"]), atol=1e-8)  # identical signs, no flipping
        np.testing.assert_allclose(y_c, np.array(other["y_scores"]), atol=1e-8)
        np.testing.assert_allclose(x_w, np.array(other["x_weights"]), atol=1e-8)
        np.testing.assert_allclose(y_w, np.array(other["y_weights"]), atol=1e-8)

    @pytest.mark.parametrize(("n", "p", "q"), [(120, 8, 6), (50, 12, 9)])
    def test_unregularized_results_agree_too_where_they_are_defined(self, n: int, p: int, q: int) -> None:
        x, y = _random_layers(n, p, q, seed=1)

        _, _, _, _, corr = integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=4, regularization=0.0)
        other = factorization_cca(_standardize(x), _standardize(y), n_components=4, regularization=0.0)

        np.testing.assert_allclose(corr, other["correlations"], rtol=1e-9)
        np.testing.assert_allclose(corr, _classical_correlations(x, y)[:4], atol=1e-9)

    @pytest.mark.parametrize(("n", "p", "q"), SHAPES[:3])
    @pytest.mark.parametrize("lam", [0.05, 2.0])
    def test_both_match_the_dense_textbook_definition(self, n: int, p: int, q: int, lam: float) -> None:
        x, y = _random_layers(n, p, q, seed=7)
        xs, ys = _standardize(x), _standardize(y)
        cxx = xs.T @ xs / (n - 1) + lam * np.eye(p)
        cyy = ys.T @ ys / (n - 1) + lam * np.eye(q)

        def inverse_sqrt(c: np.ndarray) -> np.ndarray:
            values, vectors = np.linalg.eigh(c)
            return vectors @ np.diag(values**-0.5) @ vectors.T

        expected = np.linalg.svd(inverse_sqrt(cxx) @ (xs.T @ ys / (n - 1)) @ inverse_sqrt(cyy), compute_uv=False)[:3]

        _, _, _, _, corr = integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=3, regularization=lam)
        other = factorization_cca(xs, ys, n_components=3, regularization=lam)

        np.testing.assert_allclose(corr, expected, rtol=1e-9)
        np.testing.assert_allclose(other["correlations"], expected, rtol=1e-9)

    def test_the_factorization_function_centres_but_does_not_standardize(self) -> None:
        x, y = _random_layers(60, 6, 5, seed=2)
        scaled = x * np.array([1.0, 10.0, 100.0, 0.1, 5.0, 2.0]) + 7.0  # per-feature scale and offset

        by_raw = factorization_cca(scaled, y, n_components=3, regularization=0.5)["correlations"]
        by_standardized = factorization_cca(_standardize(scaled), _standardize(y), n_components=3, regularization=0.5)[
            "correlations"
        ]
        centred_only = factorization_cca(x * 1.0, y, n_components=3, regularization=0.5)["correlations"]

        assert not np.allclose(by_raw, by_standardized, atol=1e-3)  # the ridge is on the covariance scale
        assert not np.allclose(centred_only, by_raw, atol=1e-3)
        # but a pure offset changes nothing
        np.testing.assert_allclose(
            factorization_cca(x + 123.0, y - 4.0, n_components=3, regularization=0.5)["correlations"],
            factorization_cca(x, y, n_components=3, regularization=0.5)["correlations"],
            rtol=1e-9,
        )

    def test_results_do_not_change_from_run_to_run_and_the_sign_convention_is_fixed(self) -> None:
        x, y = _random_layers(40, 15, 9, seed=9)

        first = factorization_cca(x, y, n_components=3, regularization=0.3)
        second = factorization_cca(x, y, n_components=3, regularization=0.3)

        assert first == second
        weights = np.array(first["x_weights"])
        assert np.all(weights[np.argmax(np.abs(weights), axis=0), np.arange(3)] > 0)  # largest weight positive

    def test_scores_are_the_centred_data_times_the_weights_and_loadings_are_correlations(self) -> None:
        x, y = _random_layers(35, 7, 5, seed=11)

        result = factorization_cca(x, y, n_components=2, regularization=0.4)

        x_scores, y_scores = np.array(result["x_scores"]), np.array(result["y_scores"])
        np.testing.assert_allclose(x_scores, (x - x.mean(axis=0)) @ np.array(result["x_weights"]), atol=1e-10)
        np.testing.assert_allclose(y_scores, (y - y.mean(axis=0)) @ np.array(result["y_weights"]), atol=1e-10)
        expected = np.array([[np.corrcoef(x[:, j], x_scores[:, c])[0, 1] for c in range(2)] for j in range(7)])
        np.testing.assert_allclose(np.array(result["x_loadings"]), expected, atol=1e-10)
        expected_y = np.array([[np.corrcoef(y[:, j], y_scores[:, c])[0, 1] for c in range(2)] for j in range(5)])
        np.testing.assert_allclose(np.array(result["y_loadings"]), expected_y, atol=1e-10)

    def test_a_constant_feature_has_an_undefined_loading_not_a_warning_or_a_made_up_number(self) -> None:
        x, y = _random_layers(30, 5, 4, seed=5)
        x[:, 2] = 3.0

        result = factorization_cca(x, y, n_components=2, regularization=0.2)  # a RuntimeWarning here would fail

        loadings = np.array(result["x_loadings"])
        assert np.all(np.isnan(loadings[2])) and not np.any(np.isnan(np.delete(loadings, 2, axis=0)))


class TestMissingValuesAndArguments:
    def test_the_integration_version_imputes_each_feature_with_its_own_mean(self) -> None:
        x, y = _random_layers(40, 6, 5, seed=13)
        x[3, 1] = np.nan
        x[10, 4] = np.nan
        y[7, 0] = np.nan
        layers = _layers(x, y)
        explicitly_imputed = {name: frame.fillna(frame.mean()) for name, frame in layers.items()}

        mine = integration_cca(layers, layers=("rna", "protein"), n_components=3, regularization=0.5)
        expected = integration_cca(explicitly_imputed, layers=("rna", "protein"), n_components=3, regularization=0.5)

        for got, want in zip(mine, expected):
            np.testing.assert_allclose(got, want, atol=1e-12)

    def test_imputation_does_not_use_the_mean_of_the_whole_layer(self) -> None:
        # a feature on a large scale next to one on a small scale: the old fill (mean of every entry) is far off
        rng = np.random.default_rng(14)
        x = np.column_stack([rng.normal(1000.0, 10.0, 40), rng.normal(0.0, 1.0, 40), rng.normal(5.0, 2.0, 40)])
        y = x[:, :2] @ rng.standard_normal((2, 3)) + rng.standard_normal((40, 3))
        x_missing = x.copy()
        x_missing[5, 1] = np.nan
        global_mean_fill = x.copy()
        global_mean_fill[5, 1] = np.nanmean(x_missing)  # what the previous implementation put in the gap

        kwargs: dict[str, Any] = {"layers": ("rna", "protein"), "n_components": 2, "regularization": 0.5}
        own_mean = integration_cca(_layers(x_missing, y), **kwargs)[4]
        wrong_fill = integration_cca(_layers(global_mean_fill, y), **kwargs)[4]
        column_fill = x.copy()
        column_fill[5, 1] = np.nanmean(x_missing[:, 1])
        right_fill = integration_cca(_layers(column_fill, y), **kwargs)[4]

        np.testing.assert_allclose(own_mean, right_fill, atol=1e-12)
        assert not np.allclose(own_mean, wrong_fill, atol=1e-6)

    def test_a_feature_with_no_observed_value_cannot_be_imputed(self) -> None:
        x, y = _random_layers(20, 4, 3)
        x[:, 2] = np.nan

        with pytest.raises(ValueError, match=r"layer 'rna' has 1 feature\(s\) with no observed value"):
            integration_cca(_layers(x, y), layers=("rna", "protein"), regularization=0.5)

    def test_the_factorization_version_refuses_missing_values_instead_of_guessing(self) -> None:
        x, y = _random_layers(20, 4, 3)
        x[1, 1] = np.nan

        with pytest.raises(ValueError, match=r"NaN or infinite"):
            factorization_cca(x, y, regularization=0.5)

    @pytest.mark.parametrize("bad", [-0.1, float("nan"), float("inf"), float("-inf"), "0.1", None])
    def test_the_factorization_version_validates_the_ridge(self, bad: Any) -> None:
        x, y = _random_layers(20, 4, 3)

        with pytest.raises(ValueError, match=r"regularization must be a finite number >= 0"):
            factorization_cca(x, y, regularization=bad)

    def test_numpy_scalar_ridges_are_accepted(self) -> None:
        x, y = _random_layers(20, 4, 3)

        for lam in (np.float32(0.5), np.float64(0.5), np.int64(1), 1):
            assert factorization_cca(x, y, n_components=1, regularization=lam)["correlations"]

    @pytest.mark.parametrize("n_samples", [0, 1])
    def test_fewer_than_two_samples_raise(self, n_samples: int) -> None:
        with pytest.raises(ValueError, match=r"at least 2 samples"):
            factorization_cca(np.ones((n_samples, 3)), np.ones((n_samples, 2)))

    def test_shape_and_component_errors_keep_their_messages(self) -> None:
        x, y = _random_layers(20, 4, 3)

        with pytest.raises(ValueError, match="n_components"):
            factorization_cca(x, y, n_components=0)
        with pytest.raises(ValueError, match="mismatch"):
            factorization_cca(x, y[:15])
        with pytest.raises(ValueError, match="2-D"):
            factorization_cca(x[:, 0], y)

    def test_more_components_than_the_ranks_support_are_clamped_for_both_functions(self) -> None:
        x, y = _random_layers(40, 3, 4)

        assert len(factorization_cca(x, y, n_components=50)["correlations"]) == 3
        assert len(integration_cca(_layers(x, y), layers=("rna", "protein"), n_components=50)[4]) == 3

    def test_lists_of_lists_are_accepted_by_the_factorization_version(self) -> None:
        x, y = _random_layers(20, 4, 3)

        from_lists = factorization_cca(x.tolist(), y.tolist(), n_components=2, regularization=0.3)
        from_arrays = factorization_cca(x, y, n_components=2, regularization=0.3)

        assert from_lists == from_arrays
