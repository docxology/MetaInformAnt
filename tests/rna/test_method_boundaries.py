"""Numerical counterexamples and independent oracles for RNA methods."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control
from sklearn.decomposition import PCA

from metainformant.rna.analysis.expression_analysis import (
    adjust_pvalues,
    differential_expression,
    pca_analysis,
)
from metainformant.rna.analysis.expression_core import (
    estimate_size_factors,
    normalize_counts,
)


@pytest.mark.parametrize("values", [[-0.1, 0.2], [0.1, 1.01], [0.1, np.inf], [[0.1, 0.2]]])
def test_adjustment_rejects_invalid_probabilities(values: list[float]) -> None:
    # Given invalid probabilities; when adjusted; then they cannot become evidence.
    with pytest.raises(ValueError):
        adjust_pvalues(np.asarray(values))


def test_adjustment_matches_independent_bh_oracle() -> None:
    # Given a complete family including ties and endpoints.
    p = np.array([0.01, 0.01, 0.08, 1.0, 0.0, 0.4, 0.03])
    # When corrected; then compare against SciPy's separate implementation.
    np.testing.assert_allclose(adjust_pvalues(p), false_discovery_control(p))


@pytest.mark.parametrize("method", ["ttest", "wilcoxon", "deseq2_like"])
def test_de_rejects_unreplicated_groups(method: str) -> None:
    # Given one sample per group.
    counts = pd.DataFrame({"a": [20.0, 100.0], "b": [30.0, 110.0]})
    # When tested; then no significance can be manufactured from unestimated variance.
    with pytest.raises(ValueError, match="at least two"):
        differential_expression(counts, ["A", "B"], method=method)


@pytest.mark.parametrize("function", [normalize_counts, estimate_size_factors])
@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_normalization_rejects_nonfinite_counts(function, bad: float) -> None:
    # Given a non-finite abundance.
    counts = pd.DataFrame({"a": [bad, 2.0], "b": [1.0, 3.0]})
    # When normalized; then missingness cannot be interpreted as abundance.
    with pytest.raises(ValueError):
        function(counts)


@pytest.mark.parametrize(
    "values",
    [
        [[1.0, 1.0], [2.0, 2.0]],
        [[1.0], [2.0]],
        [[np.nan, 2.0], [1.0, 3.0]],
        [[np.inf, 2.0], [1.0, 3.0]],
    ],
)
def test_pca_rejects_unidentifiable_or_nonfinite_input(
    values: list[list[float]],
) -> None:
    # Given undefined variation or unreported missingness.
    with pytest.raises(ValueError):
        pca_analysis(pd.DataFrame(values))


@pytest.mark.parametrize("n", [0, -1, 1.5, True])
def test_pca_rejects_invalid_component_count(n: int) -> None:
    # Given an invalid requested dimension.
    with pytest.raises(ValueError):
        pca_analysis(pd.DataFrame([[1.0, 2.0, 4.0], [3.0, 8.0, 2.0]]), n_components=n)


def test_pca_matches_sklearn_variance_and_reconstructs_centered_matrix() -> None:
    # Given finite expression with distinct sample and feature labels.
    frame = pd.DataFrame(np.random.default_rng(37).normal(size=(7, 5)))
    # When decomposed without scaling.
    result = pca_analysis(frame, n_components=5, scale=False)
    oracle = PCA(n_components=5, svd_solver="full").fit(frame.T)
    # Then variance and the reconstructed centered matrix match independently.
    np.testing.assert_allclose(result["explained_variance_ratio"], oracle.explained_variance_ratio_, atol=1e-12)
    reconstructed = result["transformed"].to_numpy() @ result["components"]
    np.testing.assert_allclose(reconstructed, frame.T.to_numpy() - frame.T.to_numpy().mean(axis=0), atol=1e-12)


def test_pca_mean_imputation_is_explicit_and_recorded() -> None:
    # Given one explicitly accepted missing entry.
    frame = pd.DataFrame([[1.0, np.nan, 3.0], [2.0, 4.0, 1.0]])
    # When the caller requests mean imputation.
    result = pca_analysis(frame, missing="mean")
    # Then its count and policy accompany the coordinates.
    assert result["preprocessing"]["imputed_values"] == 1
    assert result["preprocessing"]["missing"] == "mean"
    assert frame.isna().sum().sum() == 1


def test_size_factors_use_genes_positive_in_every_sample() -> None:
    # Given two stable genes plus a gene missing from one library.
    counts = pd.DataFrame({"a": [10.0, 0.0], "b": [40.0, 1000.0]})
    # When standard median-of-ratios factors are estimated.
    factors = estimate_size_factors(counts)
    # Then the incomplete gene cannot distort the reference geometric means.
    np.testing.assert_allclose(factors.to_numpy(), [0.5, 2.0])


@pytest.mark.parametrize("method", ["tpm", "rpkm"])
def test_length_normalization_rejects_infinite_lengths(method: str) -> None:
    # Given an infinite feature length.
    frame = pd.DataFrame({"a": [10.0, 20.0], "b": [40.0, 80.0]})
    # When length-normalized; then the feature cannot silently disappear.
    with pytest.raises(ValueError):
        normalize_counts(frame, method=method, gene_lengths=pd.Series([np.inf, 1000.0]))


def test_welch_one_constant_group_uses_variable_group_degrees_of_freedom() -> None:
    # Given unequal group sizes and one constant group; stable genes fix size factors at one.
    counts = pd.DataFrame(np.full((11, 5), 100.0), columns=["r1", "r2", "r3", "t1", "t2"])
    counts.loc[0] = [10.0, 20.0, 40.0, 80.0, 80.0]
    # When Welch's test is computed.
    result = differential_expression(counts, ["r", "r", "r", "t", "t"], method="ttest").set_index("gene").loc[0]
    # Then the variable reference group's n-1 degrees of freedom determine the p-value.
    from scipy.stats import t

    reference = np.log2(np.array([10.0, 20.0, 40.0]) + 1)
    treatment = np.log2(81.0)
    statistic = (treatment - reference.mean()) / np.sqrt(reference.var(ddof=1) / 3)
    np.testing.assert_allclose(result["stat"], statistic)
    np.testing.assert_allclose(result["p_value"], 2 * t.sf(abs(statistic), df=2))
