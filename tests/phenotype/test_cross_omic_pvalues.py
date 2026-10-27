"""Regression tests: cross-omic association p-values come from the exact t and F distributions.

Defect: ``cross_omic._t_cdf`` and ``_f_cdf`` were hand-written approximations ("Very rough approximation for small
df"): the t CDF fell back to the normal distribution for more than 30 degrees of freedom (a two-sided p of 0.045
where Student's t gives 0.054 at 31 degrees of freedom) and used a formula that is not the t CDF below that, and
``_f_cdf`` is exact only for two numerator degrees of freedom while the regression test uses one. Both are replaced
by ``scipy.stats.t.sf`` and ``scipy.stats.f.sf``, evaluated in the survival form so that tiny p-values are
not rounded to 0 by ``1 - cdf``.

Oracles: ``scipy.stats.pearsonr``, ``scipy.stats.linregress`` and ``scipy.stats.spearmanr`` on the same data.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np
import pytest
from scipy import stats

from metainformant.phenotype.integration.cross_omic import (
    _linear_regression_stats,
    _pearson_correlation,
    _spearman_correlation,
    phenotype_environment_interaction,
    phenotype_genotype_association,
    trait_expression_correlation,
)

SAMPLE_SIZES = [3, 4, 5, 8, 12, 20, 31, 32, 33, 60, 150]  # degrees of freedom 1 .. 148: both old code branches


def _correlated_pair(n: int, rho: float, seed: int) -> Tuple[List[float], List[float]]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    y = rho * x + math.sqrt(1.0 - rho**2) * rng.normal(size=n)
    return x.tolist(), y.tolist()


@pytest.mark.parametrize("n", SAMPLE_SIZES)
@pytest.mark.parametrize("rho", [0.0, 0.3, 0.7])
def test_pearson_p_value_matches_scipy_at_every_sample_size(n: int, rho: float) -> None:
    x, y = _correlated_pair(n, rho, seed=1000 + n)

    r, p = _pearson_correlation(x, y)

    oracle = stats.pearsonr(x, y)
    assert r == pytest.approx(oracle.statistic, rel=1e-9, abs=1e-12)
    assert p == pytest.approx(oracle.pvalue, rel=1e-6, abs=1e-12)


def test_pearson_p_value_at_31_degrees_of_freedom_is_students_t_not_the_normal() -> None:
    """The old code switched to the normal distribution above 30 degrees of freedom; the two differ by ~20%."""
    n, df = 33, 31
    r_target = 2.0 / math.sqrt(4.0 + df)  # the correlation that gives a t statistic of exactly 2.0
    x = np.arange(n, dtype=float)
    noise = np.random.default_rng(5).normal(size=n)
    noise = noise - np.polyval(np.polyfit(x, noise, 1), x)  # orthogonal to x and to the constant
    # Two orthogonal unit-variance vectors mixed with weights r and sqrt(1 - r^2) have correlation exactly r.
    y = r_target * (x - x.mean()) / x.std() + math.sqrt(1 - r_target**2) * noise / noise.std()

    r, p = _pearson_correlation(x.tolist(), y.tolist())

    assert r == pytest.approx(r_target, rel=1e-9)
    assert p == pytest.approx(2 * stats.t.sf(2.0, df), rel=1e-6)  # 0.0544
    assert p == pytest.approx(0.0544, abs=5e-4)
    assert abs(p - 2 * stats.norm.sf(2.0)) > 0.005  # the normal approximation (0.0455) is visibly different


def test_pearson_p_value_is_not_rounded_to_zero_in_the_far_tail() -> None:
    """``1 - cdf`` loses everything below ~1e-16; the survival function keeps the tail."""
    rng = np.random.default_rng(9)
    x = rng.normal(size=40)
    y = x + 0.05 * rng.normal(size=40)

    r, p = _pearson_correlation(x.tolist(), y.tolist())

    oracle = stats.pearsonr(x, y)
    assert 0.0 < oracle.pvalue < 1e-30
    assert p == pytest.approx(oracle.pvalue, rel=1e-6)


@pytest.mark.parametrize("n", SAMPLE_SIZES)
@pytest.mark.parametrize("rho", [0.0, 0.4, -0.8])
def test_linear_regression_p_value_matches_scipy_at_every_sample_size(n: int, rho: float) -> None:
    x, y = _correlated_pair(n, rho, seed=2000 + n)

    slope, intercept, r_squared, p = _linear_regression_stats(x, y)

    oracle = stats.linregress(x, y)
    assert slope == pytest.approx(oracle.slope, rel=1e-9, abs=1e-12)
    assert intercept == pytest.approx(oracle.intercept, rel=1e-9, abs=1e-12)
    assert r_squared == pytest.approx(oracle.rvalue**2, rel=1e-9, abs=1e-12)
    assert p == pytest.approx(oracle.pvalue, rel=1e-6, abs=1e-12)


def test_linear_regression_p_value_is_not_rounded_to_zero_in_the_far_tail() -> None:
    rng = np.random.default_rng(10)
    x = rng.normal(size=40)
    y = 2.0 * x + 0.05 * rng.normal(size=40)

    _, _, _, p = _linear_regression_stats(x.tolist(), y.tolist())

    oracle = stats.linregress(x, y)
    assert 0.0 < oracle.pvalue < 1e-30
    assert p == pytest.approx(oracle.pvalue, rel=1e-6)


def test_degenerate_inputs_keep_their_documented_results() -> None:
    assert _pearson_correlation([1.0, 2.0], [2.0, 4.0]) == (0.0, 1.0)  # fewer than three points
    assert _pearson_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == (0.0, 1.0)  # no variance
    r, p = _pearson_correlation([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0])  # perfect correlation
    assert r == pytest.approx(1.0) and p == 0.0
    assert _linear_regression_stats([1.0, 2.0], [1.0, 2.0]) == (0.0, 0.0, 0.0, 1.0)
    assert _linear_regression_stats([2.0, 2.0, 2.0], [1.0, 2.0, 4.0])[3] == 1.0  # constant predictor


def test_p_values_stay_within_zero_and_one() -> None:
    for n in SAMPLE_SIZES:
        x, y = _correlated_pair(n, 0.5, seed=n)
        assert 0.0 <= _pearson_correlation(x, y)[1] <= 1.0
        assert 0.0 <= _linear_regression_stats(x, y)[3] <= 1.0


@pytest.mark.parametrize("n", [5, 12, 40])
def test_spearman_p_value_matches_scipy_for_data_without_ties(n: int) -> None:
    x, y = _correlated_pair(n, 0.6, seed=3000 + n)

    rho, p = _spearman_correlation(x, y)

    oracle = stats.spearmanr(x, y)
    assert rho == pytest.approx(oracle.statistic, rel=1e-9, abs=1e-12)
    assert p == pytest.approx(oracle.pvalue, rel=1e-6, abs=1e-12)


# ---------------------------------------------------------------------------
# The public functions that report these p-values
# ---------------------------------------------------------------------------


def _association_inputs(n: int = 40) -> Tuple[Dict[str, float], Dict[str, List[int]]]:
    rng = np.random.default_rng(77)
    samples = [f"s{i}" for i in range(n)]
    genotypes = {f"v{k}": [int(g) for g in rng.integers(0, 3, n)] for k in range(3)}
    phenotype_values = 0.5 * np.array(genotypes["v0"]) + rng.normal(size=n)
    return dict(zip(samples, phenotype_values.tolist())), genotypes


def test_genotype_association_p_values_match_scipy_for_both_methods() -> None:
    phenotypes, genotypes = _association_inputs()
    y = list(phenotypes.values())

    regression = phenotype_genotype_association(phenotypes, genotypes, method="linear_regression")
    correlation = phenotype_genotype_association(phenotypes, genotypes, method="correlation")

    for variant, calls in genotypes.items():
        fit = stats.linregress([float(c) for c in calls], y)
        pearson = stats.pearsonr([float(c) for c in calls], y)
        assert regression["associations"][variant]["beta"] == pytest.approx(fit.slope, rel=1e-9)
        assert regression["associations"][variant]["p_value"] == pytest.approx(fit.pvalue, rel=1e-6)
        assert correlation["associations"][variant]["correlation"] == pytest.approx(pearson.statistic, rel=1e-9)
        assert correlation["associations"][variant]["p_value"] == pytest.approx(pearson.pvalue, rel=1e-6)


def test_genotype_association_skips_missing_calls_before_testing() -> None:
    phenotypes, genotypes = _association_inputs()
    calls = list(genotypes["v0"])
    calls[3] = calls[11] = -1  # two missing genotype calls
    y = list(phenotypes.values())
    kept = [i for i, call in enumerate(calls) if call >= 0]

    result = phenotype_genotype_association(phenotypes, {"v0": calls}, method="correlation")

    oracle = stats.pearsonr([float(calls[i]) for i in kept], [y[i] for i in kept])
    assert result["associations"]["v0"]["n_valid"] == 38
    assert result["associations"]["v0"]["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)


def test_trait_expression_correlation_p_values_match_scipy() -> None:
    rng = np.random.default_rng(31)
    samples = [f"s{i}" for i in range(25)]
    trait = rng.normal(size=25)
    traits = dict(zip(samples, trait.tolist()))
    expression = {
        "geneA": dict(zip(samples, (0.8 * trait + rng.normal(size=25)).tolist())),
        "geneB": dict(zip(samples, rng.normal(size=25).tolist())),
    }

    pearson = trait_expression_correlation(traits, expression, method="pearson")
    spearman = trait_expression_correlation(traits, expression, method="spearman")

    for gene, values in expression.items():
        x = [traits[s] for s in samples]
        y = [values[s] for s in samples]
        assert pearson["correlations"][gene]["p_value"] == pytest.approx(stats.pearsonr(x, y).pvalue, rel=1e-6)
        assert spearman["correlations"][gene]["p_value"] == pytest.approx(stats.spearmanr(x, y).pvalue, rel=1e-6)
    assert pearson["correlations"]["geneA"]["p_value"] < pearson["correlations"]["geneB"]["p_value"]


def test_gene_by_environment_interaction_p_value_matches_scipy() -> None:
    rng = np.random.default_rng(55)
    n = 36
    samples = [f"s{i}" for i in range(n)]
    calls = [int(g) for g in rng.integers(0, 3, n)]
    environment = rng.normal(size=n)
    outcome = 0.6 * np.array(calls) * environment + rng.normal(size=n)

    result = phenotype_environment_interaction(
        dict(zip(samples, outcome.tolist())),
        {"v1": calls},
        dict(zip(samples, environment.tolist())),
        interaction_model="multiplicative",
    )

    interaction = np.array(calls, dtype=float) * environment
    oracle = stats.pearsonr(interaction, outcome)
    assert result["interactions"]["v1"]["gxe_correlation"] == pytest.approx(oracle.statistic, rel=1e-9)
    assert result["interactions"]["v1"]["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)
