"""Regression tests: Spearman correlation in ``cross_omic`` averages tied ranks.

Defect: ``_to_ranks`` gave tied values consecutive distinct ranks in list order (three zeros became ranks 1, 2, 3).
The result then depended on the order in which the samples happened to be listed, and on data with many ties --
expression matrices with many zero counts, genotype dosages 0/1/2 -- the coefficient was not Spearman's at all:
``scipy.stats.spearmanr`` gives the tied values the mean of the ranks they span.

Oracles: ranks worked out by hand, ``scipy.stats.rankdata`` and ``scipy.stats.spearmanr`` on the same data, and
the defining property that a rank statistic cannot depend on the order of the samples.
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pytest
from scipy import stats

from metainformant.phenotype.integration.cross_omic import (
    _spearman_correlation,
    _to_ranks,
    phenotype_genotype_association,
    trait_expression_correlation,
)


def zero_heavy_pair(n: int, zero_fraction: float, seed: int) -> Tuple[List[float], List[float]]:
    """A continuous trait and an expression vector in which about ``zero_fraction`` of the samples are exactly zero.

    The first four samples are never zero, so the vector is never constant (scipy refuses a constant input).
    """
    rng = np.random.default_rng(seed)
    trait = rng.normal(size=n)
    expression = np.exp(0.7 * trait + rng.normal(size=n))
    zeroed = rng.random(n) < zero_fraction
    zeroed[:4] = False
    expression[zeroed] = 0.0
    return trait.tolist(), expression.tolist()


# ---------------------------------------------------------------------------
# The ranks themselves
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([3.0, 1.0, 3.0, 2.0, 3.0], [4.0, 1.0, 4.0, 2.0, 4.0]),  # three 3s span ranks 3, 4, 5 -> 4 each
        ([5.0, 5.0, 5.0, 5.0], [2.5, 2.5, 2.5, 2.5]),  # all tied: mean of 1..4
        ([2.0, 1.0, 2.0, 1.0], [3.5, 1.5, 3.5, 1.5]),  # two pairs of ties
        ([0.0, 0.0, 0.0, 7.0, 0.0, 9.0], [2.5, 2.5, 2.5, 5.0, 2.5, 6.0]),  # four zeros span ranks 1-4
        ([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]),  # no ties: plain ranks
        ([3.0, 2.0, 1.0], [3.0, 2.0, 1.0]),
        ([7.0], [1.0]),
        ([], []),
    ],
)
def test_ranks_match_the_hand_computed_average_ranks(values: List[float], expected: List[float]) -> None:
    assert _to_ranks(values) == pytest.approx(expected)


@pytest.mark.parametrize("zero_fraction", [0.0, 0.3, 0.55, 0.8])
def test_ranks_equal_scipy_rankdata_and_always_sum_to_n_times_n_plus_one_over_two(zero_fraction: float) -> None:
    _, expression = zero_heavy_pair(60, zero_fraction, seed=5)

    ranks = _to_ranks(expression)

    assert ranks == pytest.approx(stats.rankdata(expression).tolist())
    assert sum(ranks) == pytest.approx(60 * 61 / 2)  # the rank total does not change when values tie


def test_a_value_has_the_same_rank_wherever_it_sits_in_the_list() -> None:
    values = [0.0, 4.0, 0.0, 2.0, 0.0, 9.0, 4.0]
    expected = {0.0: 2.0, 2.0: 4.0, 4.0: 5.5, 9.0: 7.0}
    shuffled = [values[i] for i in (4, 6, 2, 0, 5, 1, 3)]

    for ordering in (values, shuffled, values[::-1]):
        assert _to_ranks(ordering) == pytest.approx([expected[v] for v in ordering])


# ---------------------------------------------------------------------------
# Spearman's coefficient and its p-value
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("zero_fraction", [0.0, 0.3, 0.55, 0.8])
@pytest.mark.parametrize("n", [12, 40, 120])
def test_spearman_matches_scipy_on_zero_heavy_expression(n: int, zero_fraction: float) -> None:
    trait, expression = zero_heavy_pair(n, zero_fraction, seed=100 + n)

    rho, p = _spearman_correlation(trait, expression)

    oracle = stats.spearmanr(trait, expression)
    assert rho == pytest.approx(oracle.statistic, rel=1e-9, abs=1e-12)
    assert p == pytest.approx(oracle.pvalue, rel=1e-6, abs=1e-12)


def test_spearman_matches_scipy_when_both_variables_have_many_ties() -> None:
    """Genotype dosages 0/1/2 against zero-heavy expression: ties on both sides."""
    rng = np.random.default_rng(8)
    dosage = rng.integers(0, 3, size=50).astype(float)
    expression = np.where(rng.random(50) < 0.5, 0.0, np.exp(0.6 * dosage + rng.normal(size=50)))

    rho, p = _spearman_correlation(dosage.tolist(), expression.tolist())

    oracle = stats.spearmanr(dosage, expression)
    assert rho == pytest.approx(oracle.statistic, rel=1e-9)
    assert p == pytest.approx(oracle.pvalue, rel=1e-6)


def test_spearman_does_not_depend_on_the_order_of_the_samples() -> None:
    trait, expression = zero_heavy_pair(40, 0.55, seed=21)
    reference = _spearman_correlation(trait, expression)

    for seed in range(5):
        order = np.random.default_rng(seed).permutation(40)
        permuted = _spearman_correlation([trait[i] for i in order], [expression[i] for i in order])

        assert permuted == pytest.approx(reference, rel=1e-12)


def test_a_constant_variable_has_no_rank_correlation() -> None:
    assert _spearman_correlation([1.0, 2.0, 3.0, 4.0], [5.0, 5.0, 5.0, 5.0]) == (0.0, 1.0)
    assert _spearman_correlation([0.0] * 6, [0.0] * 6) == (0.0, 1.0)


def test_a_perfect_monotone_relationship_with_ties_is_one() -> None:
    x = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0]
    y = [10.0, 10.0, 20.0, 20.0, 30.0, 30.0, 40.0]  # the same ties in the same places

    rho, p = _spearman_correlation(x, y)

    assert rho == pytest.approx(1.0)
    assert p == 0.0


# ---------------------------------------------------------------------------
# The public function that reports it
# ---------------------------------------------------------------------------


def test_trait_expression_correlation_reports_scipy_values_for_zero_heavy_genes() -> None:
    rng = np.random.default_rng(31)
    samples = [f"s{i}" for i in range(45)]
    trait = rng.normal(size=45)
    traits = dict(zip(samples, trait.tolist()))
    genes = {}
    for name, zero_fraction in (("sparse", 0.8), ("moderate", 0.4), ("dense", 0.0)):
        values = np.exp(0.8 * trait + rng.normal(size=45))
        values[rng.random(45) < zero_fraction] = 0.0
        genes[name] = dict(zip(samples, values.tolist()))

    result = trait_expression_correlation(traits, genes, method="spearman")

    for name, values in genes.items():
        oracle = stats.spearmanr([traits[s] for s in samples], [values[s] for s in samples])
        reported = result["correlations"][name]
        assert reported["correlation"] == pytest.approx(oracle.statistic, rel=1e-9)
        assert reported["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)


def test_pearson_is_unaffected_by_the_ranking_change() -> None:
    trait, expression = zero_heavy_pair(30, 0.5, seed=3)
    samples = [f"s{i}" for i in range(30)]

    result = trait_expression_correlation(
        dict(zip(samples, trait)), {"g": dict(zip(samples, expression))}, method="pearson"
    )

    oracle = stats.pearsonr(trait, expression)
    assert result["correlations"]["g"]["correlation"] == pytest.approx(oracle.statistic, rel=1e-9)
    assert result["correlations"]["g"]["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)


def test_genotype_association_methods_do_not_use_ranks() -> None:
    """``phenotype_genotype_association`` offers Pearson and regression only; the dosage ties do not touch it."""
    rng = np.random.default_rng(2)
    phenotypes = {f"s{i}": float(v) for i, v in enumerate(rng.normal(size=30))}
    calls = [int(g) for g in rng.integers(0, 3, size=30)]

    result = phenotype_genotype_association(phenotypes, {"v": calls}, method="correlation")

    oracle = stats.pearsonr([float(c) for c in calls], list(phenotypes.values()))
    assert result["associations"]["v"]["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)
