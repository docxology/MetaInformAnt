"""Regression tests for the log-factorials behind ``over_representation_analysis``.

``_log_factorial`` summed logarithms up to n = 20 and switched to Stirling's approximation above that. The
approximation is off by about ``1 / (12 n)`` in the logarithm, and a hypergeometric term needs several of them, so
every over-representation p-value of a gene set with more than 20 genes was a percent or so away from the exact value,
and the curve had a step at n = 20 / 21. The function is now ``lgamma(n + 1)``.

The oracles are ``math.log(math.factorial(n))`` (exact integer arithmetic), ``math.log(math.comb(n, k))`` and
``scipy.stats.hypergeom.sf``. No test doubles.
"""

from __future__ import annotations

import math
import random

import pytest
from scipy import stats

from metainformant.ontology.pathway_enrichment import enrichment
from metainformant.ontology.pathway_enrichment.enrichment import over_representation_analysis


def test_log_factorial_matches_the_logarithm_of_the_exact_factorial_up_to_300() -> None:
    for n in range(0, 301):
        assert enrichment._log_factorial(n) == pytest.approx(math.log(math.factorial(n)), rel=1e-13, abs=1e-13), n


def test_log_factorial_has_no_step_where_stirling_used_to_take_over() -> None:
    for n in (19, 20, 21, 22):
        assert enrichment._log_factorial(n + 1) - enrichment._log_factorial(n) == pytest.approx(
            math.log(n + 1), rel=1e-12
        )


def test_log_comb_matches_the_logarithm_of_the_exact_binomial_coefficient() -> None:
    rng = random.Random(3)
    for _ in range(300):
        n = rng.randint(0, 400)
        k = rng.randint(0, n)
        assert enrichment._log_comb(n, k) == pytest.approx(math.log(math.comb(n, k)), rel=1e-12, abs=1e-12), (n, k)
    assert enrichment._log_comb(5, 7) == -math.inf


def test_hypergeometric_sf_matches_scipy_for_sets_larger_than_twenty() -> None:
    rng = random.Random(11)
    for _ in range(300):
        population = rng.randint(25, 500)
        successes = rng.randint(21, population)
        draws = rng.randint(21, population)
        k = rng.randint(1, min(successes, draws))

        expected = float(stats.hypergeom.sf(k - 1, population, successes, draws))

        assert enrichment._hypergeometric_sf(k, population, successes, draws) == pytest.approx(
            expected, rel=1e-9, abs=1e-300
        ), (k, population, successes, draws)


def _oracle_p_values(
    query: list[str], gene_sets: dict[str, list[str]], background: list[str] | None
) -> dict[str, float]:
    universe = set(background) if background is not None else set().union(*gene_sets.values()) | set(query)
    drawn = set(query) & universe
    p_values = {}
    for term, genes in gene_sets.items():
        members = set(genes) & universe
        if not members:
            continue
        overlap = len(drawn & members)
        p_values[term] = (
            1.0 if overlap == 0 else float(stats.hypergeom.sf(overlap - 1, len(universe), len(members), len(drawn)))
        )
    return p_values


@pytest.mark.parametrize("with_background", [False, True])
def test_over_representation_p_values_match_scipy_for_realistic_set_sizes(with_background: bool) -> None:
    rng = random.Random(5)
    genes = [f"g{i:04d}" for i in range(2000)]
    query = rng.sample(genes, 150)
    gene_sets = {}
    for index in range(12):
        enriched = rng.sample(query, rng.randint(0, 25))
        gene_sets[f"T{index}"] = list(set(enriched + rng.sample(genes, rng.randint(25, 300))))
    background = genes if with_background else None

    results = over_representation_analysis(query, gene_sets, background=background, correction="none")

    expected = _oracle_p_values(query, gene_sets, background)
    assert {r["term_id"] for r in results} == set(expected)
    for result in results:
        assert result["p_value"] == pytest.approx(expected[result["term_id"]], rel=1e-9, abs=1e-300), result["term_id"]
