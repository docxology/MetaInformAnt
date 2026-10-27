"""The hypergeometric p-value of the metabolite enrichment is exact, or an error, never a fabricated 0.0.

Regression tests for ``metainformant.metabolomics.pathways.enrichment._hypergeometric_pvalue``. It wrapped each term
of the tail sum in ``except (ValueError, ZeroDivisionError): continue``. When the parameters did not describe a draw
from a population (a ``background_size`` smaller than a pathway or than the query) every term raised, was skipped
and the sum stayed at 0.0: the most significant p-value possible, reported for an invalid input.

Oracles independent of the implementation: ``scipy.stats.hypergeom.sf``, ``scipy.stats.fisher_exact`` and a
hand-computed fraction.

Real computation; no test doubles.
"""

from __future__ import annotations

import pytest
from scipy.stats import fisher_exact, hypergeom

from metainformant.metabolomics.pathways.enrichment import (
    _hypergeometric_pvalue,
    enrichment_with_fdr,
    metabolite_set_enrichment,
)


def _scipy_upper_tail(k: int, K: int, n: int, N: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws), by scipy."""
    return float(hypergeom.sf(k - 1, N, K, n))


class TestHypergeometricPvalue:
    def test_hand_computed_tail(self) -> None:
        # N=10, K=4, n=3: P(X >= 2) = (C(4,2)*C(6,1) + C(4,3)*C(6,0)) / C(10,3) = (36 + 4) / 120 = 1/3
        assert _hypergeometric_pvalue(2, 4, 3, 10) == pytest.approx(1 / 3, rel=1e-15)

    def test_matches_scipy_over_a_grid_of_parameters(self) -> None:
        checked = 0
        for N in (12, 40, 200):
            for K in sorted({0, 1, 5, N // 2, N}):
                for n in sorted({0, 1, 7, N // 3, N}):
                    for k in range(0, min(K, n) + 1):
                        ours = _hypergeometric_pvalue(k, K, n, N)
                        assert ours == pytest.approx(_scipy_upper_tail(k, K, n, N), rel=1e-12, abs=1e-300)
                        checked += 1
        assert checked > 300

    def test_matches_scipy_fisher_exact_on_two_by_two_tables(self) -> None:
        for k, K, n, N in [(4, 20, 15, 100), (0, 20, 15, 100), (9, 10, 12, 60), (5, 30, 5, 31)]:
            table = [[k, K - k], [n - k, N - K - n + k]]
            _, expected = fisher_exact(table, alternative="greater")

            assert _hypergeometric_pvalue(k, K, n, N) == pytest.approx(expected, rel=1e-10)

    @pytest.mark.parametrize(
        ("k", "K", "n", "N"),
        [(60, 300, 500, 20000), (40, 300, 500, 20000), (200, 1000, 1000, 5000)],
    )
    def test_stays_accurate_for_large_universes_and_tiny_tails(self, k: int, K: int, n: int, N: int) -> None:
        assert _hypergeometric_pvalue(k, K, n, N) == pytest.approx(_scipy_upper_tail(k, K, n, N), rel=1e-12)

    def test_k_zero_is_certain(self) -> None:
        assert _hypergeometric_pvalue(0, 5, 3, 20) == 1.0

    def test_an_overlap_that_cannot_occur_has_probability_zero(self) -> None:
        # more overlaps than pathway members: a true 0.0 (an impossible event), not a swallowed error
        assert _hypergeometric_pvalue(6, 5, 8, 20) == 0.0

    def test_an_overlap_forced_by_the_pigeonhole_principle_is_certain(self) -> None:
        # 8 draws from 10 items of which 5 are pathway members: at least 3 must be members
        assert _hypergeometric_pvalue(3, 5, 8, 10) == pytest.approx(1.0, rel=1e-15)

    @pytest.mark.parametrize(
        ("k", "K", "n", "N"),
        [
            (2, 5, 2, 3),  # more pathway members than the universe
            (1, 2, 10, 4),  # more draws than the universe
            (0, -1, 2, 5),
            (0, 2, -1, 5),
            (0, 0, 0, -1),
        ],
    )
    def test_parameters_that_are_not_a_draw_from_a_population_raise(self, k: int, K: int, n: int, N: int) -> None:
        with pytest.raises(ValueError, match=r"inconsistent hypergeometric parameters"):
            _hypergeometric_pvalue(k, K, n, N)


class TestEnrichmentBackgroundSize:
    PATHWAYS = {"p": ["a", "b", "c", "d", "e"], "q": ["f", "g"]}

    def test_a_background_smaller_than_a_pathway_raises_instead_of_reporting_p_zero(self) -> None:
        with pytest.raises(ValueError, match=r"background_size must be at least as large as the query"):
            metabolite_set_enrichment(["a", "b"], self.PATHWAYS, background_size=3)

    def test_a_background_smaller_than_the_query_raises(self) -> None:
        query = ["a", "b", "c", "f", "g", "x", "y", "z"]

        with pytest.raises(ValueError, match=r"inconsistent hypergeometric parameters"):
            metabolite_set_enrichment(query, self.PATHWAYS, background_size=6)

    def test_a_negative_background_raises(self) -> None:
        with pytest.raises(ValueError, match=r"inconsistent hypergeometric parameters"):
            metabolite_set_enrichment(["a"], self.PATHWAYS, background_size=-5)

    def test_the_fdr_wrapper_raises_too(self) -> None:
        with pytest.raises(ValueError, match=r"background_size"):
            enrichment_with_fdr(["a", "b"], self.PATHWAYS, background_size=3)

    def test_the_default_universe_always_contains_the_query_and_every_pathway(self) -> None:
        results = metabolite_set_enrichment(["a", "b", "never_in_a_pathway"], self.PATHWAYS)

        assert all(0.0 <= r.p_value <= 1.0 for r in results)

    def test_p_values_equal_fisher_exact_for_an_explicit_background(self) -> None:
        pathways = {"glycolysis": [f"g{i}" for i in range(12)], "tca": [f"t{i}" for i in range(8)]}
        query = ["g0", "g1", "g2", "g3", "t0", "x1", "x2"]
        universe = 150

        for result in metabolite_set_enrichment(query, pathways, background_size=universe):
            K = len(pathways[result.pathway_name])
            n = len(set(query))
            k = result.overlap
            table = [[k, K - k], [n - k, universe - K - n + k]]
            _, expected = fisher_exact(table, alternative="greater")

            assert result.p_value == pytest.approx(expected, rel=1e-10)

    def test_a_strongly_enriched_pathway_is_significant_and_a_disjoint_one_is_not(self) -> None:
        pathways = {"hit": [f"h{i}" for i in range(10)], "miss": [f"m{i}" for i in range(10)]}
        query = [f"h{i}" for i in range(8)]

        by_name = {r.pathway_name: r for r in metabolite_set_enrichment(query, pathways, background_size=500)}

        assert by_name["hit"].p_value < 1e-10
        assert by_name["miss"].p_value == pytest.approx(1.0)
