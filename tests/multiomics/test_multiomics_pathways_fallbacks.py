"""The scipy-free p-value helpers of the pathway module do what their docstrings say.

Regression tests for ``metainformant.multiomics.pathways.enrichment``. ``_py_chi2_sf`` is the fallback that
Fisher's combined-probability test uses when scipy is missing. Its docstring promised an exact treatment of even
degrees of freedom (all that Fisher's method needs: ``df = 2k``) but the code applied the Wilson-Hilferty cube-root
approximation to every ``df``, which is off by a factor of several for small ``df`` and small p-values: for one
p-value (``df = 2``) the combined p-value was not even equal to that p-value. The unmeasured-omic convention of
``multi_omic_enrichment`` (``p = 1.0``, still one of the combined omics) was undocumented.

The oracle is scipy's chi-squared survival function. Real computation; no test doubles (real-implementation
policy).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from metainformant.multiomics.pathways.enrichment import _fisher_combine, _py_chi2_sf, multi_omic_enrichment


class TestChiSquaredSurvivalWithoutScipy:
    @pytest.mark.parametrize("df", [2, 4, 6, 10, 40])
    @pytest.mark.parametrize("x", [0.01, 0.5, 2.0, 5.99, 20.0, 80.0, 300.0, 1380.0])
    def test_even_degrees_of_freedom_are_exact(self, x: float, df: int) -> None:
        expected = stats.chi2.sf(x, df)

        assert _py_chi2_sf(x, df) == pytest.approx(expected, rel=1e-10, abs=0.0)

    def test_one_p_value_combines_to_itself(self) -> None:
        # -2 ln p is chi-squared with 2 degrees of freedom, whose survival function is exp(-x/2) = p
        for p in (0.9, 0.5, 0.05, 1e-4, 1e-12):
            assert _py_chi2_sf(-2.0 * math.log(p), 2) == pytest.approx(p, rel=1e-12)

    def test_the_fallback_agrees_with_the_scipy_route_of_fisher_combine(self) -> None:
        p_values = [0.003, 0.2, 0.04, 0.9]
        statistic = -2.0 * sum(math.log(p) for p in p_values)

        assert _fisher_combine(p_values) == pytest.approx(_py_chi2_sf(statistic, 2 * len(p_values)), rel=1e-10)

    def test_a_huge_statistic_underflows_to_zero_without_overflowing(self) -> None:
        assert _py_chi2_sf(5000.0, 4) == 0.0
        assert 0.0 < _py_chi2_sf(1400.0, 4) < 1e-300

    def test_non_positive_statistics_and_degrees_of_freedom_keep_their_conventions(self) -> None:
        assert _py_chi2_sf(0.0, 4) == 1.0 and _py_chi2_sf(-3.0, 4) == 1.0
        assert _py_chi2_sf(5.0, 0) == 0.0

    @pytest.mark.parametrize("df", [1, 3, 5, 9])
    def test_odd_degrees_of_freedom_are_the_documented_approximation(self, df: int) -> None:
        for x in np.linspace(0.1, 40.0, 80):
            got, expected = _py_chi2_sf(float(x), df), float(stats.chi2.sf(x, df))
            assert abs(got - expected) < 0.01  # Wilson-Hilferty: close in absolute terms, not exact
            assert 0.0 <= got <= 1.0

    def test_the_docstring_says_which_degrees_of_freedom_are_exact(self) -> None:
        doc = " ".join((_py_chi2_sf.__doc__ or "").split())

        assert "Exact for even" in doc and "Wilson-Hilferty" in doc


class TestUnmeasuredOmicsAreEnteredAsNoEvidence:
    GENE_SETS = {"pathway": ["g1", "g2"]}

    def test_an_omic_without_the_pathway_genes_counts_as_p_one_in_the_combination(self) -> None:
        results = multi_omic_enrichment(
            self.GENE_SETS, {"rna": {"g1": 0.01, "g2": 0.01}, "protein": {"other": 0.5}}, method="fisher_combined"
        )

        (row,) = results
        rna_p = float(stats.chi2.sf(-2.0 * 2 * math.log(0.01), 4))  # Fisher over the pathway's two genes
        assert row["per_omic_p"] == {"rna": pytest.approx(rna_p), "protein": 1.0}
        # the unmeasured omic adds nothing to the statistic but does add two degrees of freedom
        assert row["combined_p"] == pytest.approx(float(stats.chi2.sf(-2.0 * math.log(rna_p), 4)), rel=1e-10)
        assert row["combined_p"] > rna_p

    def test_the_documentation_states_the_convention(self) -> None:
        doc = " ".join((multi_omic_enrichment.__doc__ or "").split())

        assert "entered as ``p = 1.0``" in doc and "Bonferroni factor" in doc

    def test_min_p_counts_the_unmeasured_omic_in_the_bonferroni_factor(self) -> None:
        (row,) = multi_omic_enrichment(
            self.GENE_SETS, {"rna": {"g1": 0.001, "g2": 0.5}, "protein": {"other": 0.5}}, method="min_p"
        )

        rna_p = float(stats.chi2.sf(-2.0 * (math.log(0.001) + math.log(0.5)), 4))
        assert row["combined_p"] == pytest.approx(min(1.0, min(rna_p, 1.0) * 2))
        assert np.isfinite(row["combined_p"])
