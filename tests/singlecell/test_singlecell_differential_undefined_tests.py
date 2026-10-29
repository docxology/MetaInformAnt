"""A differential-expression test that cannot be computed yields NaN and a status, never p = 1.

Before, a NaN p-value from scipy (constant data, one value in a group) was returned as 1.0, so was any
exception, and the invented 1.0 then entered the Benjamini-Hochberg family. Real data and real scipy calls
only; nothing is replaced by a test double, and scipy is only the oracle that the exact formulas are
checked against.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from scipy import stats as scipy_stats

from metainformant.singlecell.differential import expression
from metainformant.singlecell.differential.expression import (
    differential_expression,
    gene_set_scoring,
    pseudobulk_de,
    volcano_data,
)

GROUPS = [0] * 10 + [1] * 10


def _two_genes_plus_a_constant_one() -> tuple[np.ndarray, list[str]]:
    """Cells x genes: one gene that differs between the groups, one noise gene, one gene equal in every cell."""
    rng = np.random.default_rng(0)
    varied = rng.normal(5.0, 1.0, (20, 2))
    varied[10:, 0] += 3.0
    constant = np.full((20, 1), 4.0)
    return np.hstack([varied, constant]), ["g_de", "g_noise", "g_const"]


def _by_gene(results: list[dict]) -> dict[str, dict]:
    return {row["gene"]: row for row in results}


class TestUndefinedGenesGetNaNAndAStatus:
    @pytest.mark.parametrize("method", ["wilcoxon", "t_test"])
    def test_a_gene_that_is_constant_in_every_cell_has_no_p_value(self, method: str) -> None:
        matrix, names = _two_genes_plus_a_constant_one()

        rows = _by_gene(differential_expression(matrix, GROUPS, names, method=method))

        constant = rows["g_const"]
        assert math.isnan(constant["p_value"])
        assert math.isnan(constant["adjusted_p"])
        assert constant["test_status"].startswith("zero_variance")
        assert rows["g_de"]["test_status"] == "ok"
        assert 0.0 <= rows["g_de"]["p_value"] < 0.01

    @pytest.mark.parametrize("method", ["wilcoxon", "t_test"])
    def test_the_undefined_gene_is_outside_the_benjamini_hochberg_family(self, method: str) -> None:
        matrix, names = _two_genes_plus_a_constant_one()

        with_constant = _by_gene(differential_expression(matrix, GROUPS, names, method=method))
        without = _by_gene(differential_expression(matrix[:, :2], GROUPS, names[:2], method=method))

        for gene in ("g_de", "g_noise"):
            assert with_constant[gene]["adjusted_p"] == without[gene]["adjusted_p"]
            assert with_constant[gene]["p_value"] == without[gene]["p_value"]

    @pytest.mark.parametrize("method", ["wilcoxon", "t_test"])
    def test_genes_without_a_p_value_sort_after_every_gene_that_has_one(self, method: str) -> None:
        matrix, names = _two_genes_plus_a_constant_one()
        # Two undefined genes around a tested one, so the order cannot come out right by accident
        matrix = np.hstack([np.full((20, 1), 1.0), matrix, np.full((20, 1), 7.0)])
        names = ["c_first", *names, "c_last"]

        results = differential_expression(matrix, GROUPS, names, method=method)

        statuses = [row["test_status"] == "ok" for row in results]
        assert statuses == sorted(statuses, reverse=True)
        assert sum(not ok for ok in statuses) == 3
        defined = [row["adjusted_p"] for row in results if row["test_status"] == "ok"]
        assert defined == sorted(defined)

    def test_a_group_with_one_cell_has_no_t_test(self) -> None:
        rng = np.random.default_rng(1)
        matrix = rng.normal(5.0, 1.0, (10, 3))
        groups = [0] + [1] * 9

        results = differential_expression(matrix, groups, ["a", "b", "c"], method="t_test")

        assert len(results) == 3
        for row in results:
            assert math.isnan(row["p_value"]) and math.isnan(row["adjusted_p"])
            assert row["test_status"].startswith("insufficient_observations")

    def test_the_log_reports_the_genes_that_could_not_be_tested(self, caplog: pytest.LogCaptureFixture) -> None:
        matrix, names = _two_genes_plus_a_constant_one()

        with caplog.at_level("WARNING", logger=expression.logger.name):
            differential_expression(matrix, GROUPS, names, method="t_test")

        assert any("1 of 3 genes have no p-value" in record.getMessage() for record in caplog.records)

    def test_pseudobulk_with_one_sample_per_group_has_no_p_values(self) -> None:
        rng = np.random.default_rng(2)
        matrix = rng.poisson(5.0, (12, 4)).astype(float)
        sample_labels = ["S1"] * 6 + ["S2"] * 6
        groups = [0] * 6 + [1] * 6

        results = pseudobulk_de(matrix, ["T"] * 12, sample_labels, groups, gene_names=list("abcd"))

        assert len(results) == 4
        for row in results:
            assert math.isnan(row["p_value"]) and math.isnan(row["adjusted_p"])
            assert row["test_status"].startswith("insufficient_observations")
            assert row["normalization"] == "cpm"


class TestWilcoxonHelper:
    def test_an_empty_group_is_undefined(self) -> None:
        outcome = expression._wilcoxon_rank_sum([], [1.0, 2.0, 3.0])

        assert math.isnan(outcome.p_value) and math.isnan(outcome.statistic)
        assert outcome.status.startswith("empty_group")

    def test_all_tied_values_are_undefined_not_a_p_value_of_one(self) -> None:
        outcome = expression._wilcoxon_rank_sum([2.0] * 5, [2.0] * 4)

        assert math.isnan(outcome.p_value)
        assert outcome.status.startswith("zero_variance")

    def test_non_finite_values_are_undefined(self) -> None:
        outcome = expression._wilcoxon_rank_sum([1.0, math.nan, 3.0], [1.0, 2.0, 3.0])

        assert math.isnan(outcome.p_value)
        assert outcome.status.startswith("non_finite_values")

    def test_a_defined_test_reports_scipys_statistic_and_p_value(self) -> None:
        a = [5.0, 3.0, 8.0, 1.0, 9.0, 2.0, 7.0]
        b = [4.0, 6.0, 10.0, 12.0, 11.0, 13.0]

        outcome = expression._wilcoxon_rank_sum(a, b)

        expected = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
        assert outcome.status == "ok"
        assert outcome.p_value == pytest.approx(float(expected.pvalue), rel=1e-12)
        assert outcome.statistic == pytest.approx(float(expected.statistic))


class TestWelchHelper:
    @pytest.mark.parametrize(
        "n_a, n_b, shift, spread", [(5, 5, 0.0, 1.0), (6, 9, 1.0, 3.0), (30, 4, 0.5, 0.2), (200, 150, 0.1, 2.0)]
    )
    def test_matches_scipys_welch_test(self, n_a: int, n_b: int, shift: float, spread: float) -> None:
        rng = np.random.default_rng(n_a * 1000 + n_b)
        a = rng.normal(0.0, 1.0, n_a).tolist()
        b = rng.normal(shift, spread, n_b).tolist()

        outcome = expression._welch_t_test(a, b)

        expected = scipy_stats.ttest_ind(a, b, equal_var=False)
        assert outcome.status == "ok"
        assert outcome.statistic == pytest.approx(float(expected.statistic), rel=1e-10)
        assert outcome.p_value == pytest.approx(float(expected.pvalue), rel=1e-9)

    def test_one_constant_group_is_still_a_defined_welch_test(self) -> None:
        # The standard error comes from the other group alone; scipy computes it too, but warns about precision
        a = [2.0] * 6
        b = [1.0, 2.0, 3.5, 4.0, 5.0, 6.5]

        outcome = expression._welch_t_test(a, b)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            expected = scipy_stats.ttest_ind(a, b, equal_var=False)
        assert outcome.status == "ok"
        assert outcome.p_value == pytest.approx(float(expected.pvalue), rel=1e-9)

    @pytest.mark.parametrize("constant", [0.0, 1.0, 0.1, 0.3, 123456.789])
    def test_two_constant_groups_are_undefined_whatever_their_means(self, constant: float) -> None:
        same = expression._welch_t_test([constant] * 5, [constant] * 7)
        different = expression._welch_t_test([constant] * 5, [constant + 1.0] * 7)

        for outcome in (same, different):
            assert math.isnan(outcome.p_value) and math.isnan(outcome.statistic)
            assert outcome.status.startswith("zero_variance")

    @pytest.mark.parametrize("a, b", [([1.0], [1.0, 2.0, 3.0]), ([1.0, 2.0], []), ([], [])])
    def test_fewer_than_two_values_in_a_group_is_undefined(self, a: list[float], b: list[float]) -> None:
        outcome = expression._welch_t_test(a, b)

        assert math.isnan(outcome.p_value)
        assert outcome.status.startswith("insufficient_observations")

    def test_non_finite_values_are_undefined(self) -> None:
        outcome = expression._welch_t_test([1.0, 2.0, math.inf], [1.0, 2.0, 3.0])

        assert math.isnan(outcome.p_value)
        assert outcome.status.startswith("non_finite_values")


class TestExactTails:
    @pytest.mark.parametrize("df", [1.0, 2.0, 3.7, 10.0, 50.0, 1000.0, 1.0e5])
    @pytest.mark.parametrize("t", [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 40.0])
    def test_the_scipy_free_student_t_tail_matches_scipy(self, t: float, df: float) -> None:
        """The incomplete-beta evaluation used without scipy is exact, not a normal approximation."""
        x = df / (df + t * t)
        ours = expression._regularized_incomplete_beta(df / 2.0, 0.5, x)

        expected = float(2.0 * scipy_stats.t.sf(t, df))
        assert ours == pytest.approx(expected, rel=1e-8, abs=1e-300)

    @pytest.mark.parametrize("z", [0.0, 1.0, 1.96, 4.0, 6.0, 8.0, 12.0])
    def test_the_normal_cdf_is_accurate_in_the_far_tail(self, z: float) -> None:
        """A rational erf approximation is only good to ~1e-7 absolute, which is useless for tiny p-values."""
        assert expression._standard_normal_cdf(-z) == pytest.approx(float(scipy_stats.norm.cdf(-z)), rel=1e-12)


class TestBenjaminiHochbergWithMissingPValues:
    def test_nan_p_values_are_outside_the_family_and_stay_nan(self) -> None:
        raw = [0.01, math.nan, 0.04, 0.03, math.nan]

        adjusted = expression._benjamini_hochberg(raw)

        defined = [0, 2, 3]
        expected = scipy_stats.false_discovery_control([raw[i] for i in defined], method="bh")
        for position, value in zip(defined, expected):
            assert adjusted[position] == pytest.approx(float(value), rel=1e-12)
        assert math.isnan(adjusted[1]) and math.isnan(adjusted[4])

    def test_all_nan_gives_all_nan(self) -> None:
        adjusted = expression._benjamini_hochberg([math.nan, math.nan])

        assert len(adjusted) == 2 and all(math.isnan(v) for v in adjusted)

    def test_without_nan_it_is_the_ordinary_correction(self) -> None:
        rng = np.random.default_rng(3)
        raw = rng.uniform(0.0, 1.0, 40).tolist()

        adjusted = expression._benjamini_hochberg(raw)

        np.testing.assert_allclose(adjusted, scipy_stats.false_discovery_control(raw, method="bh"), rtol=1e-12)


class TestVolcanoDataWithoutAPValue:
    def test_a_nan_adjusted_p_is_not_significant_and_has_no_y_value(self) -> None:
        rows = [
            {"gene": "tested", "log2fc": 3.0, "adjusted_p": 0.001},
            {"gene": "untested", "log2fc": 3.0, "adjusted_p": math.nan},
        ]

        result = volcano_data(rows)

        assert result["classification"] == ["up", "ns"]
        assert result["neg_log10_p"][0] == pytest.approx(3.0)
        assert math.isnan(result["neg_log10_p"][1])
        assert (result["n_up"], result["n_down"], result["n_ns"]) == (1, 0, 1)

    def test_a_row_without_an_adjusted_p_is_not_given_a_p_value_of_one(self) -> None:
        result = volcano_data([{"gene": "g", "log2fc": 5.0}])

        assert result["classification"] == ["ns"]
        assert math.isnan(result["neg_log10_p"][0])  # not -log10(1) = 0.0, which would claim a p-value of exactly 1


class TestGeneSetsThatCannotBeScored:
    def test_a_gene_set_missing_from_the_matrix_scores_nan_for_every_method(self) -> None:
        matrix = np.arange(24, dtype=float).reshape(6, 4)

        for method in ("mean", "sum"):
            result = gene_set_scoring(
                matrix, {"absent": ["X1", "X2"], "present": ["a", "b"]}, list("abcd"), method=method, seed=0
            )

            assert result["gene_set_sizes"] == {"absent": 0, "present": 2}
            assert len(result["scores"]["absent"]) == 6
            assert all(math.isnan(score) for score in result["scores"]["absent"])
            assert not any(math.isnan(score) for score in result["scores"]["present"])


class TestLog2FoldChangeSign:
    """``group 1`` is the smaller group label: a gene higher in the cells labelled 1 has a negative log2fc."""

    def test_the_documented_sign_convention_is_the_behaviour(self) -> None:
        rng = np.random.default_rng(4)
        matrix = rng.normal(5.0, 0.5, (20, 1))
        matrix[10:, 0] += 4.0  # higher in the cells labelled 1

        row = differential_expression(matrix, GROUPS, ["g"], method="t_test")[0]

        assert row["log2fc"] == pytest.approx(math.log2((row["mean_group1"] + 1.0) / (row["mean_group2"] + 1.0)))
        assert row["log2fc"] < 0.0
        assert row["mean_group2"] > row["mean_group1"]
        documented = differential_expression.__doc__ or ""
        assert "NEGATIVE" in documented
        assert "group 1 vs group 0" not in documented
