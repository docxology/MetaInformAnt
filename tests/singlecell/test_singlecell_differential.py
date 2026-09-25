"""Tests for single-cell differential expression analysis.

Real implementation testing for DE statistical tests, pseudobulk DE,
fold-change computation, volcano plot data, and gene set scoring.
Real implementationing used - all tests use real computational methods and data.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.singlecell.differential.expression import (
    compute_log_fold_change,
    differential_expression,
    gene_set_scoring,
    pseudobulk_de,
    volcano_data,
)
from tests._support.synth import make_de_matrix

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_de_data(
    n_cells: int = 80,
    n_genes: int = 100,
    n_de_genes: int = 10,
    seed: int = 42,
) -> tuple[list[list[float]], list[int], list[str]]:
    """Expression matrix with known DE genes (shared factory)."""
    return make_de_matrix(
        n_cells=n_cells, n_genes=n_genes, n_de_genes=n_de_genes, seed=seed
    )


# ---------------------------------------------------------------------------
# differential_expression
# ---------------------------------------------------------------------------


class TestDifferentialExpression:
    """Tests for the differential_expression function."""

    def test_wilcoxon_returns_sorted_results(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(matrix, groups, gene_names, method="wilcoxon")
        assert isinstance(results, list)
        assert len(results) > 0
        # Should be sorted by adjusted_p
        p_values = [r["adjusted_p"] for r in results]
        assert p_values == sorted(p_values)

    def test_t_test_returns_results(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(matrix, groups, gene_names, method="t_test")
        assert len(results) > 0

    def test_result_dict_keys(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(matrix, groups, gene_names, method="wilcoxon")
        expected_keys = {
            "gene",
            "log2fc",
            "p_value",
            "adjusted_p",
            "pct_group1",
            "pct_group2",
            "mean_group1",
            "mean_group2",
        }
        for r in results:
            assert set(r.keys()) == expected_keys

    def test_de_genes_have_high_foldchange(self) -> None:
        """Known DE genes should appear with large absolute log2fc."""
        matrix, groups, gene_names = _make_de_data(n_de_genes=5)
        results = differential_expression(matrix, groups, gene_names, method="wilcoxon")
        de_gene_names = {f"gene_{i}" for i in range(5)}
        # Group 1 is upregulated; FC direction depends on which group is
        # the numerator in the implementation (g1_val vs g2_val).  We just
        # check absolute fold change.
        found_genes = {r["gene"] for r in results if abs(r["log2fc"]) > 0.5}
        overlap = de_gene_names & found_genes
        # At least some of the known DE genes should appear
        assert len(overlap) >= 1, (
            f"Expected DE genes {de_gene_names} in results; found {found_genes}"
        )

    def test_adjusted_p_values_bounded(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(matrix, groups, gene_names, method="wilcoxon")
        for r in results:
            assert 0.0 <= r["adjusted_p"] <= 1.0
            assert 0.0 <= r["p_value"] <= 1.0

    def test_pct_group_range(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(matrix, groups, gene_names, method="t_test")
        for r in results:
            assert 0.0 <= r["pct_group1"] <= 100.0
            assert 0.0 <= r["pct_group2"] <= 100.0

    def test_min_log2fc_filter(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(
            matrix, groups, gene_names, method="wilcoxon", min_log2fc=2.0
        )
        for r in results:
            assert abs(r["log2fc"]) >= 2.0

    def test_invalid_method_raises(self) -> None:
        matrix, groups, gene_names = _make_de_data()
        with pytest.raises(ValueError, match="Invalid method"):
            differential_expression(matrix, groups, gene_names, method="pseudobulk")

    def test_groups_length_mismatch_raises(self) -> None:
        matrix, _, gene_names = _make_de_data(n_cells=20)
        with pytest.raises(ValueError, match="must match"):
            differential_expression(matrix, [0, 1], gene_names)

    def test_gene_names_mismatch_raises(self) -> None:
        matrix, groups, _ = _make_de_data(n_cells=20)
        with pytest.raises(ValueError, match="must match"):
            differential_expression(matrix, groups, ["only_one"])

    def test_more_than_two_groups_raises(self) -> None:
        matrix = [[1.0, 2.0]] * 9
        groups = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        with pytest.raises(ValueError, match="exactly 2 unique"):
            differential_expression(matrix, groups, ["a", "b"])

    def test_numpy_array_input(self) -> None:
        matrix, groups, gene_names = _make_de_data(n_cells=40)
        arr = np.array(matrix)
        results = differential_expression(arr, groups, gene_names, method="wilcoxon")
        assert len(results) > 0

    def test_min_cells_filter(self) -> None:
        """Genes expressed in very few cells in BOTH groups are filtered.

        The filter skips a gene only when both groups have fewer expressing
        cells than min_cells.  We create a gene where neither group reaches
        the threshold.
        """
        rng = np.random.RandomState(11)
        n = 40
        # Baseline: all zeros (no expression anywhere)
        matrix = np.zeros((n, 10))
        # Genes 1-9: expressed in many cells so they pass
        matrix[:, 1:] = rng.exponential(2.0, size=(n, 9))
        # Gene 0: expressed in only 1 cell per group (total 2)
        matrix[0, 0] = 5.0  # group 0
        matrix[20, 0] = 5.0  # group 1
        groups = [0] * 20 + [1] * 20
        gene_names = [f"g{i}" for i in range(10)]
        results = differential_expression(
            matrix.tolist(), groups, gene_names, min_cells=3
        )
        gene_set = {r["gene"] for r in results}
        # gene_0 has only 1 expressing cell per group (< 3), should be skipped
        assert "g0" not in gene_set

    def test_all_zero_expression_returns_empty(self) -> None:
        """Genes with no expressing cells anywhere are filtered by min_cells."""
        matrix = np.zeros((20, 5)).tolist()
        groups = [0] * 10 + [1] * 10
        results = differential_expression(matrix, groups, [f"g{i}" for i in range(5)])
        assert results == []

    def test_adjusted_p_at_least_raw_p(self) -> None:
        """Benjamini-Hochberg adjusted p-values never fall below raw ones."""
        matrix, groups, gene_names = _make_de_data()
        results = differential_expression(matrix, groups, gene_names, method="wilcoxon")
        for r in results:
            assert r["adjusted_p"] >= r["p_value"] - 1e-12


# ---------------------------------------------------------------------------
# pseudobulk_de
# ---------------------------------------------------------------------------


class TestPseudobulkDE:
    """Tests for pseudobulk differential expression."""

    def _make_pseudobulk_data(
        self,
    ) -> tuple[list[list[float]], list[str], list[str], list[int], list[str]]:
        rng = np.random.RandomState(7)
        n_cells = 120
        n_genes = 50
        gene_names = [f"gene_{i}" for i in range(n_genes)]

        # 4 samples, 2 per group
        samples = ["S1", "S2", "S3", "S4"]
        sample_groups = {"S1": 0, "S2": 0, "S3": 1, "S4": 1}
        cells_per_sample = n_cells // 4

        matrix = rng.exponential(1.0, size=(n_cells, n_genes))
        cell_labels_list: list[str] = []
        sample_labels_list: list[str] = []
        group_labels_list: list[int] = []

        for s_idx, s_name in enumerate(samples):
            start = s_idx * cells_per_sample
            end = start + cells_per_sample
            sample_labels_list.extend([s_name] * cells_per_sample)
            cell_labels_list.extend(["T_cell"] * cells_per_sample)
            group_labels_list.extend([sample_groups[s_name]] * cells_per_sample)
            # Add group effect to first 5 genes in group 1
            if sample_groups[s_name] == 1:
                matrix[start:end, :5] += rng.exponential(
                    3.0, size=(cells_per_sample, 5)
                )

        return (
            matrix.tolist(),
            cell_labels_list,
            sample_labels_list,
            group_labels_list,
            gene_names,
        )

    def test_pseudobulk_returns_results(self) -> None:
        matrix, cell_labels, sample_labels, groups, gene_names = (
            self._make_pseudobulk_data()
        )
        results = pseudobulk_de(
            matrix,
            cell_labels,
            sample_labels,
            groups,
            gene_names=gene_names,
        )
        assert isinstance(results, list)
        assert len(results) > 0

    def test_pseudobulk_result_keys(self) -> None:
        matrix, cell_labels, sample_labels, groups, gene_names = (
            self._make_pseudobulk_data()
        )
        results = pseudobulk_de(
            matrix,
            cell_labels,
            sample_labels,
            groups,
            gene_names=gene_names,
        )
        for r in results:
            assert "gene" in r
            assert "log2fc" in r
            assert "adjusted_p" in r

    def test_pseudobulk_sorted_by_adjusted_p(self) -> None:
        matrix, cell_labels, sample_labels, groups, gene_names = (
            self._make_pseudobulk_data()
        )
        results = pseudobulk_de(
            matrix,
            cell_labels,
            sample_labels,
            groups,
            gene_names=gene_names,
        )
        adj_p = [r["adjusted_p"] for r in results]
        assert adj_p == sorted(adj_p)

    def test_dimension_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            pseudobulk_de(
                [[1.0, 2.0]],
                ["a"],
                ["s1"],
                [0, 1],  # wrong length
                gene_names=["g1", "g2"],
            )


# ---------------------------------------------------------------------------
# compute_log_fold_change
# ---------------------------------------------------------------------------


class TestLogFoldChange:
    """Tests for compute_log_fold_change."""

    def test_equal_means(self) -> None:
        assert compute_log_fold_change(5.0, 5.0) == 0.0

    def test_double_expression(self) -> None:
        # (3+1)/(1+1) = 2.0, log2(2) = 1.0
        fc = compute_log_fold_change(3.0, 1.0)
        assert abs(fc - 1.0) < 1e-9

    def test_zero_means(self) -> None:
        # (0+1)/(0+1) = 1, log2(1) = 0
        assert compute_log_fold_change(0.0, 0.0) == 0.0

    def test_custom_pseudocount(self) -> None:
        fc = compute_log_fold_change(0.0, 0.0, pseudocount=0.5)
        assert fc == 0.0

    def test_negative_fold_change(self) -> None:
        fc = compute_log_fold_change(1.0, 3.0)
        assert fc < 0


# ---------------------------------------------------------------------------
# volcano_data
# ---------------------------------------------------------------------------


class TestVolcanoData:
    """Tests for volcano_data preparation."""

    def _make_de_results(self) -> list[dict]:
        return [
            {"gene": "up1", "log2fc": 2.0, "adjusted_p": 0.001},
            {"gene": "up2", "log2fc": 1.5, "adjusted_p": 0.01},
            {"gene": "down1", "log2fc": -2.5, "adjusted_p": 0.005},
            {"gene": "ns1", "log2fc": 0.1, "adjusted_p": 0.5},
            {"gene": "ns2", "log2fc": 1.2, "adjusted_p": 0.8},
        ]

    def test_volcano_returns_expected_keys(self) -> None:
        result = volcano_data(self._make_de_results())
        for key in [
            "genes",
            "log2fc",
            "neg_log10_p",
            "classification",
            "n_up",
            "n_down",
            "n_ns",
        ]:
            assert key in result

    def test_classification_counts(self) -> None:
        result = volcano_data(
            self._make_de_results(), fc_threshold=1.0, p_threshold=0.05
        )
        assert result["n_up"] == 2  # up1, up2
        assert result["n_down"] == 1  # down1
        assert result["n_ns"] == 2  # ns1, ns2

    def test_neg_log10_p_positive(self) -> None:
        result = volcano_data(self._make_de_results())
        for val in result["neg_log10_p"]:
            assert val >= 0

    def test_gene_order_preserved(self) -> None:
        de = self._make_de_results()
        result = volcano_data(de)
        assert result["genes"] == [r["gene"] for r in de]

    def test_strict_threshold(self) -> None:
        result = volcano_data(
            self._make_de_results(), fc_threshold=3.0, p_threshold=0.001
        )
        # Only up1 is a borderline candidate but it has fc=2.0 < 3.0
        assert result["n_up"] == 0
        assert result["n_down"] == 0

    def test_empty_results(self) -> None:
        result = volcano_data([])
        assert result["n_up"] == 0
        assert result["n_down"] == 0
        assert result["n_ns"] == 0


# ---------------------------------------------------------------------------
# gene_set_scoring
# ---------------------------------------------------------------------------


class TestGeneSetScoring:
    """Tests for gene_set_scoring function."""

    def _make_scoring_data(
        self,
    ) -> tuple[list[list[float]], dict[str, list[str]], list[str]]:
        rng = np.random.RandomState(55)
        n_cells = 50
        n_genes = 80
        gene_names = [f"g{i}" for i in range(n_genes)]
        matrix = rng.exponential(1.0, size=(n_cells, n_genes))
        # Elevate genes 0-4 in first 25 cells
        matrix[:25, :5] += 5.0
        gene_sets = {"pathway_A": gene_names[:5], "pathway_B": gene_names[10:15]}
        return matrix.tolist(), gene_sets, gene_names

    def test_scores_returned_per_set(self) -> None:
        matrix, gene_sets, gene_names = self._make_scoring_data()
        result = gene_set_scoring(matrix, gene_sets, gene_names)
        assert "pathway_A" in result["scores"]
        assert "pathway_B" in result["scores"]

    def test_score_length_matches_cells(self) -> None:
        matrix, gene_sets, gene_names = self._make_scoring_data()
        result = gene_set_scoring(matrix, gene_sets, gene_names)
        for gs_name, scores in result["scores"].items():
            assert len(scores) == 50

    def test_mean_method(self) -> None:
        matrix, gene_sets, gene_names = self._make_scoring_data()
        result = gene_set_scoring(matrix, gene_sets, gene_names, method="mean", seed=0)
        assert result["n_gene_sets"] == 2
        assert result["n_cells"] == 50

    def test_sum_method(self) -> None:
        matrix, gene_sets, gene_names = self._make_scoring_data()
        result = gene_set_scoring(matrix, gene_sets, gene_names, method="sum")
        # Sum scores should be positive for cells with elevated pathway
        pathway_a_scores = result["scores"]["pathway_A"]
        avg_first_half = sum(pathway_a_scores[:25]) / 25
        avg_second_half = sum(pathway_a_scores[25:]) / 25
        assert avg_first_half > avg_second_half

    def test_gene_set_sizes(self) -> None:
        matrix, gene_sets, gene_names = self._make_scoring_data()
        result = gene_set_scoring(matrix, gene_sets, gene_names)
        assert result["gene_set_sizes"]["pathway_A"] == 5
        assert result["gene_set_sizes"]["pathway_B"] == 5

    def test_missing_genes_handled(self) -> None:
        matrix = [[1.0, 2.0]] * 10
        gene_sets = {"missing": ["NOT_A_GENE"]}
        result = gene_set_scoring(matrix, gene_sets, ["a", "b"])
        # Should return zeros for missing gene sets
        assert all(s == 0.0 for s in result["scores"]["missing"])

    def test_invalid_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid method"):
            gene_set_scoring([[1.0]], {"s": ["g"]}, ["g"], method="bad")

    def test_gene_names_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            gene_set_scoring([[1.0, 2.0]], {"s": ["g"]}, ["a"])

    def test_mean_method_separates_elevated_cells(self) -> None:
        """Mean scores must be higher for cells with elevated pathway genes."""
        matrix, gene_sets, gene_names = self._make_scoring_data()
        result = gene_set_scoring(matrix, gene_sets, gene_names, method="mean", seed=0)
        scores = result["scores"]["pathway_A"]
        assert np.mean(scores[:25]) > np.mean(scores[25:])


# ---------------------------------------------------------------------------
# Tie-aware Wilcoxon rank-sum
# ---------------------------------------------------------------------------


class TestWilcoxonTieCorrection:
    """Single-cell counts are massively tied; the test must be tie-aware."""

    def test_helper_matches_scipy_mannwhitneyu_with_ties(self) -> None:
        from scipy.stats import mannwhitneyu

        from metainformant.singlecell.differential.expression import _wilcoxon_rank_sum

        a = [5.0, 5.0, 5.0, 1.0, 2.0, 0.0]
        b = [5.0, 5.0, 0.0, 0.0, 0.0, 3.0]
        expected = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
        assert expected < 1.0
        assert _wilcoxon_rank_sum(a, b) == pytest.approx(expected)

    def test_public_wilcoxon_matches_scipy_on_tied_counts(self) -> None:
        from scipy.stats import mannwhitneyu

        group0 = [3.0, 3.0, 3.0, 1.0, 2.0, 0.0]
        group1 = [3.0, 3.0, 0.0, 0.0, 0.0, 4.0]
        matrix = [[v] for v in group0 + group1]
        groups = [0] * 6 + [1] * 6
        results = differential_expression(matrix, groups, ["g_tied"], method="wilcoxon")
        expected = float(mannwhitneyu(group0, group1, alternative="two-sided").pvalue)
        assert results[0]["p_value"] == pytest.approx(expected)


# ---------------------------------------------------------------------------
# BH family independence of the effect-size filter
# ---------------------------------------------------------------------------


class TestBHFamilyIndependenceOfEffectSizeFilter:
    """min_log2fc filters reported results, never the BH multiple-testing family."""

    @staticmethod
    def _four_gene_matrix() -> tuple[list[list[float]], list[int], list[str]]:
        # 8 cells (rows) x 4 genes (columns), 4 cells per group. g_flat has
        # near-equal means so it is removed by min_log2fc=2.0 but must stay
        # inside the BH family.
        matrix = [
            [10.0, 10.0, 10.0, 40.0],  # cell 0 (group 0)
            [10.0, 12.0, 22.0, 42.0],
            [10.0, 11.0, 21.0, 41.0],
            [11.0, 9.0, 19.0, 39.0],
            [10.0, 50.0, 100.0, 8.0],  # cell 4 (group 1)
            [10.0, 52.0, 102.0, 10.0],
            [11.0, 51.0, 101.0, 9.0],
            [10.0, 49.0, 99.0, 7.0],
        ]
        groups = [0] * 4 + [1] * 4
        names = ["g_flat", "g_up", "g_up2", "g_down"]
        return matrix, groups, names

    def test_filtered_out_gene_stays_in_bh_family(self) -> None:
        matrix, groups, names = self._four_gene_matrix()
        all_res = differential_expression(
            matrix, groups, names, method="wilcoxon", min_log2fc=0.0
        )
        filtered = differential_expression(
            matrix, groups, names, method="wilcoxon", min_log2fc=2.0
        )

        assert {r["gene"] for r in all_res} == set(names)
        assert "g_flat" not in {r["gene"] for r in filtered}

        by_gene = {r["gene"]: r["adjusted_p"] for r in all_res}
        for r in filtered:
            assert r["adjusted_p"] == by_gene[r["gene"]]

    def test_bh_multiplier_reflects_full_family(self) -> None:
        matrix, groups, names = self._four_gene_matrix()
        filtered = differential_expression(
            matrix, groups, names, method="wilcoxon", min_log2fc=2.0
        )
        p_up = next(r["p_value"] for r in filtered if r["gene"] == "g_up")
        adj_up = next(r["adjusted_p"] for r in filtered if r["gene"] == "g_up")
        # The three DE genes tie on p, so BH's step-up cascade settles all of
        # them at the rank-3 multiplier: p * 4 / 3 over the FULL 4-gene family.
        # A 3-gene family (filter applied pre-BH) would instead give p * 3 / 3 = p.
        assert adj_up == pytest.approx(p_up * 4 / 3, rel=1e-12)
        assert adj_up > p_up


# ---------------------------------------------------------------------------
# Pseudobulk group consistency and normalization
# ---------------------------------------------------------------------------


class TestPseudobulkGroupConsistency:
    """A sample whose cells span several groups has no valid pseudobulk identity."""

    def test_sample_spanning_multiple_groups_raises(self) -> None:
        from metainformant.core.utils.errors import ValidationError

        matrix = [[5.0, 3.0]] * 12
        cell_labels = ["T_cell"] * 12
        sample_labels = ["S1"] * 6 + ["S2"] * 6
        # S2 carries a stray group-0 cell; both samples pass min_cells_per_sample
        groups = [0] * 6 + [1, 1, 1, 0, 1, 1]
        with pytest.raises(ValidationError, match="multiple groups"):
            pseudobulk_de(
                matrix, cell_labels, sample_labels, groups, gene_names=["g0", "g1"]
            )

    def test_consistent_labels_do_not_raise(self) -> None:
        matrix = [[5.0, 3.0]] * 12
        cell_labels = ["T_cell"] * 12
        sample_labels = ["S1"] * 6 + ["S2"] * 6
        groups = [0] * 6 + [1] * 6
        results = pseudobulk_de(
            matrix, cell_labels, sample_labels, groups, gene_names=["g0", "g1"]
        )
        assert len(results) == 2


class TestPseudobulkCpmNormalization:
    """Pseudobulk sums are CPM-normalized per sample before testing."""

    @staticmethod
    def _depth_design() -> tuple[
        list[list[float]], list[str], list[str], list[int], list[str]
    ]:
        # 4 samples x 5 cells; group-1 samples are sequenced 2x deeper but
        # carry the same per-gene proportions as their group-0 counterparts.
        # Proportions differ slightly BETWEEN the two samples of a group so
        # the Welch t-test keeps nonzero within-group variance after CPM
        # (identical per-sample CPM vectors would give a degenerate test).
        proportions = {
            "S1": [0.50, 0.25, 0.25],
            "S2": [0.48, 0.26, 0.26],
            "S3": [0.50, 0.25, 0.25],
            "S4": [0.48, 0.26, 0.26],
        }
        totals = {"S1": 100.0, "S2": 120.0, "S3": 200.0, "S4": 240.0}
        matrix: list[list[float]] = []
        sample_labels: list[str] = []
        cell_labels: list[str] = []
        groups: list[int] = []
        for sample, group in (("S1", 0), ("S2", 0), ("S3", 1), ("S4", 1)):
            for _ in range(5):
                total = totals[sample]
                matrix.append([total * p for p in proportions[sample]])
                sample_labels.append(sample)
                cell_labels.append("T_cell")
                groups.append(group)
        return matrix, cell_labels, sample_labels, groups, ["g0", "g1", "g2"]

    def test_cpm_removes_library_depth_effect(self) -> None:
        matrix, cell_labels, sample_labels, groups, gene_names = self._depth_design()
        results = pseudobulk_de(
            matrix, cell_labels, sample_labels, groups, gene_names=gene_names
        )
        assert len(results) == 3
        for r in results:
            assert r["normalization"] == "cpm"
            # Same per-gene proportions in both groups: raw per-sample sums
            # (500/600 vs 1000/1200 total) would be called significant, but
            # CPM values are identical in both groups, so p is exactly 1.
            assert r["p_value"] == pytest.approx(1.0)
            assert r["log2fc"] == pytest.approx(0.0, abs=1e-9)
            assert r["mean_group1"] == pytest.approx(r["mean_group2"], rel=1e-9)

    def test_cpm_log2fc_tracks_proportion_change(self) -> None:
        import math

        matrix, cell_labels, sample_labels, groups, gene_names = self._depth_design()
        for i, sample in enumerate(sample_labels):
            if sample == "S3":
                matrix[i] = [150.0, 25.0, 25.0]  # 75% of the cell total, vs 50%
            elif sample == "S4":
                matrix[i] = [172.8, 33.6, 33.6]  # 72% / 14% / 14%
        results = pseudobulk_de(
            matrix, cell_labels, sample_labels, groups, gene_names=gene_names
        )
        g0 = next(r for r in results if r["gene"] == "g0")
        # CPM means: group0 (500000+480000)/2 = 490000; group1 (750000+720000)/2
        # = 735000. compute_log_fold_change takes (group0, group1) and adds a
        # pseudocount of 1, so the sign is negative for an up-in-group1 gene.
        assert g0["log2fc"] == pytest.approx(
            math.log2((490000 + 1) / (735000 + 1)), rel=1e-12
        )
        assert g0["p_value"] < 0.2
        g1 = next(r for r in results if r["gene"] == "g1")
        # CPM means: group0 (250000+260000)/2 = 255000; group1 (125000+140000)/2
        # = 132500.
        assert g1["log2fc"] == pytest.approx(
            math.log2((255000 + 1) / (132500 + 1)), rel=1e-12
        )
        assert g1["p_value"] < 0.2
