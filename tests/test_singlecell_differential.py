"""Tests for single-cell differential expression analysis.

Real implementation testing for DE statistical tests, pseudobulk DE,
fold-change computation, volcano plot data, and gene set scoring.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.singlecell.differential.expression import (
    compute_log_fold_change,
    differential_expression,
    gene_set_scoring,
    pseudobulk_de,
    volcano_data,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_de_data(
    n_cells: int = 80,
    n_genes: int = 100,
    n_de_genes: int = 10,
    seed: int = 42,
) -> tuple[list[list[float]], list[int], list[str]]:
    """Create expression matrix with known differentially expressed genes.

    Group 0 (first half) has baseline expression; group 1 (second half)
    has elevated expression in the first n_de_genes genes.
    """
    rng = np.random.RandomState(seed)
    half = n_cells // 2
    gene_names = [f"gene_{i}" for i in range(n_genes)]

    matrix = rng.exponential(1.0, size=(n_cells, n_genes))
    # Upregulate first n_de_genes in group 1
    matrix[half:, :n_de_genes] += rng.exponential(5.0, size=(half, n_de_genes))

    groups = [0] * half + [1] * half
    return matrix.tolist(), groups, gene_names


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
        assert len(overlap) >= 1, f"Expected DE genes {de_gene_names} in results; found {found_genes}"

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
        results = differential_expression(matrix, groups, gene_names, method="wilcoxon", min_log2fc=2.0)
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
        results = differential_expression(matrix.tolist(), groups, gene_names, min_cells=3)
        gene_set = {r["gene"] for r in results}
        # gene_0 has only 1 expressing cell per group (< 3), should be skipped
        assert "g0" not in gene_set


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
                matrix[start:end, :5] += rng.exponential(3.0, size=(cells_per_sample, 5))

        return (
            matrix.tolist(),
            cell_labels_list,
            sample_labels_list,
            group_labels_list,
            gene_names,
        )

    def test_pseudobulk_returns_results(self) -> None:
        matrix, cell_labels, sample_labels, groups, gene_names = self._make_pseudobulk_data()
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
        matrix, cell_labels, sample_labels, groups, gene_names = self._make_pseudobulk_data()
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
        matrix, cell_labels, sample_labels, groups, gene_names = self._make_pseudobulk_data()
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
        result = volcano_data(self._make_de_results(), fc_threshold=1.0, p_threshold=0.05)
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
        result = volcano_data(self._make_de_results(), fc_threshold=3.0, p_threshold=0.001)
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
