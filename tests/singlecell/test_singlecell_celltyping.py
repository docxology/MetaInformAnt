"""Tests for single-cell cell type annotation and classification.

Real implementation testing for marker-based annotation, label transfer,
novel cell type detection, and cell type composition analysis.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest

from metainformant.singlecell.celltyping.annotation import (
    annotate_by_markers,
    cell_type_composition,
    find_novel_types,
    score_cell_type,
    transfer_labels,
)

# ---------------------------------------------------------------------------
# Fixtures: realistic single-cell data generators
# ---------------------------------------------------------------------------


def _make_expression_matrix(
    n_cells: int = 80,
    n_genes: int = 200,
    n_types: int = 3,
    seed: int = 42,
) -> tuple[list[list[float]], list[str], dict[str, list[str]]]:
    """Build a synthetic expression matrix with known marker structure.

    Returns matrix (cells x genes), gene names, and marker dict.
    Each cell type has 5 marker genes with elevated expression, making
    the ground truth recoverable by overlap / scoring methods.
    """
    rng = np.random.RandomState(seed)
    gene_names = [f"gene_{i}" for i in range(n_genes)]

    cells_per_type = n_cells // n_types
    marker_genes: dict[str, list[str]] = {}
    type_names = [f"type_{t}" for t in range(n_types)]
    matrix = rng.exponential(0.5, size=(n_cells, n_genes))

    for t_idx, t_name in enumerate(type_names):
        start_gene = t_idx * 5
        markers = gene_names[start_gene : start_gene + 5]
        marker_genes[t_name] = markers

        start_cell = t_idx * cells_per_type
        end_cell = start_cell + cells_per_type
        for g in range(start_gene, start_gene + 5):
            matrix[start_cell:end_cell, g] += rng.exponential(5.0, size=cells_per_type)

    # Apply dropout
    dropout = rng.random(matrix.shape) < 0.2
    matrix[dropout] = 0.0

    return matrix.tolist(), gene_names, marker_genes


# ---------------------------------------------------------------------------
# annotate_by_markers
# ---------------------------------------------------------------------------


class TestAnnotateByMarkers:
    """Tests for annotate_by_markers function."""

    def test_overlap_method_returns_expected_keys(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix()
        result = annotate_by_markers(matrix, marker_genes, gene_names=gene_names, method="overlap")
        assert "cell_labels" in result
        assert "confidence_scores" in result
        assert "ambiguous_cells" in result
        assert "marker_stats" in result

    def test_overlap_labels_have_correct_length(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix(n_cells=60)
        result = annotate_by_markers(matrix, marker_genes, gene_names=gene_names, method="overlap")
        assert len(result["cell_labels"]) == 60
        assert len(result["confidence_scores"]) == 60

    def test_scoring_method_assigns_known_types(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix(seed=99)
        result = annotate_by_markers(matrix, marker_genes, gene_names=gene_names, method="scoring")
        labels = result["cell_labels"]
        # At least some cells should be assigned to each type
        unique_labels = set(labels)
        for t_name in marker_genes:
            assert t_name in unique_labels, f"Expected {t_name} to appear in labels"

    def test_correlation_method_runs(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix(n_cells=40)
        result = annotate_by_markers(matrix, marker_genes, gene_names=gene_names, method="correlation")
        assert len(result["cell_labels"]) == 40

    def test_threshold_produces_unassigned(self) -> None:
        """High threshold should leave some cells unassigned."""
        matrix, gene_names, marker_genes = _make_expression_matrix()
        result = annotate_by_markers(
            matrix,
            marker_genes,
            gene_names=gene_names,
            method="scoring",
            threshold=999.0,
        )
        assert "unassigned" in result["cell_labels"]

    def test_marker_stats_coverage(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix()
        result = annotate_by_markers(matrix, marker_genes, gene_names=gene_names, method="overlap")
        for t_name in marker_genes:
            stats = result["marker_stats"][t_name]
            assert "mean_score" in stats
            assert "n_cells_assigned" in stats
            assert "marker_coverage" in stats
            assert stats["marker_coverage"] == 1.0  # all markers present

    def test_invalid_method_raises(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix()
        with pytest.raises(ValueError, match="Invalid method"):
            annotate_by_markers(
                matrix,
                marker_genes,
                gene_names=gene_names,
                method="bogus",
            )

    def test_empty_matrix_raises(self) -> None:
        with pytest.raises(ValueError, match="at least one cell"):
            annotate_by_markers([], {"A": ["g"]}, gene_names=[], method="overlap")

    def test_no_marker_genes_found_raises(self) -> None:
        matrix = [[1.0, 2.0], [3.0, 4.0]]
        gene_names = ["x", "y"]
        with pytest.raises(ValueError, match="None of the provided marker"):
            annotate_by_markers(
                matrix,
                {"T": ["DOES_NOT_EXIST"]},
                gene_names=gene_names,
                method="overlap",
            )

    def test_gene_names_length_mismatch_raises(self) -> None:
        matrix = [[1.0, 2.0]]
        with pytest.raises(ValueError, match="does not match"):
            annotate_by_markers(
                matrix,
                {"A": ["g"]},
                gene_names=["a", "b", "c"],
                method="overlap",
            )

    def test_ambiguous_cells_detected(self) -> None:
        """When two types have similar marker expression, cells are flagged."""
        n_cells = 40
        n_genes = 20
        rng = np.random.RandomState(7)
        matrix = rng.exponential(2.0, size=(n_cells, n_genes)).tolist()
        gene_names = [f"g{i}" for i in range(n_genes)]
        # Markers overlap significantly
        markers = {"A": gene_names[:5], "B": gene_names[:5]}
        result = annotate_by_markers(matrix, markers, gene_names=gene_names, method="overlap")
        # All cells should be ambiguous because markers are identical
        assert len(result["ambiguous_cells"]) == n_cells

    def test_numpy_array_input(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix(n_cells=30)
        arr = np.array(matrix)
        result = annotate_by_markers(arr, marker_genes, gene_names=gene_names, method="scoring")
        assert len(result["cell_labels"]) == 30

    def test_auto_generated_gene_names(self) -> None:
        """When gene_names is None, automatic names are generated."""
        matrix, _, marker_genes_raw = _make_expression_matrix(n_cells=30)
        # Remap markers to auto-generated names
        marker_genes = {"A": ["gene_0", "gene_1"]}
        result = annotate_by_markers(matrix, marker_genes, method="scoring")
        assert len(result["cell_labels"]) == 30


# ---------------------------------------------------------------------------
# score_cell_type
# ---------------------------------------------------------------------------


class TestScoreCellType:
    """Tests for score_cell_type function."""

    def test_positive_score_for_high_markers(self) -> None:
        rng = random.Random(0)
        gene_names = [f"g{i}" for i in range(100)]
        expression = [rng.uniform(0, 1) for _ in range(100)]
        # Set marker genes to high expression
        for idx in range(5):
            expression[idx] = 10.0
        markers = gene_names[:5]
        score = score_cell_type(expression, gene_names, markers, seed=42)
        assert score > 0

    def test_zero_score_when_markers_absent(self) -> None:
        gene_names = [f"g{i}" for i in range(50)]
        expression = [1.0] * 50
        score = score_cell_type(expression, gene_names, ["MISSING_1", "MISSING_2"], seed=0)
        assert score == 0.0

    def test_reproducible_with_seed(self) -> None:
        gene_names = [f"g{i}" for i in range(80)]
        expression = [float(i) for i in range(80)]
        markers = gene_names[:3]
        s1 = score_cell_type(expression, gene_names, markers, seed=123)
        s2 = score_cell_type(expression, gene_names, markers, seed=123)
        assert s1 == s2

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            score_cell_type([1.0, 2.0], ["a", "b", "c"], ["a"])


# ---------------------------------------------------------------------------
# transfer_labels
# ---------------------------------------------------------------------------


class TestTransferLabels:
    """Tests for transfer_labels (kNN label transfer)."""

    def _make_ref_query(self, seed: int = 42) -> tuple[list[list[float]], list[str], list[list[float]]]:
        rng = np.random.RandomState(seed)
        n_ref = 60
        n_query = 20
        n_genes = 50

        ref_data = rng.randn(n_ref, n_genes)
        ref_labels: list[str] = []
        for i in range(n_ref):
            cluster = i % 3
            ref_data[i, cluster * 10 : cluster * 10 + 10] += 4.0
            ref_labels.append(f"type_{cluster}")

        query_data = rng.randn(n_query, n_genes)
        for i in range(n_query):
            cluster = i % 3
            query_data[i, cluster * 10 : cluster * 10 + 10] += 4.0

        return ref_data.tolist(), ref_labels, query_data.tolist()

    def test_returns_expected_keys(self) -> None:
        ref, labels, query = self._make_ref_query()
        result = transfer_labels(ref, labels, query, n_neighbors=5)
        assert "predicted_labels" in result
        assert "prediction_scores" in result
        assert "mapping_quality" in result

    def test_predicted_labels_length(self) -> None:
        ref, labels, query = self._make_ref_query()
        result = transfer_labels(ref, labels, query, n_neighbors=5)
        assert len(result["predicted_labels"]) == 20

    def test_prediction_scores_range(self) -> None:
        ref, labels, query = self._make_ref_query()
        result = transfer_labels(ref, labels, query, n_neighbors=10)
        for s in result["prediction_scores"]:
            assert 0.0 <= s <= 1.0

    def test_mapping_quality_stats(self) -> None:
        ref, labels, query = self._make_ref_query()
        result = transfer_labels(ref, labels, query, n_neighbors=5)
        mq = result["mapping_quality"]
        assert "mean_confidence" in mq
        assert "n_high_confidence" in mq
        assert "n_query_cells" in mq
        assert mq["n_query_cells"] == 20

    def test_dimension_mismatch_raises(self) -> None:
        ref = [[1.0, 2.0], [3.0, 4.0]]
        query = [[1.0, 2.0, 3.0]]
        with pytest.raises(ValueError, match="same number of genes"):
            transfer_labels(ref, ["a", "b"], query)

    def test_label_length_mismatch_raises(self) -> None:
        ref = [[1.0, 2.0], [3.0, 4.0]]
        query = [[1.0, 2.0]]
        with pytest.raises(ValueError, match="must match"):
            transfer_labels(ref, ["only_one"], query)

    def test_empty_data_raises(self) -> None:
        with pytest.raises(ValueError, match="at least one cell"):
            transfer_labels([], [], [[1.0]])

    def test_numpy_array_input(self) -> None:
        ref, labels, query = self._make_ref_query()
        result = transfer_labels(np.array(ref), labels, np.array(query), n_neighbors=5)
        assert len(result["predicted_labels"]) == 20


# ---------------------------------------------------------------------------
# find_novel_types
# ---------------------------------------------------------------------------


class TestFindNovelTypes:
    """Tests for find_novel_types function."""

    def test_unassigned_labels_flagged(self) -> None:
        matrix = [[1.0, 2.0]] * 10
        labels = ["T_cell"] * 7 + ["unassigned"] * 3
        result = find_novel_types(matrix, labels)
        assert result["n_novel"] == 3
        assert len(result["novel_cell_indices"]) == 3

    def test_no_novel_when_all_assigned(self) -> None:
        matrix = [[1.0, 2.0]] * 5
        labels = ["A", "B", "A", "B", "A"]
        result = find_novel_types(matrix, labels)
        assert result["n_novel"] == 0

    def test_rescoring_with_markers(self) -> None:
        matrix, gene_names, marker_genes = _make_expression_matrix(n_cells=30)
        labels = ["type_0"] * 10 + ["type_1"] * 10 + ["type_2"] * 10
        result = find_novel_types(
            matrix,
            labels,
            marker_genes=marker_genes,
            gene_names=gene_names,
            threshold=100.0,  # very high: most cells become novel
        )
        assert result["n_novel"] > 0
        assert result["fraction_novel"] > 0.0

    def test_cluster_summary_present(self) -> None:
        matrix = np.random.RandomState(0).randn(20, 10).tolist()
        labels = ["assigned"] * 10 + ["unknown"] * 10
        result = find_novel_types(matrix, labels)
        assert "cluster_summary" in result
        assert len(result["cluster_summary"]) >= 1

    def test_label_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            find_novel_types([[1.0]], ["a", "b"])


# ---------------------------------------------------------------------------
# cell_type_composition
# ---------------------------------------------------------------------------


class TestCellTypeComposition:
    """Tests for cell_type_composition function."""

    def test_overall_proportions_sum_to_one(self) -> None:
        labels = ["A"] * 30 + ["B"] * 20 + ["C"] * 50
        result = cell_type_composition(labels)
        total = sum(result["overall"].values())
        assert abs(total - 1.0) < 1e-9

    def test_overall_proportions_correct(self) -> None:
        labels = ["A"] * 40 + ["B"] * 60
        result = cell_type_composition(labels)
        assert abs(result["overall"]["A"] - 0.4) < 1e-9
        assert abs(result["overall"]["B"] - 0.6) < 1e-9

    def test_per_group_with_groups(self) -> None:
        labels = ["A", "B", "A", "B"]
        groups = ["g1", "g1", "g2", "g2"]
        result = cell_type_composition(labels, groups=groups)
        assert "per_group" in result
        assert result["n_groups"] == 2
        assert abs(result["per_group"]["g1"]["A"] - 0.5) < 1e-9

    def test_confidence_intervals_present(self) -> None:
        labels = ["X"] * 50 + ["Y"] * 50
        result = cell_type_composition(labels)
        assert "confidence_intervals" in result
        for key, (lo, hi) in result["confidence_intervals"].items():
            assert 0.0 <= lo <= hi <= 1.0

    def test_confidence_intervals_with_groups(self) -> None:
        labels = ["A", "B"] * 20
        groups = ["g1"] * 20 + ["g2"] * 20
        result = cell_type_composition(labels, groups=groups)
        ci = result["confidence_intervals"]
        assert any("g1::" in k for k in ci)
        assert any("g2::" in k for k in ci)

    def test_empty_labels_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            cell_type_composition([])

    def test_groups_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            cell_type_composition(["A", "B"], groups=["g1"])

    def test_n_cells_and_n_types(self) -> None:
        labels = ["A"] * 10 + ["B"] * 10 + ["C"] * 10
        result = cell_type_composition(labels)
        assert result["n_cells"] == 30
        assert result["n_types"] == 3

    def test_different_confidence_levels(self) -> None:
        labels = ["A"] * 50 + ["B"] * 50
        r90 = cell_type_composition(labels, confidence_level=0.90)
        r99 = cell_type_composition(labels, confidence_level=0.99)
        # 99% interval should be wider than 90%
        ci_90 = r90["confidence_intervals"]
        ci_99 = r99["confidence_intervals"]
        for key in ci_90:
            lo_90, hi_90 = ci_90[key]
            lo_99, hi_99 = ci_99[key]
            width_90 = hi_90 - lo_90
            width_99 = hi_99 - lo_99
            assert width_99 >= width_90 - 1e-9
