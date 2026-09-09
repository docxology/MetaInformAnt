"""Tests for synthetic benchmark dataset generators.

Tests cover generate_benchmark_dataset (all tasks + errors), synthetic GWAS
variants, synthetic RNA-seq expression, evaluate_benchmark metrics, and the
benchmark_suite comparison runner.

All tests use real implementations (real-implementation policy), seeded for
reproducibility, and no network.
"""

from __future__ import annotations

import pytest

from metainformant.simulation.benchmark.generators import (
    benchmark_suite,
    evaluate_benchmark,
    generate_benchmark_dataset,
    generate_synthetic_expression,
    generate_synthetic_variants,
)


class TestGenerateBenchmarkDataset:
    """Test general benchmark dataset generation."""

    @pytest.mark.parametrize("task", ["classification", "regression", "clustering", "de_genes"])
    def test_dataset_structure(self, task):
        """Test that every task returns the documented keys and shapes."""
        data = generate_benchmark_dataset(task=task, n_samples=20, n_features=40, seed=1)
        assert len(data["X"]) == 20
        assert all(len(row) == 40 for row in data["X"])
        assert len(data["y"]) == 20
        assert len(data["true_labels"]) == 20
        assert data["metadata"]["task"] == task
        assert data["metadata"]["seed"] == 1
        assert isinstance(data["task_description"], str)

    def test_classification_labels_balanced(self):
        """Test binary classification labels."""
        data = generate_benchmark_dataset(task="classification", n_samples=10, n_features=20, seed=3)
        assert sorted(set(data["y"])) == [0, 1]

    def test_unknown_task_raises(self):
        """Test that an unknown task raises ValueError."""
        with pytest.raises(ValueError, match="Unknown task"):
            generate_benchmark_dataset(task="spam", n_samples=10, n_features=20)

    def test_deterministic_with_seed(self):
        """Test reproducibility with a fixed seed."""
        data_a = generate_benchmark_dataset(task="classification", n_samples=10, n_features=20, seed=5)
        data_b = generate_benchmark_dataset(task="classification", n_samples=10, n_features=20, seed=5)
        assert data_a["X"] == data_b["X"]
        assert data_a["y"] == data_b["y"]

    def test_difficulty_changes_signal_scale(self):
        """Test that difficulty maps to different noise/signal metadata."""
        easy = generate_benchmark_dataset(task="classification", n_samples=10, n_features=20, difficulty="easy", seed=2)
        hard = generate_benchmark_dataset(task="classification", n_samples=10, n_features=20, difficulty="hard", seed=2)
        assert easy["metadata"]["signal_scale"] > hard["metadata"]["signal_scale"]
        assert easy["metadata"]["noise_scale"] < hard["metadata"]["noise_scale"]


class TestGenerateSyntheticVariants:
    """Test synthetic GWAS data generation."""

    def test_structure(self):
        """Test documented output keys and shapes."""
        result = generate_synthetic_variants(n_variants=30, n_causal=5, n_samples=10, seed=1)
        assert set(result) == {"genotypes", "phenotypes", "causal_variants", "true_effects", "heritability", "maf"}
        assert len(result["genotypes"]) == 10
        assert all(len(row) == 30 for row in result["genotypes"])
        assert all(g in (0, 1, 2) for row in result["genotypes"] for g in row)
        assert len(result["phenotypes"]) == 10
        assert len(result["causal_variants"]) == 5
        assert len(result["true_effects"]) == 5
        assert len(result["maf"]) == 30

    def test_explicit_effect_sizes_respected(self):
        """Test that provided effect sizes are used for causal variants."""
        result = generate_synthetic_variants(n_variants=10, n_causal=2, effect_sizes=[0.5, -0.5], n_samples=5, seed=1)
        assert result["true_effects"][:2] == [0.5, -0.5]

    def test_deterministic_with_seed(self):
        """Test reproducibility with a fixed seed."""
        result_a = generate_synthetic_variants(n_variants=15, n_causal=3, n_samples=6, seed=7)
        result_b = generate_synthetic_variants(n_variants=15, n_causal=3, n_samples=6, seed=7)
        assert result_a["genotypes"] == result_b["genotypes"]
        assert result_a["phenotypes"] == result_b["phenotypes"]
        assert result_a["causal_variants"] == result_b["causal_variants"]

    def test_heritability_with_zero_variance_genetic_values(self):
        """Test heritability edge case when genetic values are all zero."""
        result = generate_synthetic_variants(n_variants=5, n_causal=0, n_samples=6, seed=1)
        assert result["heritability"] == 0.0
        assert result["causal_variants"] == []


class TestGenerateSyntheticExpression:
    """Test synthetic RNA-seq count generation."""

    def test_structure(self):
        """Test documented output keys and shapes."""
        result = generate_synthetic_expression(n_genes=20, n_samples=8, n_de_genes=4, seed=1)
        assert set(result) == {"counts", "groups", "de_genes", "true_fold_changes"}
        assert len(result["counts"]) == 20
        assert all(len(row) == 8 for row in result["counts"])
        assert all(isinstance(c, int) and c >= 0 for row in result["counts"] for c in row)
        assert result["groups"] == [0, 0, 0, 0, 1, 1, 1, 1]
        assert len(result["de_genes"]) == 4
        assert len(result["true_fold_changes"]) == 4

    def test_deterministic_with_seed(self):
        """Test reproducibility with a fixed seed."""
        result_a = generate_synthetic_expression(n_genes=10, n_samples=4, n_de_genes=2, seed=9)
        result_b = generate_synthetic_expression(n_genes=10, n_samples=4, n_de_genes=2, seed=9)
        assert result_a["counts"] == result_b["counts"]
        assert result_a["de_genes"] == result_b["de_genes"]

    def test_fold_changes_override(self):
        """Test that user-supplied fold changes are used verbatim."""
        result = generate_synthetic_expression(n_genes=10, n_samples=4, n_de_genes=2, fold_changes=[1.0, -2.0], seed=1)
        assert result["true_fold_changes"] == [1.0, -2.0]


class TestEvaluateBenchmark:
    """Test prediction evaluation."""

    def test_perfect_classification(self):
        """Test that perfect predictions score 1.0 on all metrics."""
        truth = [0, 0, 1, 1]
        result = evaluate_benchmark(truth, truth, task="classification")
        assert result["metrics"]["accuracy"] == 1.0
        assert result["metrics"]["precision"] == 1.0
        assert result["metrics"]["recall"] == 1.0
        assert result["metrics"]["f1_score"] == 1.0
        assert result["confusion_matrix"] == [[2, 0], [0, 2]]

    def test_imperfect_classification(self):
        """Test confusion matrix counts for one misclassification."""
        result = evaluate_benchmark([0, 1, 1, 1], [0, 0, 1, 1], task="classification")
        assert result["metrics"]["accuracy"] == 0.75
        assert result["confusion_matrix"] == [[1, 1], [0, 2]]

    def test_perfect_regression(self):
        """Test that perfect predictions give R2 of 1.0 and RMSE of 0."""
        result = evaluate_benchmark([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], task="regression")
        assert result["metrics"]["mse"] == 0.0
        assert result["metrics"]["rmse"] == 0.0
        assert result["metrics"]["r_squared"] == 1.0

    def test_regression_metrics(self):
        """Test MSE/MAE computation on known residuals."""
        result = evaluate_benchmark([1.0, 3.0], [2.0, 2.0], task="regression")
        assert result["metrics"]["mae"] == 1.0
        assert result["metrics"]["mse"] == 1.0

    def test_clustering_ari_bounds(self):
        """Test ARI is 1.0 for identical partitions and <= 1 in general."""
        truth = [0, 0, 1, 1]
        result = evaluate_benchmark(truth, truth, task="clustering")
        assert result["metrics"]["adjusted_rand_index"] == pytest.approx(1.0)

    def test_length_mismatch_raises(self):
        """Test that mismatched lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            evaluate_benchmark([0, 1], [0], task="classification")

    def test_unknown_task_raises(self):
        """Test that an unknown task raises ValueError."""
        with pytest.raises(ValueError, match="Unknown task"):
            evaluate_benchmark([0], [0], task="spam")


class TestBenchmarkSuite:
    """Test the multi-method comparison runner."""

    def _dataset(self):
        """Build a small deterministic classification dataset."""
        return generate_benchmark_dataset(task="classification", n_samples=8, n_features=10, seed=1)

    def test_rankings_and_best_method(self):
        """Test that rankings sort by primary metric and pick the best method."""
        dataset = self._dataset()

        def perfect(X, y):
            return list(y)

        def constant_zero(X, y):
            return [0] * len(y)

        result = benchmark_suite({"perfect": perfect, "constant_zero": constant_zero}, dataset, n_repeats=2)
        assert set(result["results_per_method"]) == {"perfect", "constant_zero"}
        assert result["best_method"] == "perfect"
        assert result["rankings"][0][0] == "perfect"
        assert result["rankings"][0][1] >= result["rankings"][1][1]
        assert "Benchmark comparison" in result["summary"]
        for method_results in result["results_per_method"].values():
            assert len(method_results) == 2

    def test_failing_method_records_error(self):
        """Test that a method raising an exception yields a failed result entry."""
        dataset = self._dataset()

        def broken(X, y):
            raise RuntimeError("boom")

        result = benchmark_suite({"broken": broken}, dataset, n_repeats=1)
        entry = result["results_per_method"]["broken"][0]
        assert "Failed: boom" in entry["summary"]
        assert result["best_method"] == "broken"
