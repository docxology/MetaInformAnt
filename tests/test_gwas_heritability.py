"""Tests for GWAS SNP heritability estimation module."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.analysis.heritability import (
    estimate_heritability,
    heritability_bar_chart,
    partition_heritability_by_chromosome,
)


class TestEstimateHeritability:
    """Tests for the estimate_heritability function."""

    def test_basic_heritability_estimation(self) -> None:
        """Estimate h2 from synthetic kinship with known genetic component."""
        np.random.seed(42)
        n = 20

        # Build a kinship matrix with population structure
        A = np.random.randn(n, 5)
        K = A @ A.T / 5
        # Normalize to have unit diagonal on average
        diag_mean = np.mean(np.diag(K))
        K = K / diag_mean

        # Generate phenotypes with a genetic component
        sigma_g_true = 1.0
        sigma_e_true = 1.0
        L = np.linalg.cholesky(K + 1e-6 * np.eye(n))
        g = L @ np.random.randn(n) * sigma_g_true
        e = np.random.randn(n) * sigma_e_true
        y = list(g + e + 5.0)

        result = estimate_heritability(K, y)
        assert result["status"] == "success"
        assert 0.0 <= result["h2"] <= 1.0
        assert result["h2_se"] > 0
        assert result["sigma_g"] >= 0
        assert result["sigma_e"] >= 0
        assert result["n_samples"] == n
        assert result["method"] == "reml"
        assert "log_likelihood" in result

    def test_h2_bounded_zero_one(self) -> None:
        """h2 estimate should always be in [0, 1] regardless of inputs."""
        np.random.seed(99)
        n = 20
        K = np.eye(n) + 0.2 * np.ones((n, n))
        phenotypes = list(np.random.randn(n) * 3.0 + 10.0)

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "success"
        assert 0.0 <= result["h2"] <= 1.0

    def test_zero_variance_phenotype(self) -> None:
        """Zero-variance phenotype should return error gracefully."""
        n = 20
        K = np.eye(n)
        phenotypes = [5.0] * n  # constant phenotype

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "error"
        assert "zero" in result["message"].lower() or "variance" in result["message"].lower()

    def test_too_few_samples(self) -> None:
        """Fewer than 3 samples should return an error."""
        K = np.eye(2)
        phenotypes = [1.0, 2.0]

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "error"
        assert "3" in result["message"]

    def test_minimum_samples(self) -> None:
        """3 samples should work (minimum viable)."""
        np.random.seed(7)
        n = 3
        K = np.eye(n) + 0.1 * np.ones((n, n))
        phenotypes = [1.0, 5.0, 3.0]

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "success"
        assert 0.0 <= result["h2"] <= 1.0

    def test_identity_kinship_no_relatedness(self) -> None:
        """With identity kinship (no relatedness), h2 should be estimable."""
        np.random.seed(12)
        n = 20
        K = np.eye(n)
        phenotypes = list(np.random.randn(n) * 2.0)

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "success"
        assert 0.0 <= result["h2"] <= 1.0

    def test_mismatched_kinship_shape(self) -> None:
        """Kinship matrix that does not match sample count should error."""
        K = np.eye(5)
        phenotypes = [1.0, 2.0, 3.0]

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "error"
        assert "shape" in result["message"].lower() or "match" in result["message"].lower()

    def test_return_dict_structure(self) -> None:
        """Verify complete return dictionary structure on success."""
        np.random.seed(55)
        n = 15
        K = np.eye(n) + 0.05 * np.ones((n, n))
        phenotypes = list(np.random.randn(n))

        result = estimate_heritability(K, phenotypes)
        assert result["status"] == "success"
        required_keys = {"status", "h2", "h2_se", "sigma_g", "sigma_e", "log_likelihood", "n_samples", "method"}
        assert required_keys.issubset(result.keys())


class TestPartitionHeritability:
    """Tests for partition_heritability_by_chromosome."""

    def test_three_chromosomes(self) -> None:
        """Partition h2 across 3 chromosomes with distinct kinship matrices."""
        np.random.seed(42)
        n = 20

        kinship_matrices = {}
        for chrom in [1, 2, 3]:
            A = np.random.randn(n, 3)
            K_chr = A @ A.T / 3
            diag_mean = np.mean(np.diag(K_chr))
            K_chr = K_chr / max(diag_mean, 1e-6)
            kinship_matrices[chrom] = K_chr

        phenotypes = list(np.random.randn(n) * 2.0 + 5.0)

        result = partition_heritability_by_chromosome(kinship_matrices, phenotypes)
        assert result["status"] == "success"
        assert result["n_chromosomes"] == 3
        assert len(result["per_chromosome"]) == 3
        assert result["total_h2"] >= 0
        assert result["total_h2"] <= 1.0

        for chrom_key in ["1", "2", "3"]:
            assert chrom_key in result["per_chromosome"]
            assert "h2" in result["per_chromosome"][chrom_key]
            assert "h2_se" in result["per_chromosome"][chrom_key]

    def test_empty_kinship_matrices(self) -> None:
        """Empty kinship_matrices dict should return error."""
        result = partition_heritability_by_chromosome({}, [1.0, 2.0, 3.0])
        assert result["status"] == "error"

    def test_too_few_samples_partition(self) -> None:
        """Partition with fewer than 3 samples should error."""
        K = {1: np.eye(2)}
        result = partition_heritability_by_chromosome(K, [1.0, 2.0])
        assert result["status"] == "error"


class TestHeritabilityBarChart:
    """Tests for heritability_bar_chart visualization."""

    def test_bar_chart_saves_file(self, tmp_path: Path) -> None:
        """Bar chart should save a PNG file when output_file is given."""
        h2_data = {
            "per_chromosome": {
                "1": {"h2": 0.10, "h2_se": 0.03},
                "2": {"h2": 0.08, "h2_se": 0.02},
                "3": {"h2": 0.12, "h2_se": 0.04},
            },
            "total_h2": 0.30,
        }
        output_file = tmp_path / "heritability.png"

        result = heritability_bar_chart(h2_data, output_file=output_file)

        if result["status"] == "skipped":
            pytest.skip("matplotlib not available")

        assert result["status"] == "success"
        assert output_file.exists()
        assert output_file.stat().st_size > 0

    def test_bar_chart_no_output_file(self) -> None:
        """Bar chart without output_file should succeed without saving."""
        h2_data = {
            "per_chromosome": {
                "1": {"h2": 0.05, "h2_se": 0.02},
                "2": {"h2": 0.07, "h2_se": 0.03},
            },
            "total_h2": 0.12,
        }

        result = heritability_bar_chart(h2_data)

        if result["status"] == "skipped":
            pytest.skip("matplotlib not available")

        assert result["status"] == "success"
        assert result["output_path"] is None

    def test_bar_chart_empty_data(self) -> None:
        """Empty per_chromosome data should return failed status."""
        h2_data = {"per_chromosome": {}, "total_h2": 0.0}
        result = heritability_bar_chart(h2_data)
        assert result["status"] in ("failed", "skipped")

    def test_bar_chart_return_structure(self, tmp_path: Path) -> None:
        """Verify the return dict always has status and output_path."""
        h2_data = {
            "per_chromosome": {"1": {"h2": 0.15, "h2_se": 0.05}},
            "total_h2": 0.15,
        }
        output_file = tmp_path / "chart.png"
        result = heritability_bar_chart(h2_data, output_file=output_file)
        assert "status" in result
        assert "output_path" in result
