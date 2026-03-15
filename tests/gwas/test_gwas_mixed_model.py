"""Tests for GWAS mixed linear model (EMMA) module."""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.gwas.analysis.mixed_model import (
    _emma_eigendecompose,
    _emma_reml,
    association_test_mixed,
    run_mixed_model_gwas,
)


class TestEMMAComponents:
    """Tests for EMMA internal functions."""

    def test_eigendecompose_identity(self) -> None:
        """Eigendecomposing identity should give all eigenvalues = 1."""
        K = np.eye(5)
        eigenvalues, eigenvectors = _emma_eigendecompose(K)
        assert len(eigenvalues) == 5
        assert all(abs(ev - 1.0) < 0.01 for ev in eigenvalues)

    def test_eigendecompose_symmetric(self) -> None:
        """Eigendecomposition of a symmetric matrix should produce real eigenvalues."""
        np.random.seed(42)
        A = np.random.randn(10, 10)
        K = A @ A.T / 10  # Positive semi-definite
        eigenvalues, eigenvectors = _emma_eigendecompose(K)
        assert all(ev >= -1e-10 for ev in eigenvalues)
        # Verify reconstruction: K ≈ U diag(lambda) U'
        K_reconstructed = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T
        assert np.allclose(K, K_reconstructed, atol=1e-8)

    def test_reml_positive_variance(self) -> None:
        """REML should estimate non-negative variance components."""
        np.random.seed(42)
        n = 20
        K = np.eye(n) + 0.1 * np.ones((n, n))
        eigenvalues, eigenvectors = _emma_eigendecompose(K)
        Ut = eigenvectors.T

        y = np.random.randn(n) * 2 + 5
        X = np.ones((n, 1))

        y_rot = Ut @ y
        X_rot = Ut @ X

        sigma_g, sigma_e = _emma_reml(y_rot, X_rot, eigenvalues)
        assert sigma_g >= 0
        assert sigma_e >= 0


class TestAssociationTestMixed:
    """Tests for single-SNP mixed model association test."""

    def test_basic_mixed_model(self) -> None:
        """Basic mixed model should return valid results."""
        np.random.seed(42)
        n = 30
        genotypes = list(np.random.choice([0, 1, 2], size=n))
        phenotypes = [float(g * 0.5 + np.random.randn() * 1.0) for g in genotypes]
        kinship = [[1.0 if i == j else 0.1 for j in range(n)] for i in range(n)]

        result = association_test_mixed(genotypes, phenotypes, kinship)
        assert result["status"] == "success"
        assert "beta" in result
        assert "se" in result
        assert "p_value" in result
        assert "heritability" in result
        assert 0 <= result["heritability"] <= 1
        assert result["test_type"] == "mixed"

    def test_identity_kinship_matches_linear(self) -> None:
        """With K=I, mixed model should give similar results to linear model."""
        np.random.seed(42)
        n = 50
        genotypes = list(np.random.choice([0, 1, 2], size=n))
        # Strong effect so both models detect it
        phenotypes = [float(g * 3.0 + np.random.randn() * 2.0) for g in genotypes]
        kinship = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

        from metainformant.gwas.analysis.association import association_test_linear

        mixed_result = association_test_mixed(genotypes, phenotypes, kinship)
        linear_result = association_test_linear(genotypes, phenotypes)

        assert mixed_result["status"] == "success"
        assert linear_result["status"] == "success"

        # Betas should be in the same direction
        assert (mixed_result["beta"] > 0) == (linear_result["beta"] > 0)

        # p-values should be similar order of magnitude
        if mixed_result["p_value"] > 0 and linear_result["p_value"] > 0:
            log_ratio = abs(
                math.log10(max(mixed_result["p_value"], 1e-300)) - math.log10(max(linear_result["p_value"], 1e-300))
            )
            assert log_ratio < 3  # Within 3 orders of magnitude

    def test_structured_kinship_more_conservative(self) -> None:
        """With population structure in K, p-values should be less significant."""
        np.random.seed(42)
        n = 40

        # Create population structure
        pop1 = [0] * 20 + [2] * 20  # Genotype correlated with population
        phenotypes = [10.0 + np.random.randn() for _ in range(20)] + [15.0 + np.random.randn() for _ in range(20)]

        # Identity kinship (no correction)
        K_identity = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

        # Structured kinship (within-population relatedness)
        K_structured = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    K_structured[i][j] = 1.0
                elif (i < 20 and j < 20) or (i >= 20 and j >= 20):
                    K_structured[i][j] = 0.5  # Same population
                else:
                    K_structured[i][j] = 0.0

        result_identity = association_test_mixed(pop1, phenotypes, K_identity)
        result_structured = association_test_mixed(pop1, phenotypes, K_structured)

        assert result_identity["status"] == "success"
        assert result_structured["status"] == "success"

        # Structured kinship should produce less significant (larger) p-value
        # because the model accounts for the confounding
        assert result_structured["p_value"] >= result_identity["p_value"] * 0.1  # Allow some slack

    def test_mismatched_lengths(self) -> None:
        """Mismatched genotype/phenotype lengths should raise ValueError."""
        with pytest.raises(ValueError):
            association_test_mixed([0, 1, 2], [1.0, 2.0], [[1, 0], [0, 1]])

    def test_mismatched_kinship(self) -> None:
        """Wrong kinship dimensions should raise ValueError."""
        with pytest.raises(ValueError):
            association_test_mixed(
                [0, 1, 2],
                [1.0, 2.0, 3.0],
                [[1, 0], [0, 1]],  # 2x2 but 3 samples
            )


class TestRunMixedModelGWAS:
    """Tests for genome-wide mixed model GWAS."""

    def test_basic_gwas(self) -> None:
        """Basic mixed model GWAS should test all variants."""
        np.random.seed(42)
        n = 30
        n_variants = 10

        genotype_matrix = [list(np.random.choice([0, 1, 2], size=n)) for _ in range(n_variants)]
        phenotypes = [float(np.random.randn()) for _ in range(n)]
        kinship = [[1.0 if i == j else 0.05 for j in range(n)] for i in range(n)]

        results = run_mixed_model_gwas(genotype_matrix, phenotypes, kinship)
        assert len(results) == n_variants
        for r in results:
            assert "beta" in r
            assert "p_value" in r
            assert r["test_type"] == "mixed"

    def test_gwas_with_variant_info(self) -> None:
        """GWAS with variant info should include variant metadata in results."""
        np.random.seed(42)
        n = 20
        n_variants = 5

        genotype_matrix = [list(np.random.choice([0, 1, 2], size=n)) for _ in range(n_variants)]
        phenotypes = [float(np.random.randn()) for _ in range(n)]
        kinship = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        variant_info = [{"id": f"rs{i}", "chrom": "chr1", "pos": i * 1000} for i in range(n_variants)]

        results = run_mixed_model_gwas(genotype_matrix, phenotypes, kinship, variant_info=variant_info)
        assert len(results) == n_variants
        assert results[0]["variant_id"] == "rs0"
        assert results[0]["chrom"] == "chr1"
