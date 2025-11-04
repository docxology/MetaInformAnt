"""Tests for kinship matrix mathematical properties."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.gwas import compute_kinship_matrix


class TestKinshipMatrixProperties:
    """Test mathematical properties of kinship matrices."""

    def test_kinship_matrix_symmetric(self):
        """Test that kinship matrix is symmetric."""
        genotypes = [
            [0, 1, 2, 0, 1],
            [0, 1, 2, 0, 1],
            [2, 1, 0, 2, 1],
            [1, 1, 1, 1, 1],
        ]

        for method in ["vanraden", "astle", "yang"]:
            result = compute_kinship_matrix(genotypes, method=method)
            assert result["status"] == "success"
            
            kinship = np.array(result["kinship_matrix"])
            
            # Check symmetry: K = K^T
            assert np.allclose(kinship, kinship.T), f"Kinship matrix not symmetric for {method} method"

    def test_kinship_matrix_positive_semi_definite(self):
        """Test that kinship matrix is positive semi-definite."""
        genotypes = [
            [0, 1, 2, 0, 1, 2],
            [0, 1, 2, 0, 1, 2],
            [2, 1, 0, 2, 1, 0],
            [1, 1, 1, 1, 1, 1],
            [2, 0, 1, 2, 0, 1],
        ]

        for method in ["vanraden", "astle", "yang"]:
            result = compute_kinship_matrix(genotypes, method=method)
            assert result["status"] == "success"
            
            kinship = np.array(result["kinship_matrix"])
            
            # Check positive semi-definite: all eigenvalues >= 0
            eigenvalues = np.linalg.eigvals(kinship)
            assert np.all(eigenvalues >= -1e-10), (
                f"Kinship matrix not positive semi-definite for {method} method. "
                f"Min eigenvalue: {np.min(eigenvalues)}"
            )

    def test_kinship_matrix_diagonal_range(self):
        """Test that diagonal elements (self-kinship) are in reasonable range."""
        genotypes = [
            [0, 1, 2, 0, 1],
            [0, 1, 2, 0, 1],
            [2, 1, 0, 2, 1],
        ]

        for method in ["vanraden", "astle", "yang"]:
            result = compute_kinship_matrix(genotypes, method=method)
            assert result["status"] == "success"
            
            kinship = np.array(result["kinship_matrix"])
            diagonal = np.diag(kinship)
            
            # Diagonal should be non-negative (self-kinship)
            assert np.all(diagonal >= -1e-10), f"Negative diagonal values for {method} method"
            
            # Diagonal should typically be <= 1 (though not strictly required)
            # For most methods, self-kinship should be reasonable

    def test_kinship_matrix_identical_samples(self):
        """Test that identical samples have high kinship."""
        genotypes = [
            [0, 1, 2, 0, 1],
            [0, 1, 2, 0, 1],  # Identical to first
            [2, 1, 0, 2, 1],  # Different
        ]

        for method in ["vanraden", "astle", "yang"]:
            result = compute_kinship_matrix(genotypes, method=method)
            assert result["status"] == "success"
            
            kinship = np.array(result["kinship_matrix"])
            
            # Samples 0 and 1 should have high kinship (they're identical)
            kinship_0_1 = kinship[0, 1]
            
            # For identical samples, kinship should be high
            # Exact value depends on method, but should be > 0.5
            assert kinship_0_1 > 0.5, (
                f"Identical samples have low kinship ({kinship_0_1:.4f}) for {method} method"
            )

    def test_kinship_method_comparison(self):
        """Test that different methods give consistent but not identical results."""
        genotypes = [
            [0, 1, 2, 0, 1, 2],
            [0, 1, 2, 0, 1, 2],
            [2, 1, 0, 2, 1, 0],
            [1, 1, 1, 1, 1, 1],
        ]

        results = {}
        for method in ["vanraden", "astle", "yang"]:
            result = compute_kinship_matrix(genotypes, method=method)
            assert result["status"] == "success"
            results[method] = np.array(result["kinship_matrix"])

        # Methods should give similar but not identical results
        # Check that all methods identify same pairs as related
        vanraden = results["vanraden"]
        astle = results["astle"]
        yang = results["yang"]

        # High kinship pairs should be consistent across methods
        # (samples 0 and 1 are identical)
        assert vanraden[0, 1] > 0.5
        assert astle[0, 1] > 0.5
        assert yang[0, 1] > 0.5

        # Methods may differ in exact values but should agree on relative relationships
        # Check that correlation between methods is positive
        vanraden_flat = vanraden[np.triu_indices_from(vanraden, k=1)]
        astle_flat = astle[np.triu_indices_from(astle, k=1)]
        yang_flat = yang[np.triu_indices_from(yang, k=1)]

        corr_va = np.corrcoef(vanraden_flat, astle_flat)[0, 1]
        corr_vy = np.corrcoef(vanraden_flat, yang_flat)[0, 1]
        corr_ay = np.corrcoef(astle_flat, yang_flat)[0, 1]

        # Methods should be positively correlated (agree on relationships)
        assert corr_va > 0.5, "VanRaden and Astle methods should be correlated"
        assert corr_vy > 0.5, "VanRaden and Yang methods should be correlated"
        assert corr_ay > 0.5, "Astle and Yang methods should be correlated"

    def test_kinship_with_missing_data(self):
        """Test kinship computation with missing data."""
        genotypes = [
            [0, 1, -1, 2, 0],  # Sample 0 with missing
            [0, 1, 2, 0, 1],
            [2, -1, 0, 2, 1],  # Sample 2 with missing
        ]

        for method in ["vanraden", "astle", "yang"]:
            result = compute_kinship_matrix(genotypes, method=method)
            assert result["status"] == "success"
            
            kinship = np.array(result["kinship_matrix"])
            
            # Should still be symmetric
            assert np.allclose(kinship, kinship.T)
            
            # Should have missing data stats
            assert "missing_data_stats" in result
            assert result["missing_data_stats"]["total_missing"] > 0

    def test_kinship_matrix_shape(self):
        """Test that kinship matrix has correct shape."""
        genotypes = [
            [0, 1, 2],
            [0, 1, 2],
            [2, 1, 0],
            [1, 1, 1],
        ]

        result = compute_kinship_matrix(genotypes, method="vanraden")
        assert result["status"] == "success"
        
        kinship = np.array(result["kinship_matrix"])
        num_samples = len(genotypes)
        
        # Should be square matrix with shape (num_samples, num_samples)
        assert kinship.shape == (num_samples, num_samples)
        assert result["num_samples"] == num_samples

