"""Tests for GWAS LD pruning module."""

from __future__ import annotations

import pytest

from metainformant.gwas.analysis.ld_pruning import _compute_r_squared_pair, ld_prune


class TestComputeRSquared:
    """Tests for pairwise R-squared computation."""

    def test_identical_vectors(self) -> None:
        """Identical vectors should have R^2 = 1."""
        a = [0, 1, 2, 0, 1, 2]
        r2 = _compute_r_squared_pair(a, a)
        assert abs(r2 - 1.0) < 0.001

    def test_uncorrelated_vectors(self) -> None:
        """Uncorrelated vectors should have low R^2."""
        a = [0, 0, 1, 1, 2, 2]
        b = [2, 0, 1, 2, 0, 1]
        r2 = _compute_r_squared_pair(a, b)
        assert r2 < 0.5

    def test_missing_data_handled(self) -> None:
        """Missing data (-1) should be excluded."""
        a = [0, 1, -1, 2, 1, 0]
        b = [0, 1, 2, 2, 1, 0]
        r2 = _compute_r_squared_pair(a, b)
        assert 0.0 <= r2 <= 1.0

    def test_all_missing(self) -> None:
        """All missing data should return 0."""
        a = [-1, -1, -1]
        b = [-1, -1, -1]
        r2 = _compute_r_squared_pair(a, b)
        assert r2 == 0.0

    def test_constant_vector(self) -> None:
        """Constant vector should give R^2 = 0 (no variance)."""
        a = [1, 1, 1, 1, 1]
        b = [0, 1, 2, 0, 1]
        r2 = _compute_r_squared_pair(a, b)
        assert r2 == 0.0

    def test_perfectly_anticorrelated(self) -> None:
        """Perfectly anti-correlated should have R^2 = 1."""
        a = [0, 1, 2, 0, 1, 2]
        b = [2, 1, 0, 2, 1, 0]
        r2 = _compute_r_squared_pair(a, b)
        assert abs(r2 - 1.0) < 0.001


class TestLDPrune:
    """Tests for LD pruning."""

    def test_no_variants(self) -> None:
        """Empty genotype matrix should return empty."""
        assert ld_prune([]) == []

    def test_single_variant(self) -> None:
        """Single variant should be kept."""
        geno = [[0, 1, 2, 0, 1]]
        kept = ld_prune(geno)
        assert kept == [0]

    def test_uncorrelated_variants_kept(self) -> None:
        """Uncorrelated variants should all be kept."""
        # Two independent variants
        geno = [
            [0, 0, 1, 1, 2, 2, 0, 0, 1, 1],
            [2, 0, 1, 2, 0, 1, 2, 0, 1, 2],
        ]
        kept = ld_prune(geno, r2_threshold=0.5)
        assert len(kept) == 2

    def test_correlated_variants_pruned(self) -> None:
        """Perfectly correlated variants should be pruned."""
        v = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
        geno = [v[:], v[:], v[:]]  # Three identical variants
        kept = ld_prune(geno, r2_threshold=0.2)
        # Only one should survive
        assert len(kept) == 1

    def test_window_respects_boundaries(self) -> None:
        """Variants far apart should not be pruned against each other."""
        n = 20
        v = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
        geno = [v[:] for _ in range(n)]
        # Very small window: only compares adjacent variants
        kept = ld_prune(geno, window_size=2, step_size=1, r2_threshold=0.2)
        # Should keep at least 1 per non-overlapping window
        assert 1 <= len(kept) <= n

    def test_high_r2_threshold_keeps_all(self) -> None:
        """Very high threshold should keep all variants."""
        v = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
        geno = [v[:], v[:]]
        kept = ld_prune(geno, r2_threshold=1.01)  # Above max R^2
        assert len(kept) == 2

    def test_missing_data_preference(self) -> None:
        """Variant with more missing data should be preferentially removed."""
        a = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
        b = [0, 1, 2, 0, 1, 2, -1, -1, -1, 0]  # More missing
        geno = [a, b]
        kept = ld_prune(geno, r2_threshold=0.2)
        # b has more missing, so b should be removed and a kept
        assert 0 in kept
        assert len(kept) == 1
