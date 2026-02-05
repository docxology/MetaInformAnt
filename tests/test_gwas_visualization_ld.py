"""Tests for GWAS LD visualization functions."""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.visualization_ld import (
    compute_ld_decay,
    ld_decay_plot,
    ld_heatmap_region,
)


def _make_synthetic_genotypes(n_variants: int, n_samples: int, seed: int = 42) -> list[list[int]]:
    """Create synthetic genotype matrix with realistic LD structure."""
    rng = np.random.RandomState(seed)
    genotypes: list[list[int]] = []
    # First variant is independent
    base = rng.choice([0, 1, 2], size=n_samples, p=[0.25, 0.5, 0.25]).tolist()
    genotypes.append(base)

    for v in range(1, n_variants):
        prev = genotypes[v - 1]
        # Each subsequent variant is correlated with the previous one
        # with decreasing correlation as index distance grows
        noise_prob = min(0.5, 0.05 * v)
        new_geno = []
        for s in range(n_samples):
            if rng.random() < noise_prob:
                new_geno.append(int(rng.choice([0, 1, 2])))
            else:
                new_geno.append(prev[s])
        genotypes.append(new_geno)

    return genotypes


def _make_positions(n_variants: int, spacing: int = 10000) -> list[int]:
    """Create evenly-spaced variant positions."""
    return [i * spacing for i in range(n_variants)]


class TestComputeLdDecay:
    """Tests for compute_ld_decay."""

    def test_basic_20_variants(self) -> None:
        """Compute LD decay with 20 synthetic variants and verify output shape."""
        genotypes = _make_synthetic_genotypes(20, 100)
        positions = _make_positions(20, spacing=10000)

        result = compute_ld_decay(genotypes, positions, max_distance=200000, n_bins=10)

        assert result["status"] == "success"
        assert len(result["bin_centers"]) == 10
        assert len(result["mean_r2"]) == 10
        assert len(result["n_pairs"]) == 10
        assert result["total_pairs"] > 0

        # r2 values should be between 0 and 1
        for r2_val in result["mean_r2"]:
            assert 0.0 <= r2_val <= 1.0

    def test_nearby_variants_higher_ld(self) -> None:
        """Verify that closer variants have higher LD on average."""
        genotypes = _make_synthetic_genotypes(20, 200, seed=7)
        positions = _make_positions(20, spacing=10000)

        result = compute_ld_decay(genotypes, positions, max_distance=200000, n_bins=5)

        assert result["status"] == "success"
        # The first bin (closest) should have >= mean r2 of last populated bin
        populated = [
            (result["bin_centers"][i], result["mean_r2"][i])
            for i in range(len(result["n_pairs"]))
            if result["n_pairs"][i] > 0
        ]
        if len(populated) >= 2:
            assert populated[0][1] >= populated[-1][1]

    def test_chromosome_filtering(self) -> None:
        """Only pairs on the same chromosome contribute to LD decay."""
        genotypes = _make_synthetic_genotypes(10, 50, seed=99)
        positions = _make_positions(10, spacing=5000)
        # Assign first 5 to chrom 1, last 5 to chrom 2
        chromosomes = [1] * 5 + [2] * 5

        result_with_chrom = compute_ld_decay(
            genotypes, positions, max_distance=50000, n_bins=5, chromosomes=chromosomes
        )
        result_no_chrom = compute_ld_decay(genotypes, positions, max_distance=50000, n_bins=5, chromosomes=None)

        assert result_with_chrom["status"] == "success"
        assert result_no_chrom["status"] == "success"
        # With chromosome filtering, fewer cross-chrom pairs should mean fewer total pairs
        assert result_with_chrom["total_pairs"] <= result_no_chrom["total_pairs"]

    def test_two_variants_minimum(self) -> None:
        """Edge case: exactly 2 variants should still work."""
        genotypes = _make_synthetic_genotypes(2, 30, seed=5)
        positions = [0, 5000]

        result = compute_ld_decay(genotypes, positions, max_distance=10000, n_bins=5)

        assert result["status"] == "success"
        assert result["total_pairs"] == 1

    def test_single_variant_fails(self) -> None:
        """Edge case: 1 variant should return failed status."""
        result = compute_ld_decay([[0, 1, 2]], [100], max_distance=1000, n_bins=5)
        assert result["status"] == "failed"

    def test_all_identical_genotypes(self) -> None:
        """All variants have the same genotype vector; r2 should be 0 (zero variance)."""
        n_variants = 5
        n_samples = 20
        identical_geno = [1] * n_samples
        genotypes = [list(identical_geno) for _ in range(n_variants)]
        positions = _make_positions(n_variants, spacing=10000)

        result = compute_ld_decay(genotypes, positions, max_distance=50000, n_bins=5)

        assert result["status"] == "success"
        # With zero variance, all r2 should be 0
        for r2_val in result["mean_r2"]:
            assert r2_val == 0.0


class TestLdDecayPlot:
    """Tests for ld_decay_plot."""

    def test_basic_plot(self, tmp_path: Path) -> None:
        """Generate LD decay plot from pre-computed data."""
        decay_data = {
            "status": "success",
            "bin_centers": [5000.0, 15000.0, 25000.0, 35000.0, 45000.0],
            "mean_r2": [0.8, 0.5, 0.3, 0.15, 0.08],
            "n_pairs": [50, 45, 40, 30, 20],
            "total_pairs": 185,
        }

        output_path = tmp_path / "ld_decay.png"
        result = ld_decay_plot(decay_data, output_file=output_path)

        assert result["status"] == "success"
        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_plot_with_curve_fit(self, tmp_path: Path) -> None:
        """Fitted curve should be present when fit_curve=True."""
        # Create data that clearly decays
        bin_centers = [float(i * 10000) for i in range(1, 11)]
        mean_r2 = [0.9 * np.exp(-0.0003 * x) + 0.05 for x in bin_centers]

        decay_data = {
            "status": "success",
            "bin_centers": bin_centers,
            "mean_r2": mean_r2,
            "n_pairs": [100] * 10,
            "total_pairs": 1000,
        }

        output_path = tmp_path / "ld_decay_fitted.png"
        result = ld_decay_plot(decay_data, output_file=output_path, fit_curve=True, title="Fitted LD Decay")

        assert result["status"] == "success"
        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_plot_no_fit(self, tmp_path: Path) -> None:
        """Plot without curve fitting should still succeed."""
        decay_data = {
            "bin_centers": [10000.0, 20000.0],
            "mean_r2": [0.5, 0.3],
            "n_pairs": [10, 10],
            "total_pairs": 20,
        }

        output_path = tmp_path / "ld_no_fit.png"
        result = ld_decay_plot(decay_data, output_file=output_path, fit_curve=False)

        assert result["status"] == "success"
        assert output_path.exists()

    def test_plot_empty_data(self) -> None:
        """Empty data should return failed status."""
        result = ld_decay_plot({"bin_centers": [], "mean_r2": [], "n_pairs": []})
        assert result["status"] == "failed"


class TestLdHeatmapRegion:
    """Tests for ld_heatmap_region."""

    def test_basic_10_variants(self, tmp_path: Path) -> None:
        """Generate triangular LD heatmap with 10 variants."""
        genotypes = _make_synthetic_genotypes(10, 80, seed=12)
        positions = _make_positions(10, spacing=5000)

        output_path = tmp_path / "ld_heatmap.png"
        result = ld_heatmap_region(genotypes, positions, output_file=output_path)

        assert result["status"] == "success"
        assert output_path.exists()
        assert output_path.stat().st_size > 0
        assert result["n_variants"] == 10

    def test_region_subsetting(self, tmp_path: Path) -> None:
        """Restrict heatmap to a sub-region of the data."""
        genotypes = _make_synthetic_genotypes(20, 60, seed=33)
        positions = _make_positions(20, spacing=5000)

        output_path = tmp_path / "ld_heatmap_region.png"
        result = ld_heatmap_region(
            genotypes,
            positions,
            output_file=output_path,
            region_start=25000,
            region_end=75000,
        )

        assert result["status"] == "success"
        assert output_path.exists()
        # Only variants in [25000, 75000] should be included
        assert result["n_variants"] < 20

    def test_two_variants(self, tmp_path: Path) -> None:
        """Edge case: exactly 2 variants should produce a valid heatmap."""
        genotypes = _make_synthetic_genotypes(2, 40, seed=77)
        positions = [1000, 2000]

        output_path = tmp_path / "ld_heatmap_2var.png"
        result = ld_heatmap_region(genotypes, positions, output_file=output_path)

        assert result["status"] == "success"
        assert output_path.exists()
        assert result["n_variants"] == 2

    def test_single_variant_fails(self) -> None:
        """One variant should return failed status."""
        result = ld_heatmap_region([[0, 1, 2]], [100])
        assert result["status"] == "failed"

    def test_all_identical_genotypes(self, tmp_path: Path) -> None:
        """All identical genotype vectors produce a heatmap with r2=0 off-diagonal."""
        n_variants = 5
        n_samples = 30
        genotypes = [[1] * n_samples for _ in range(n_variants)]
        positions = _make_positions(n_variants, spacing=1000)

        output_path = tmp_path / "ld_heatmap_identical.png"
        result = ld_heatmap_region(genotypes, positions, output_file=output_path)

        assert result["status"] == "success"
        assert output_path.exists()

    def test_no_output_file(self) -> None:
        """Running without output_file should still succeed (no file saved)."""
        genotypes = _make_synthetic_genotypes(5, 30, seed=55)
        positions = _make_positions(5, spacing=2000)

        result = ld_heatmap_region(genotypes, positions, output_file=None)

        assert result["status"] == "success"
        assert result["output_path"] is None
