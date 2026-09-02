"""Zero-mocks tests for gwas.analysis.ldsr (LD score computation + regression)."""

from __future__ import annotations

import random

import pytest

from metainformant.gwas.analysis import ldsr


class TestComputeLDScores:
    def test_identical_variants_full_ld_score(self) -> None:
        # Two identical dosage rows: r²=1 with each other → each LD score = 2 (self + other)
        geno = [[0, 1, 2, 1, 0], [0, 1, 2, 1, 0]]
        scores = ldsr.compute_ld_scores(geno)
        assert len(scores) == 2
        assert scores[0] == pytest.approx(2.0)
        assert scores[1] == pytest.approx(2.0)

    def test_uncorrelated_variants_ld_score_one(self) -> None:
        rng = random.Random(42)
        geno = [rng.choice([0, 1, 2]) for _ in range(500)]
        rows = [list(geno), [rng.choice([0, 1, 2]) for _ in range(500)]]
        scores = ldsr.compute_ld_scores(rows)
        # r²≈0 between rows → score ≈ 1.0 (self term only)
        assert scores[0] < 1.1

    def test_window_excludes_distant_variants(self) -> None:
        geno = [[0, 1, 2, 1, 0], [0, 1, 2, 1, 0]]
        positions = [1000, 5_000_000]  # 5 Mb apart
        scores = ldsr.compute_ld_scores(geno, positions=positions, window_kb=1000)
        # Outside window: only the self term contributes
        assert scores[0] == pytest.approx(1.0)
        assert scores[1] == pytest.approx(1.0)

    def test_single_variant_score_one(self) -> None:
        assert ldsr.compute_ld_scores([[0, 1, 2]]) == [1.0]

    def test_scores_nonnegative(self) -> None:
        rng = random.Random(7)
        geno = [[rng.choice([0, 1, 2]) for _ in range(30)] for _ in range(10)]
        scores = ldsr.compute_ld_scores(geno)
        assert all(s >= 1.0 - 1e-9 for s in scores)


class TestLDSRRegression:
    def test_recovers_h2_from_synthetic_chi2(self) -> None:
        """Generate chi2 = 1 + (N h² / M) ℓ_j + noise; regression must recover h²."""
        rng = random.Random(13)
        n_samples, n_variants = 10_000, 2000
        h2_true = 0.3
        ld_scores = [rng.uniform(1.0, 30.0) for _ in range(n_variants)]
        slope = n_samples * h2_true / n_variants
        chi2 = [1.0 + slope * s + rng.gauss(0, 0.5) for s in ld_scores]
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=n_samples, n_variants=n_variants)
        assert result["method"] != "insufficient_data"
        assert abs(result["h2_snp"] - h2_true) < 0.05
        assert abs(result["intercept"] - 1.0) < 0.2
        assert result["h2_snp_se"] > 0

    def test_null_statistics_intercept_one_zero_h2(self) -> None:
        rng = random.Random(17)
        n_variants = 500
        ld_scores = [rng.uniform(1.0, 20.0) for _ in range(n_variants)]
        chi2 = [1.0 + rng.gauss(0, 0.3) for _ in range(n_variants)]  # pure confound-free null
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=5000, n_variants=n_variants)
        assert abs(result["h2_snp"]) < 0.05
        assert abs(result["intercept"] - 1.0) < 0.15

    def test_intercept_constraint(self) -> None:
        rng = random.Random(23)
        n_variants = 300
        ld_scores = [rng.uniform(1.0, 15.0) for _ in range(n_variants)]
        chi2 = [1.2 + 0.01 * s for s in ld_scores]
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=8000, n_variants=n_variants, intercept_constraint=1.2)
        assert result["intercept"] == pytest.approx(1.2)

    def test_insufficient_data(self) -> None:
        result = ldsr.ldsr_regression([2.0, 1.0], [5.0, 3.0], n_samples=1000, n_variants=100)
        assert result["method"] == "insufficient_data"
        assert result["h2_snp"] == 0.0

    def test_output_keys(self) -> None:
        rng = random.Random(29)
        ld = [rng.uniform(1.0, 10.0) for _ in range(50)]
        chi2 = [1.0 + 0.02 * s + rng.gauss(0, 0.2) for s in ld]
        result = ldsr.ldsr_regression(chi2, ld, n_samples=2000, n_variants=50)
        for key in ("h2_snp", "intercept", "mean_chi2", "method"):
            assert key in result
        assert result["mean_chi2"] == pytest.approx(sum(chi2) / len(chi2), rel=1e-3)
