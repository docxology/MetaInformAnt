"""Zero-mocks tests for gwas.analysis.prs (C+T clumping, scoring, R²)."""

from __future__ import annotations

import random

import pytest

from metainformant.gwas.analysis import prs


class TestClumpVariants:
    def test_keeps_independent_variants(self) -> None:
        # 3 variants in perfect LD (r²=1) plus 2 independent ones.
        results = [{"p_value": p} for p in [1e-8, 1e-7, 1e-9, 1e-3, 1e-4]]
        ld = [
            [1.0, 1.0, 1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]
        reps = prs.clump_variants(results, ld, r2_threshold=0.1)
        # Best in the LD block is index 2 (p=1e-9); 3 (1e-3) and 4 (1e-4) are
        # independent. Representatives come back in p-ascending order.
        assert reps == [2, 4, 3]

    def test_all_in_ld_keeps_single_best(self) -> None:
        n = 5
        results = [{"p_value": 10 ** (-i - 4)} for i in range(n)]  # index 4 best (p=1e-8)
        ld = [[1.0] * n for _ in range(n)]
        reps = prs.clump_variants(results, ld, r2_threshold=0.1)
        assert reps == [4]

    def test_p_threshold_prefilter(self) -> None:
        results = [{"p_value": p} for p in [1e-9, 0.5, 0.9]]
        ld = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        reps = prs.clump_variants(results, ld, p_threshold=0.05)
        assert reps == [0]

    def test_empty_input(self) -> None:
        assert prs.clump_variants([], []) == []

    def test_none_pass_threshold(self) -> None:
        results = [{"p_value": 0.5}, {"p_value": 0.6}]
        reps = prs.clump_variants(results, [[1.0, 0.0], [0.0, 1.0]], p_threshold=0.05)
        assert reps == []


class TestComputePRS:
    def test_single_variant_prs_equals_scaled_dosage(self) -> None:
        geno = [[0.0, 1.0, 2.0, 1.0]]
        betas = [0.5]
        scores = prs.compute_prs(geno, betas, [0], n_samples=4)
        # Normalized by number of used variants (1): score = beta * dosage
        assert scores == pytest.approx([0.0, 0.5, 1.0, 0.5])

    def test_two_variant_additive(self) -> None:
        geno = [[0.0, 2.0], [1.0, 1.0]]
        betas = [1.0, 2.0]
        scores = prs.compute_prs(geno, betas, [0, 1], n_samples=2)
        # sample0: (1*0 + 2*1)/2 = 1.0 ; sample1: (1*2 + 2*1)/2 = 2.0
        assert scores == pytest.approx([1.0, 2.0])

    def test_missing_dosages_skipped_not_averaged_in(self) -> None:
        geno = [[float("nan"), 2.0]]
        betas = [1.0]
        scores = prs.compute_prs(geno, betas, [0], n_samples=2)
        assert scores[0] == 0.0
        assert scores[1] == pytest.approx(2.0)

    def test_out_of_range_indices_skipped(self) -> None:
        geno = [[0.0, 1.0]]
        scores = prs.compute_prs(geno, [0.5], [5], n_samples=2)
        assert scores == [0.0, 0.0]

    def test_empty_selection_zero_scores(self) -> None:
        assert prs.compute_prs([[0, 1, 2]], [1.0], [], n_samples=3) == [0.0, 0.0, 0.0]


class TestPrsR2:
    def test_perfect_prediction_r2_one(self) -> None:
        prs_scores = [1.0, 2.0, 3.0, 4.0, 5.0]
        phenotypes = [10.0, 20.0, 30.0, 40.0, 50.0]
        assert prs.prs_r2(prs_scores, phenotypes) == pytest.approx(1.0, abs=1e-9)

    def test_anticorrelated_still_positive_r2(self) -> None:
        prs_scores = [1.0, 2.0, 3.0, 4.0, 5.0]
        phenotypes = [-10.0, -20.0, -30.0, -40.0, -50.0]
        assert prs.prs_r2(prs_scores, phenotypes) == pytest.approx(1.0, abs=1e-9)

    def test_uncorrelated_low_r2(self) -> None:
        rng = random.Random(3)
        prs_scores = [rng.random() for _ in range(200)]
        phenotypes = [rng.random() for _ in range(200)]
        r2 = prs.prs_r2(prs_scores, phenotypes)
        assert 0.0 <= r2 < 0.05

    def test_degenerate_inputs_zero(self) -> None:
        assert prs.prs_r2([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == 0.0
        assert prs.prs_r2([1.0], [1.0]) == 0.0
        assert prs.prs_r2([], []) == 0.0


class TestPrsFullAnalysis:
    def test_end_to_end_on_constructed_data(self) -> None:
        """PRS on a causal variant should predict the phenotype it generated."""
        rng = random.Random(11)
        n_samples = 200
        causal_dosage = [rng.choice([0, 1, 2]) for _ in range(n_samples)]
        beta_true = 1.5
        phenotypes = [beta_true * g + rng.gauss(0, 0.5) for g in causal_dosage]
        # Two variants: index 0 is the causal one; index 1 is noise.
        geno = [causal_dosage, [rng.choice([0, 1, 2]) for _ in range(n_samples)]]
        association_results = [{"p_value": 1e-12, "beta": beta_true}, {"p_value": 0.4, "beta": 0.0}]
        ld = [[1.0, 0.0], [0.0, 1.0]]
        result = prs.prs_full_analysis(association_results, geno, phenotypes, ld, p_thresholds=[5e-8, 1e-4, 0.05])
        assert isinstance(result, dict)
        assert result["best_r2"] > 0.7
        assert result["best_threshold"] == 5e-8  # causal variant included
        assert len(result["thresholds"]) == 3
        for row in result["thresholds"]:
            assert {"p_threshold", "n_variants", "r2_phenotype", "prs_scores"} <= set(row.keys())
