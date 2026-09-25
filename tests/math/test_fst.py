"""Pinned-value tests for the moment F_ST estimator and Weir & Cockerham theta."""

from __future__ import annotations

import pytest

from metainformant.math.population_genetics.fst import (
    fst_from_allele_freq_matrix,
    fst_from_allele_freqs,
    weirs_fst,
)


class TestMomentFst:
    def test_single_locus_matches_heterozygosity_estimator(self) -> None:
        # (Ht - Hs) / Ht with Hs = 0.32, Ht = 0.5 for (0.2, 0.8)
        assert fst_from_allele_freqs([0.2, 0.8]) == pytest.approx(0.36)

    def test_single_locus_agrees_with_multi_locus_for_one_locus(self) -> None:
        # The single-list call is one locus, two populations: it must equal
        # the multi-locus path applied to the same data.
        for p1, p2 in [(0.2, 0.8), (0.4, 0.6), (0.05, 0.95), (0.5, 0.5), (1.0, 0.0)]:
            single = fst_from_allele_freqs([p1, p2])
            multi = fst_from_allele_freqs([p1], [p2])
            assert single == pytest.approx(multi)

    def test_multi_locus_hand_computed_value(self) -> None:
        pop1 = [0.6, 0.4, 0.8]
        pop2 = [0.3, 0.7, 0.2]
        # Per-locus terms (var_p; var_p + mean p(1-p)):
        # (0.6, 0.3) -> 0.0225, 0.2475; (0.4, 0.7) -> 0.0225, 0.2475;
        # (0.8, 0.2) -> 0.09, 0.25.  Pooled: 0.135 / 0.745.
        assert fst_from_allele_freqs(pop1, pop2) == pytest.approx(0.135 / 0.745)

    def test_multi_locus_not_constant_two_thirds(self) -> None:
        # Regression: the previous implementation returned 2/3 for any
        # differentiated input because the within/between variance terms
        # were both proportional to (p1 - p2)^2.
        values = {
            fst_from_allele_freqs([0.51, 0.49], [0.49, 0.51]),
            fst_from_allele_freqs([0.2, 0.8], [0.8, 0.2]),
            fst_from_allele_freqs([1.0, 0.0], [0.0, 1.0]),
        }
        assert len(values) == 3
        assert 2.0 / 3.0 not in values

    def test_three_populations_use_k_over_k_minus_one_correction(self) -> None:
        # Locus 1 (0.2, 0.5, 0.8): var_p = 0.18 / (3 - 1) = 0.09 and
        # mean p(1-p) = 0.19 -> denominator 0.28.
        # Locus 2 (identical everywhere): var_p = 0, denominator 0.25.
        pops = [[0.2, 0.5], [0.5, 0.5], [0.8, 0.5]]
        assert fst_from_allele_freq_matrix(pops) == pytest.approx(0.09 / 0.53)

    def test_identical_populations_return_zero(self) -> None:
        assert fst_from_allele_freq_matrix([[0.3, 0.6], [0.3, 0.6], [0.3, 0.6]]) == 0.0

    def test_fully_differentiated_two_populations_return_one(self) -> None:
        assert fst_from_allele_freqs([1.0, 0.0]) == 1.0

    def test_zero_variance_returns_zero(self) -> None:
        assert fst_from_allele_freq_matrix([[0.5], [0.5]]) == 0.0
        assert fst_from_allele_freqs([1.0, 1.0]) == 0.0

    def test_matrix_validation(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            fst_from_allele_freq_matrix([[0.5]])
        with pytest.raises(ValueError, match="cannot be empty"):
            fst_from_allele_freq_matrix([[], []])
        with pytest.raises(ValueError, match="expected 2"):
            fst_from_allele_freq_matrix([[0.5, 0.5], [0.5]])

    def test_allele_freqs_validation(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            fst_from_allele_freqs([0.5, 0.2], [0.5])
        with pytest.raises(ValueError, match="cannot be empty"):
            fst_from_allele_freqs([], [])
        with pytest.raises(ValueError, match="Invalid frequency"):
            fst_from_allele_freqs([0.5, 1.2], [0.5, 0.2])
        with pytest.raises(ValueError, match="exactly 2"):
            fst_from_allele_freqs([0.2, 0.3, 0.4])
        with pytest.raises(ValueError, match="Invalid frequency"):
            fst_from_allele_freqs([0.5, -0.1])


class TestWeirsFst:
    def test_symmetric_two_population_hand_value(self) -> None:
        counts = {"pop1": {"A": 8, "G": 2}, "pop2": {"A": 2, "G": 8}}
        # n = (10, 10); allele A: p = (0.8, 0.2), p_bar = 0.5, s2 = 0.18
        # a = 0.18 - (1/9) * (0.25 - 0.5 * 0.18) = 1.46 / 9
        # b = (10/9) * (0.25 - 0.5 * 0.18) = 1.6 / 9, c = 0 (haploid counts)
        # theta = a / (a + b) = 1.46 / 3.06; allele G contributes identically.
        assert weirs_fst(counts) == pytest.approx(1.46 / 3.06)

    def test_fixed_different_populations(self) -> None:
        counts = {"pop1": {"A": 5}, "pop2": {"G": 5}}
        # s2 = 0.5, het = 0.25 -> a = 0.5, b = 0, c = 0 -> theta = 1
        assert weirs_fst(counts) == pytest.approx(1.0)

    def test_identical_populations_zero(self) -> None:
        counts = {"pop1": {"A": 4, "G": 4}, "pop2": {"A": 4, "G": 4}}
        assert weirs_fst(counts) == 0.0

    def test_three_populations(self) -> None:
        counts = {
            "p1": {"A": 10, "G": 2},
            "p2": {"A": 6, "G": 6},
            "p3": {"A": 2, "G": 10},
        }
        # Equal n = (12, 12, 12): n_bar = 12, n_c = (36 - 432/36) / 2 = 12.
        # Allele A: p = (10/12, 1/2, 2/12), p_bar = 0.5
        # s2 = [12*(1/3)^2 + 12*(1/3)^2] / (2 * 12) = 1/9
        # a = 1/9 - (1/11) * (1/4 - (2/3) * (1/9)) = 113/1188
        # b = (12/11) * (1/4 - (2/3) * (1/9)) = 228/1188, c = 0
        # theta = a / (a + b) = 113/341 (allele G contributes identically)
        assert weirs_fst(counts) == pytest.approx(113 / 341)
        # The value must not depend on population or haplotype key order.
        reordered = {
            "p3": {"G": 10, "A": 2},
            "p1": {"G": 2, "A": 10},
            "p2": {"G": 6, "A": 6},
        }
        assert weirs_fst(reordered) == pytest.approx(113 / 341)

    def test_unequal_sample_sizes_hand_value(self) -> None:
        counts = {"pop1": {"A": 8, "G": 2}, "pop2": {"A": 15, "G": 15}}
        # n = (10, 30), n_bar = 20, n_c = (40 - (100 + 900) / 40) / 1 = 15
        # allele A: p = (0.8, 0.5), p_bar = (10 * 0.8 + 30 * 0.5) / 40 = 0.575
        # s2 = (10 * 0.225^2 + 30 * 0.075^2) / ((2 - 1) * 20) = 0.03375
        # het = 0.575 * 0.425 = 0.244375
        # a = (20/15) * (0.03375 - (1/19) * (0.244375 - 0.5 * 0.03375))
        # b = (20/19) * (0.244375 - 0.5 * 0.03375)
        a = (20 / 15) * (0.03375 - (1 / 19) * (0.244375 - 0.5 * 0.03375))
        b = (20 / 19) * (0.244375 - 0.5 * 0.03375)
        assert weirs_fst(counts) == pytest.approx(a / (a + b))

    def test_edge_cases_return_zero(self) -> None:
        assert weirs_fst({}) == 0.0
        assert weirs_fst({"pop1": {"A": 3}}) == 0.0
        assert weirs_fst({"pop1": {"A": 0}, "pop2": {"A": 0}}) == 0.0
        assert weirs_fst({"pop1": {"A": 1}, "pop2": {"A": 1}}) == 0.0  # n_bar <= 1

    def test_negative_counts_raise(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            weirs_fst({"pop1": {"A": -1}, "pop2": {"A": 3}})
