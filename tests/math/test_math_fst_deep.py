"""Deep tests for metainformant.math.population_genetics.fst (real computation, no test doubles)."""

import math

import numpy as np
import pytest

from metainformant.math.population_genetics.fst import (
    fst_confidence_interval,
    fst_from_allele_freqs,
    fst_from_heterozygosity,
    pairwise_fst_matrix,
    weirs_fst,
)


class TestFstFromAlleleFreqs:
    def test_single_locus_known_value(self) -> None:
        # p1=0.2, p2=0.8: Hs = (2*0.16 + 2*0.16)/2 = 0.32; Ht = 2*0.5*0.5 = 0.5
        fst = fst_from_allele_freqs([0.2, 0.8])
        assert fst == pytest.approx((0.5 - 0.32) / 0.5)

    def test_single_locus_identical_populations_zero(self) -> None:
        assert fst_from_allele_freqs([0.5, 0.5]) == 0.0

    def test_single_locus_fixed_different_fst_one(self) -> None:
        # p1=1.0, p2=0.0: Hs=0, Ht=1.0 -> Fst=1
        assert fst_from_allele_freqs([1.0, 0.0]) == 1.0

    def test_single_locus_fixed_identical_returns_zero(self) -> None:
        # Both populations fixed for same allele: Ht = 0 -> guard returns 0.0
        assert fst_from_allele_freqs([1.0, 1.0]) == 0.0

    def test_single_locus_wrong_length_raises(self) -> None:
        with pytest.raises(ValueError, match="exactly 2"):
            fst_from_allele_freqs([0.2, 0.3, 0.4])

    def test_multi_locus_symmetric_matrix_consistency(self) -> None:
        pop1 = [0.6, 0.4, 0.8]
        pop2 = [0.3, 0.7, 0.2]
        fst = fst_from_allele_freqs(pop1, pop2)
        fst_swapped = fst_from_allele_freqs(pop2, pop1)
        assert fst == pytest.approx(fst_swapped)
        assert 0.0 <= fst <= 1.0

    def test_multi_locus_identical_populations_zero(self) -> None:
        assert fst_from_allele_freqs([0.5, 0.2], [0.5, 0.2]) == 0.0

    def test_multi_locus_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            fst_from_allele_freqs([0.5, 0.2], [0.5])

    def test_multi_locus_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="cannot be empty"):
            fst_from_allele_freqs([], [])

    def test_multi_locus_invalid_frequency_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid frequency"):
            fst_from_allele_freqs([0.5, 1.2], [0.5, 0.2])


class TestPairwiseFstMatrix:
    def test_shape_and_symmetry(self) -> None:
        pops = [[0.6, 0.4, 0.8], [0.3, 0.7, 0.2], [0.5, 0.5, 0.6]]
        m = pairwise_fst_matrix(pops)
        assert m.shape == (3, 3)
        assert np.allclose(m, m.T)
        assert np.allclose(np.diag(m), 0.0)

    def test_diagonal_zero_offdiagonal_matches_direct(self) -> None:
        pops = [[0.6, 0.4], [0.3, 0.7], [0.9, 0.1]]
        m = pairwise_fst_matrix(pops)
        for i in range(3):
            for j in range(i + 1, 3):
                expected = fst_from_allele_freqs(pops[i], pops[j])
                assert m[i, j] == pytest.approx(expected)

    def test_fully_different_populations(self) -> None:
        # Populations fixed for opposite alleles: Hs = 0, Ht = 1 -> Fst = 1.
        # (Regression: the removed variance-ratio implementation returned 2/3.)
        m = pairwise_fst_matrix([[1.0, 0.0], [0.0, 1.0]])
        assert m[0, 1] == pytest.approx(1.0)

    def test_too_few_populations_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            pairwise_fst_matrix([[0.5, 0.5]])

    def test_unequal_loci_raises(self) -> None:
        with pytest.raises(ValueError, match="expected"):
            pairwise_fst_matrix([[0.5, 0.5], [0.5]])


class TestWeirsFst:
    def test_symmetric_two_population_hand_value(self) -> None:
        counts = {"pop1": {"A": 8, "G": 2}, "pop2": {"A": 2, "G": 8}}
        # n = (10, 10); allele A: p = (0.8, 0.2), p_bar = 0.5, s2 = 0.18
        # a = 0.18 - (1/9) * (0.25 - 0.5 * 0.18) = 1.46 / 9
        # b = (10/9) * (0.25 - 0.5 * 0.18) = 1.6 / 9, c = 0 (haploid counts)
        # theta = sum_u a_u / sum_u (a_u + b_u + c_u) = 1.46 / 3.06
        assert weirs_fst(counts) == pytest.approx(1.46 / 3.06)

    def test_fixed_different_populations(self) -> None:
        assert weirs_fst({"pop1": {"A": 5}, "pop2": {"G": 5}}) == pytest.approx(1.0)

    def test_identical_populations_zero(self) -> None:
        counts = {"pop1": {"A": 4, "G": 4}, "pop2": {"A": 4, "G": 4}}
        assert weirs_fst(counts) == 0.0

    def test_single_population_returns_zero(self) -> None:
        assert weirs_fst({"pop1": {"AT": 10}}) == 0.0

    def test_empty_inputs_return_zero(self) -> None:
        assert weirs_fst({}) == 0.0
        assert weirs_fst({"pop1": {"AT": 0}, "pop2": {"AT": 0}}) == 0.0

    def test_negative_counts_raise(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            weirs_fst({"pop1": {"A": -1}, "pop2": {"A": 3}})


class TestFstConfidenceInterval:
    def test_contains_point_estimate(self) -> None:
        lo, hi = fst_confidence_interval(0.2, 50)
        assert lo <= 0.2 <= hi

    def test_known_z_95_bounds(self) -> None:
        fst, n = 0.2, 50
        variance = (2 * fst**2 * (1 - fst) ** 2) / n
        se = math.sqrt(variance)
        lo, hi = fst_confidence_interval(fst, n, 0.95)
        assert lo == pytest.approx(max(0.0, fst - 1.96 * se))
        assert hi == pytest.approx(min(1.0, fst + 1.96 * se))

    def test_bounds_within_unit_interval(self) -> None:
        lo, hi = fst_confidence_interval(0.999, 3)
        assert lo >= 0.0
        assert hi <= 1.0

    def test_higher_confidence_wider_interval(self) -> None:
        lo95, hi95 = fst_confidence_interval(0.3, 100, 0.95)
        lo99, hi99 = fst_confidence_interval(0.3, 100, 0.99)
        assert hi99 - lo99 >= hi95 - lo95

    def test_small_sample_correction_widens_interval(self) -> None:
        lo_big, hi_big = fst_confidence_interval(0.3, 200)
        lo_small, hi_small = fst_confidence_interval(0.3, 10)
        assert hi_small - lo_small > hi_big - lo_big

    def test_sample_size_too_small_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            fst_confidence_interval(0.2, 1)

    def test_nonstandard_confidence_uses_scipy(self) -> None:
        from scipy import stats

        lo, hi = fst_confidence_interval(0.2, 100, 0.80)
        z = stats.norm.ppf(0.90)
        assert lo == pytest.approx(0.2 - z * math.sqrt((2 * 0.04 * 0.64) / 100))


class TestFstFromHeterozygosity:
    def test_known_value(self) -> None:
        assert fst_from_heterozygosity(0.2, 0.5) == pytest.approx(0.6)

    def test_zero_ht_returns_zero(self) -> None:
        assert fst_from_heterozygosity(0.0, 0.0) == 0.0

    def test_result_clamped_to_unit_interval(self) -> None:
        # Hs > Ht would give negative; must clamp to 0
        assert fst_from_heterozygosity(0.9, 0.5) == 0.0

    def test_complete_differentiation(self) -> None:
        assert fst_from_heterozygosity(0.0, 0.5) == 1.0


class TestFstFromHeterozygosityClamping:
    def test_hs_greater_than_ht_clamps_to_zero(self) -> None:
        assert fst_from_heterozygosity(0.6, 0.5) == 0.0


class TestFstConfidenceIntervalEdges:
    def test_sample_size_below_two_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            fst_confidence_interval(0.2, 1)

    def test_ninety_percent_interval_narrower_than_ninety_five(self) -> None:
        lo90, hi90 = fst_confidence_interval(0.2, 50, confidence_level=0.90)
        lo95, hi95 = fst_confidence_interval(0.2, 50, confidence_level=0.95)
        assert lo90 <= 0.2 <= hi90
        assert (hi90 - lo90) < (hi95 - lo95)
