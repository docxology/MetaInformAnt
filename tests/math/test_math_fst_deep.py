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
        # Variance-based estimator: var_between=(1+1)/2=1, var_within=0... but per-locus
        # variance term ((1-0.5)^2+(0-0.5)^2)/2 = 0.25 -> Fst = 1/(1+... ) measured 2/3.
        m = pairwise_fst_matrix([[1.0, 0.0], [0.0, 1.0]])
        assert m[0, 1] == pytest.approx(2.0 / 3.0)

    def test_too_few_populations_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            pairwise_fst_matrix([[0.5, 0.5]])

    def test_unequal_loci_raises(self) -> None:
        with pytest.raises(ValueError, match="expected"):
            pairwise_fst_matrix([[0.5, 0.5], [0.5]])


class TestWeirsFst:
    def test_basic_two_populations(self) -> None:
        counts = {"AT": 10, "AG": 15, "GT": 8, "GG": 12}
        labels = ["pop1"] * 22 + ["pop2"] * 23
        fst = weirs_fst(counts, labels)
        assert 0.0 <= fst <= 1.0

    def test_single_population_returns_zero(self) -> None:
        assert weirs_fst({"AT": 10}, ["pop1"] * 10) == 0.0

    def test_empty_inputs_return_zero(self) -> None:
        assert weirs_fst({}, []) == 0.0
        assert weirs_fst({"AT": 5}, []) == 0.0

    def test_zero_total_counts_return_zero(self) -> None:
        assert weirs_fst({"AT": 0}, ["pop1", "pop2"]) == 0.0

    def test_fst_is_deterministic_within_process(self) -> None:
        counts = {"AT": 10, "AG": 15, "GT": 8, "GG": 12}
        labels = ["pop1"] * 22 + ["pop2"] * 23
        f1 = weirs_fst(counts, labels)
        f2 = weirs_fst(dict(counts), list(labels))
        assert f1 == pytest.approx(f2)


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
