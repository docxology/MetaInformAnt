from __future__ import annotations

import pytest

from metainformant.math.core.utilities import r_squared
from metainformant.math.population_genetics.ld import ld_coefficients


def test_ld_coefficients_and_r2():
    # pA=0.6, pa=0.4, pB=0.7, pb=0.3, pAB=0.5
    D, Dp = ld_coefficients(0.6, 0.4, 0.7, 0.3, 0.5)
    assert abs(D - (0.5 - 0.42)) < 1e-12
    # Dmax = min(pA*pb, pa*pB) = min(0.18, 0.28) = 0.18 => D' = 0.08/0.18
    assert abs(Dp - ((0.08) / 0.18)) < 1e-12

    r2 = r_squared(0.6, 0.4, 0.7, 0.3, 0.5)
    denom = 0.6 * 0.4 * 0.7 * 0.3
    assert abs(r2 - ((0.08 * 0.08) / denom)) < 1e-12


def test_ld_coefficients_genotype_mode_perfect_ld():
    # Perfect coupling: only 11 and 00 haplotypes present
    genotypes = [[1, 1], [1, 1], [0, 0], [0, 0]]
    res = ld_coefficients(genotypes)
    assert res["D"] == pytest.approx(0.5 - 0.25)  # P(AB)=0.5, pA*pB=0.25
    assert res["D_prime"] == pytest.approx(1.0)
    assert res["r_squared"] == pytest.approx(1.0)


def test_ld_coefficients_genotype_mode_negative_ld():
    # Repulsion: only 10 and 01 haplotypes present
    genotypes = [[1, 0], [0, 1], [1, 0], [0, 1]]
    res = ld_coefficients(genotypes)
    assert res["D"] == pytest.approx(0.0 - 0.25)
    assert res["D_prime"] == pytest.approx(-1.0)
    assert res["r_squared"] == pytest.approx(1.0)


def test_ld_coefficients_genotype_mode_independent_loci():
    # Half 11 / half 00 with p1 = p2 = 0.5 gives D = 0
    genotypes = [[1, 1], [1, 0], [0, 1], [0, 0]]
    res = ld_coefficients(genotypes)
    assert res["D"] == pytest.approx(0.0, abs=1e-12)
    assert res["D_prime"] == 0  # D_max > 0 but D == 0
    assert res["r_squared"] == 0


def test_ld_coefficients_genotype_mode_validation():
    with pytest.raises(ValueError, match="at least 2"):
        ld_coefficients([[1, 0]])


def test_ld_coefficients_genotype_mode_positive_d_asymmetric():
    # Hand-computed, all four haplotypes present: 11:8, 10:4, 01:3, 00:5.
    # p1 = 0.6, p2 = 0.55; D = 0.4 - 0.33 = 0.07
    # D > 0 -> Dmax = min(p1*q2, q1*p2) = min(0.27, 0.22) = 0.22
    genotypes = [[1, 1]] * 8 + [[1, 0]] * 4 + [[0, 1]] * 3 + [[0, 0]] * 5
    res = ld_coefficients(genotypes)
    assert res["D"] == pytest.approx(0.07)
    assert res["D_prime"] == pytest.approx(0.07 / 0.22)
    assert res["r_squared"] == pytest.approx(0.0049 / 0.0594)


def test_ld_coefficients_genotype_mode_negative_d_asymmetric():
    # Hand-computed, all four haplotypes present: 11:1, 10:7, 01:4, 00:8.
    # p1 = 0.4, p2 = 0.25; D = 0.05 - 0.10 = -0.05
    # D < 0 -> Dmax = min(p1*p2, q1*q2) = min(0.10, 0.45) = 0.10
    genotypes = [[1, 1]] + [[1, 0]] * 7 + [[0, 1]] * 4 + [[0, 0]] * 8
    res = ld_coefficients(genotypes)
    assert res["D"] == pytest.approx(-0.05)
    assert res["D_prime"] == pytest.approx(-0.5)
    assert res["r_squared"] == pytest.approx(0.0025 / 0.045)


def test_ld_coefficients_unphased_diploid_matches_phased():
    # No double heterozygotes: EM must reproduce the phased result exactly.
    unphased = [[[1, 1], [1, 1]]] * 3 + [[[1, 0], [1, 1]]] * 2 + [[[0, 0], [0, 0]]] * 3
    phased = [[1, 1]] * 8 + [[0, 1]] * 2 + [[0, 0]] * 6
    res_u = ld_coefficients(unphased, phased=False)
    res_p = ld_coefficients(phased)
    assert res_u["D"] == pytest.approx(0.1875)
    assert res_u["D_prime"] == pytest.approx(1.0)
    assert res_u["r_squared"] == pytest.approx(0.6)
    for key in ("D", "D_prime", "r_squared"):
        assert res_u[key] == pytest.approx(res_p[key])


def test_ld_coefficients_unphased_double_heterozygote_em():
    # 3 diploid individuals plus one double heterozygote whose EM fixed
    # point is the complete-coupling boundary: haplotypes 11:3/8, 00:3/8,
    # 10:2/8, 01:0 -> D = 9/64, D' = 1, r^2 = 0.36.
    genotypes = [[[1, 1], [1, 1]], [[0, 0], [0, 0]], [[1, 0], [1, 0]], [[1, 1], [0, 0]]]
    res = ld_coefficients(genotypes, phased=False)
    assert res["D"] == pytest.approx(9 / 64, abs=1e-9)
    assert res["D_prime"] == pytest.approx(1.0, abs=1e-9)
    assert res["r_squared"] == pytest.approx(0.36, abs=1e-9)


def test_ld_coefficients_unphased_validation():
    with pytest.raises(ValueError, match="diploid"):
        ld_coefficients([[1, 1], [1, 0]], phased=False)
    with pytest.raises(ValueError, match="0/1"):
        ld_coefficients([[[1, 1], [1, 1]], [[1, 2], [1, 1]]], phased=False)
