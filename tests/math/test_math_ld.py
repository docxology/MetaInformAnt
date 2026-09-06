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
