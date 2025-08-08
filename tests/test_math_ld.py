from __future__ import annotations

from metainformant.math import ld_coefficients, r_squared


def test_ld_coefficients_and_r2():
    # pA=0.6, pa=0.4, pB=0.7, pb=0.3, pAB=0.5
    D, Dp = ld_coefficients(0.6, 0.4, 0.7, 0.3, 0.5)
    assert abs(D - (0.5 - 0.42)) < 1e-12
    # Dmax = min(pA*pb, pa*pB) = min(0.18, 0.28) = 0.18 => D' = 0.08/0.18
    assert abs(Dp - ((0.08) / 0.18)) < 1e-12

    r2 = r_squared(0.6, 0.4, 0.7, 0.3, 0.5)
    denom = 0.6 * 0.4 * 0.7 * 0.3
    assert abs(r2 - ((0.08 * 0.08) / denom)) < 1e-12



