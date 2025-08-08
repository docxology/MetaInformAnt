from __future__ import annotations

from typing import Tuple


def ld_coefficients(
    pA: float, pa: float, pB: float, pb: float, haplotype_pAB: float
) -> Tuple[float, float]:
    """Compute linkage disequilibrium D and normalized D' given haplotype frequencies.

    Inputs should satisfy pA + pa = 1, pB + pb = 1, and 0 <= haplotype_pAB <= 1.
    Returns (D, D_prime). If impossible bounds, returns (0.0, 0.0).
    """
    if not (0 <= pA <= 1 and 0 <= pa <= 1 and abs(pA + pa - 1.0) < 1e-8):
        return 0.0, 0.0
    if not (0 <= pB <= 1 and 0 <= pb <= 1 and abs(pB + pb - 1.0) < 1e-8):
        return 0.0, 0.0
    if not (0.0 <= haplotype_pAB <= 1.0):
        return 0.0, 0.0
    D = haplotype_pAB - pA * pB
    Dmax = min(pA * pb, pa * pB) if D >= 0 else min(pA * pB, pa * pb)
    if Dmax <= 0:
        return D, 0.0
    return D, D / Dmax


def r_squared(pA: float, pa: float, pB: float, pb: float, haplotype_pAB: float) -> float:
    """Compute r^2 for two biallelic loci using haplotype frequencies."""
    D, _ = ld_coefficients(pA, pa, pB, pb, haplotype_pAB)
    denom = pA * pa * pB * pb
    if denom <= 0:
        return 0.0
    return (D * D) / denom


def ld_decay_r2(r2_initial: float, recombination_rate: float, generations: int) -> float:
    """Approximate LD decay under random mating: r2_t ≈ r2_0 * (1 - c)^{2t}.

    c is recombination rate per generation. Clamps to [0, 1].
    """
    r2 = max(0.0, min(1.0, r2_initial))
    c = max(0.0, min(1.0, recombination_rate))
    t = max(0, generations)
    return r2 * ((1.0 - c) ** (2 * t))



