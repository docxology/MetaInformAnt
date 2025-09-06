from __future__ import annotations

from typing import Tuple


def ld_coefficients(pA: float, pa: float, pB: float, pb: float, haplotype_pAB: float) -> Tuple[float, float]:
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


def haldane_d_to_c(map_distance_morgans: float) -> float:
    """Haldane mapping function: c = 0.5 (1 - exp(-2d))."""
    d = max(0.0, float(map_distance_morgans))
    from math import exp

    c = 0.5 * (1.0 - exp(-2.0 * d))
    return max(0.0, min(0.5, c))


def haldane_c_to_d(recombination_fraction: float) -> float:
    """Inverse Haldane: d = -0.5 ln(1 - 2c)."""
    c = max(0.0, min(0.5, float(recombination_fraction)))
    from math import log

    if c >= 0.5:
        return float("inf")
    return -0.5 * log(1.0 - 2.0 * c)


def kosambi_d_to_c(map_distance_morgans: float) -> float:
    """Kosambi mapping: c = 0.5 * tanh(2d)."""
    d = max(0.0, float(map_distance_morgans))
    from math import tanh

    c = 0.5 * tanh(2.0 * d)
    return max(0.0, min(0.5, c))


def kosambi_c_to_d(recombination_fraction: float) -> float:
    """Inverse Kosambi: d = 0.25 ln((1 + 2c) / (1 - 2c))."""
    c = max(0.0, min(0.5, float(recombination_fraction)))
    from math import log

    if c >= 0.5:
        return float("inf")
    return 0.25 * log((1.0 + 2.0 * c) / (1.0 - 2.0 * c))


def expected_r2_from_Ne_c(effective_population_size: float, recombination_fraction: float) -> float:
    """Approximate E[r^2] ≈ 1 / (1 + 4 Ne c) (neutral, unlinked approximation)."""
    Ne = max(0.0, float(effective_population_size))
    c = max(0.0, min(0.5, float(recombination_fraction)))
    denom = 1.0 + 4.0 * Ne * c
    if denom <= 0:
        return 0.0
    return 1.0 / denom
