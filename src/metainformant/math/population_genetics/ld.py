"""Linkage disequilibrium functions.

This module provides mathematical functions for analyzing linkage disequilibrium (LD).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def ld_coefficients(
    pA_or_genotypes: float | Sequence[Any] | np.ndarray,
    pa: float | None = None,
    pB: float | None = None,
    pb: float | None = None,
    pAB: float | None = None,
    phased: bool = True,
) -> Dict[str, float] | Tuple[float, float]:
    """Calculate linkage disequilibrium coefficients.

    Can be called with either:
    - Allele frequencies: ld_coefficients(pA, pa, pB, pb, pAB) -> (D, D')
    - Phased haplotypes: ld_coefficients(genotypes) -> dict with D, D', r²,
      where each row ``[g1, g2]`` is a single PHASED haplotype (allele g1 at
      locus 1, allele g2 at locus 2)
    - Unphased diploid genotypes: ld_coefficients(genotypes, phased=False)
      -> dict with D, D', r², where each row is ``[[a1, a2], [b1, b2]]``:
      the unordered allele pair at locus 1 and at locus 2.
      Maximum-likelihood haplotype frequencies are recovered with a
      deterministic expectation-maximization (EM) routine.

    Args:
        pA_or_genotypes: Either allele frequency for A or 2D list of genotypes
        pa: Allele frequency for a (if using frequency mode)
        pB: Allele frequency for B (if using frequency mode)
        pb: Allele frequency for b (if using frequency mode)
        pAB: Haplotype frequency for AB (if using frequency mode)
        phased: Genotype-mode interpretation. ``True`` (default) treats each
            row as one phased haplotype; ``False`` treats each row as an
            unphased diploid genotype and runs EM.

    Returns:
        If using frequencies: Tuple of (D, D')
        If using genotypes: Dictionary with LD coefficients (D, D_prime, r_squared)
    """
    # Check if we're in frequency mode
    if pa is not None and pB is not None and pb is not None and pAB is not None:
        if isinstance(pA_or_genotypes, bool) or not isinstance(pA_or_genotypes, (int, float, np.integer, np.floating)):
            raise TypeError(
                f"The first argument must be the frequency of allele A in frequency mode, "
                f"got {type(pA_or_genotypes).__name__}"
            )
        pA = float(pA_or_genotypes)
        # Calculate D = P(AB) - P(A)*P(B)
        D = pAB - pA * pB

        # Calculate D' = D / D_max
        if D >= 0:
            D_max = min(pA * pb, pa * pB)
        else:
            D_max = min(pA * pB, pa * pb)

        D_prime = D / D_max if D_max > 0 else 0

        return D, D_prime

    # Otherwise use genotype mode
    if not isinstance(pA_or_genotypes, (list, tuple, np.ndarray)):
        raise TypeError(
            f"genotypes must be a list, tuple or array of rows in genotype mode, got {type(pA_or_genotypes).__name__}"
        )
    genotypes: List[Any] = list(
        pA_or_genotypes.tolist() if isinstance(pA_or_genotypes, np.ndarray) else pA_or_genotypes
    )
    if len(genotypes) < 2 or len(genotypes[0]) != 2:
        raise ValueError("Need at least 2 samples with 2 loci each")

    n = len(genotypes)
    if phased:
        # Each row is one phased haplotype: [allele at locus 1, allele at locus 2].
        for row in genotypes:
            if len(row) != 2:
                raise ValueError("Phased haplotype rows must be [g1, g2]")
        haplotype_freqs: Dict[str, float] = {}
        for g1, g2 in genotypes:
            key = f"{g1}{g2}"
            haplotype_freqs[key] = haplotype_freqs.get(key, 0) + 1
        f_11 = haplotype_freqs.get("11", 0) / n
        p1 = sum(row[0] for row in genotypes) / n
        p2 = sum(row[1] for row in genotypes) / n
    else:
        # Unphased diploid rows: [[a1, a2], [b1, b2]]. Haplotype frequencies
        # are recovered by EM before computing the coefficients.
        freqs = _em_haplotype_frequencies(genotypes)
        f_11 = freqs[(1, 1)]
        p1 = freqs[(1, 0)] + freqs[(1, 1)]
        p2 = freqs[(0, 1)] + freqs[(1, 1)]

    q1 = 1 - p1
    q2 = 1 - p2

    # D = P(AB) - P(A)P(B)
    D = f_11 - p1 * p2

    # D' = D / D_max. For D > 0 the excess of AB is bounded by the scarcer
    # of the Ab / aB complements; for D < 0 the AB deficit is bounded by
    # the scarcer of AB / ab themselves.
    D_max = min(p1 * q2, q1 * p2) if D > 0 else min(p1 * p2, q1 * q2)
    D_prime = D / D_max if D_max > 0 else 0

    # r² = D² / (p1*q1*p2*q2)
    r_squared = (D**2) / (p1 * q1 * p2 * q2) if (p1 * q1 * p2 * q2) > 0 else 0

    return {"D": D, "D_prime": D_prime, "r_squared": r_squared}


def _em_haplotype_frequencies(
    genotypes: List[List[Any]],
    max_iter: int = 500,
    tol: float = 1e-12,
) -> Dict[Tuple[int, int], float]:
    """Maximum-likelihood haplotype frequencies from unphased diploid genotypes.

    Expectation-maximization for two biallelic loci. Each genotype row is
    ``[[a1, a2], [b1, b2]]``: the unordered allele pair at locus 1 and the
    unordered allele pair at locus 2. Only double heterozygotes
    (heterozygous at both loci) are ambiguous; the split between the
    coupling and repulsion configurations is iterated from a deterministic
    product-of-marginals initialization until convergence.

    Args:
        genotypes: Unphased diploid two-locus genotypes with 0/1 alleles.
        max_iter: EM iteration cap.
        tol: Convergence threshold on the maximum frequency change.

    Returns:
        Mapping (allele at locus 1, allele at locus 2) -> frequency, with
        the four frequencies summing to 1. The iteration stops after ``max_iter`` rounds even if the
        frequencies are still changing by more than ``tol``; no warning is given.

    Raises:
        ValueError: If rows are not diploid two-locus pairs or alleles are
            not 0/1.
    """
    for row in genotypes:
        if len(row) != 2 or not isinstance(row[0], (list, tuple)) or not isinstance(row[1], (list, tuple)):
            raise ValueError(
                "Unphased diploid rows must be [[a1, a2], [b1, b2]]: "
                "the unordered allele pair at each of the two loci"
            )
        if len(row[0]) != 2 or len(row[1]) != 2:
            raise ValueError("Unphased diploid rows must hold exactly two alleles per locus")
        for allele in (row[0][0], row[0][1], row[1][0], row[1][1]):
            if allele not in (0, 1):
                raise ValueError("Unphased diploid mode requires 0/1 allele coding")

    n = len(genotypes)
    locus1 = [row[0] for row in genotypes]
    locus2 = [row[1] for row in genotypes]
    p1 = sum(sum(pair) for pair in locus1) / (2 * n)
    p2 = sum(sum(pair) for pair in locus2) / (2 * n)

    # Deterministic product-of-marginals initialization.
    freqs = {(a, b): (p1 if a else 1 - p1) * (p2 if b else 1 - p2) for a in (0, 1) for b in (0, 1)}

    for _ in range(max_iter):
        expected = dict.fromkeys(freqs, 0.0)
        for (l1a, l1b), (l2a, l2b) in zip(locus1, locus2):
            if l1a == l1b:
                if l2a == l2b:
                    expected[(l1a, l2a)] += 2.0
                else:
                    expected[(l1a, l2a)] += 1.0
                    expected[(l1a, l2b)] += 1.0
            elif l2a == l2b:
                expected[(l1a, l2a)] += 1.0
                expected[(l1b, l2a)] += 1.0
            else:
                # Double heterozygote: split coupling vs repulsion by the
                # current haplotype frequencies.
                coupling = freqs[(l1a, l2a)] * freqs[(l1b, l2b)]
                repulsion = freqs[(l1a, l2b)] * freqs[(l1b, l2a)]
                total = coupling + repulsion
                weight = coupling / total if total > 0 else 0.5
                expected[(l1a, l2a)] += weight
                expected[(l1b, l2b)] += weight
                expected[(l1a, l2b)] += 1.0 - weight
                expected[(l1b, l2a)] += 1.0 - weight

        updated = {key: value / (2 * n) for key, value in expected.items()}
        delta = max(abs(updated[key] - freqs[key]) for key in freqs)
        freqs = updated
        if delta < tol:
            break

    return freqs


def _fit_exponential_decay(distances: np.ndarray, r_squared: np.ndarray) -> Tuple[float, float]:
    """Least-squares fit of r² = a * exp(-k * d) with a >= 0 and k >= 0.

    Every (distance, r²) pair enters the fit with equal weight; nothing is binned.

    Args:
        distances: Distances between SNP pairs (non-negative, at least three distinct values)
        r_squared: Corresponding r² values

    Returns:
        (a, k): the fitted r² at distance 0 and the decay rate per unit of distance

    Raises:
        RuntimeError: If the optimiser does not converge
    """
    from scipy.optimize import least_squares

    # Work in units of the longest distance so the optimiser is well scaled whatever the unit.
    scale = float(distances.max())
    x = distances / scale

    # Start from the mean r² of the closest third of the pairs and the slope of a straight line
    # through log(r²); the bounded fit below refines both.
    closest = np.argsort(x)[: max(1, len(x) // 3)]
    a0 = max(float(np.mean(r_squared[closest])), 1e-6)
    positive = r_squared > 0
    if np.count_nonzero(positive) >= 2 and np.ptp(x[positive]) > 0:
        k0 = max(-float(np.polyfit(x[positive], np.log(r_squared[positive]), 1)[0]), 1e-3)
    else:
        k0 = 1.0

    def residuals(params: np.ndarray) -> np.ndarray:
        return np.asarray(params[0] * np.exp(-params[1] * x) - r_squared, dtype=float)

    result = least_squares(residuals, x0=[a0, k0], bounds=([0.0, 0.0], [np.inf, np.inf]))
    if not result.success:
        raise RuntimeError(f"Exponential LD decay fit did not converge: {result.message}")
    return float(result.x[0]), float(result.x[1]) / scale


def ld_decay_r2(
    distances: List[float] | float,
    r_squared_values: List[float] | None = None,
    max_distance: Optional[float] = None,
    *,
    recombination_rate: float | None = None,
    generations: int | None = None,
) -> Dict[str, float] | float:
    """Analyze LD decay with distance.

    Fits the exponential model r² = r²₀ * exp(-k * d) to the (distance, r²) pairs by bounded
    nonlinear least squares (r²₀ >= 0, k >= 0), using every pair with equal weight. It is a
    descriptive decay constant: the model has no finite-sample floor (a sample of n haplotypes
    has E[r²] near 1/n at long range) and it is not the Hill and Weir (1988) expectation, so k
    does not estimate 4Nₑc.

    With a scalar ``distances`` it instead returns the expected r² after ``generations`` of
    recombination at rate ``recombination_rate``: ``distances * (1 - c) ** (2 * generations)``.

    Args:
        distances: Distances between SNP pairs (any unit; the rate comes back per that unit), or
            a scalar starting r² for the closed-form decay
        r_squared_values: Corresponding r² values
        max_distance: Maximum distance to consider
        recombination_rate: Recombination fraction per generation (scalar form only)
        generations: Number of generations (scalar form only)

    Returns:
        A dictionary for the fit, or a float for the scalar form. The fit gives

        - ``decay_rate``: k, the exponential rate constant per unit of distance
        - ``half_decay_distance``: ln(2) / k, the distance at which the fitted r² is half its
          value at distance 0 (infinite when there is no decay)
        - ``initial_r2``: r²₀, the fitted r² at distance 0

        When no pair is left or fewer than three distinct distances remain, there is nothing to
        fit and the dictionary holds ``decay_rate`` 0.0 and an infinite ``half_decay_distance``
        only.

    Raises:
        ValueError: If the inputs differ in length, hold a non-finite value, or a distance is negative
        RuntimeError: If the optimiser does not converge
    """
    if isinstance(distances, (int, float)):
        if recombination_rate is None or generations is None:
            raise ValueError("recombination_rate and generations are required for scalar LD decay")
        return float(distances) * ((1.0 - recombination_rate) ** (2 * generations))

    if r_squared_values is None or len(distances) != len(r_squared_values):
        raise ValueError("Distances and r² values must have same length")

    d = np.asarray(distances, dtype=float)
    r2 = np.asarray(r_squared_values, dtype=float)
    if not (np.all(np.isfinite(d)) and np.all(np.isfinite(r2))):
        raise ValueError("Distances and r² values must be finite")
    if np.any(d < 0):
        raise ValueError("Distances must be non-negative")

    # Filter by max distance if specified
    if max_distance is not None:
        keep = d <= max_distance
        d, r2 = d[keep], r2[keep]

    # Two parameters need at least three distinct distances to be fitted with any redundancy.
    if np.unique(d).size < 3:
        return {"decay_rate": 0.0, "half_decay_distance": float("inf")}

    if not np.any(r2 > 0):
        return {"decay_rate": 0.0, "half_decay_distance": float("inf"), "initial_r2": 0.0}

    initial_r2, decay_rate = _fit_exponential_decay(d, r2)
    # The model contains the constant curve (k = 0), so it never fits worse than the mean. When it
    # fits no better, the data show no decay (flat or rising): report none, not a vanishing rate.
    flat_sse = float(np.sum((r2 - r2.mean()) ** 2))
    fit_sse = float(np.sum((initial_r2 * np.exp(-decay_rate * d) - r2) ** 2))
    if fit_sse >= flat_sse * (1.0 - 1e-9):
        decay_rate, initial_r2 = 0.0, float(r2.mean())

    return {
        "decay_rate": decay_rate,
        "half_decay_distance": math.log(2.0) / decay_rate if decay_rate > 0 else float("inf"),
        "initial_r2": initial_r2,
    }


def haldane_c_to_d(recombination_fraction: float) -> float:
    """Convert recombination fraction to genetic distance using Haldane's mapping function.

    Args:
        recombination_fraction: Recombination fraction (c) between 0 and 0.5

    Returns:
        Genetic distance in Morgans (d)
    """
    if not (0 <= recombination_fraction <= 0.5):
        raise ValueError("Recombination fraction must be between 0 and 0.5")

    if recombination_fraction == 0:
        return 0.0
    elif recombination_fraction == 0.5:
        return float("inf")

    # Haldane's mapping function: d = -0.5 * ln(1 - 2c)
    return float(-0.5 * np.log(1 - 2 * recombination_fraction))


def haldane_d_to_c(genetic_distance: float) -> float:
    """Convert genetic distance (Morgans) to recombination fraction using Haldane's mapping function.

    Haldane's mapping function: c = 0.5(1 - exp(-2d))
    where d is genetic distance in Morgans and c is recombination fraction.

    Args:
        genetic_distance: Genetic distance in Morgans, non-negative

    Returns:
        Recombination fraction (0 <= c <= 0.5)

    Raises:
        ValueError: If genetic_distance is negative
    """
    if genetic_distance < 0:
        raise ValueError("Genetic distance cannot be negative")

    return 0.5 * (1 - math.exp(-2 * genetic_distance))


def kosambi_c_to_d(recombination_fraction: float) -> float:
    """Convert recombination fraction to genetic distance using Kosambi mapping function.

    The Kosambi mapping function: c = 0.5 * tanh(2d)
    where c is recombination fraction and d is genetic distance in Morgans.

    Args:
        recombination_fraction: Recombination fraction (0 <= c <= 0.5)

    Returns:
        Genetic distance in Morgans
    """
    if not (0 <= recombination_fraction <= 0.5):
        raise ValueError("Recombination fraction must be between 0 and 0.5")

    if recombination_fraction == 0.5:
        return float("inf")  # Infinite distance
    elif recombination_fraction == 0:
        return 0.0

    return 0.25 * math.log((1 + 2 * recombination_fraction) / (1 - 2 * recombination_fraction))


def kosambi_d_to_c(genetic_distance: float) -> float:
    """Convert genetic distance to recombination fraction using Kosambi mapping function.

    The Kosambi mapping function: c = 0.5 * tanh(2d)
    where d is genetic distance in Morgans and c is recombination fraction.

    Args:
        genetic_distance: Genetic distance in Morgans

    Returns:
        Recombination fraction (0 <= c <= 0.5)
    """
    if genetic_distance < 0:
        raise ValueError("Genetic distance cannot be negative")

    if genetic_distance == 0:
        return 0.0

    return 0.5 * math.tanh(2 * genetic_distance)
