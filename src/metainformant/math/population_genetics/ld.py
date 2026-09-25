"""Linkage disequilibrium functions.

This module provides mathematical functions for analyzing linkage disequilibrium (LD).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def ld_coefficients(
    pA_or_genotypes: float | List[List[int]],
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
        pA = pA_or_genotypes
        assert isinstance(pA, float)
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
    genotypes = pA_or_genotypes
    assert isinstance(genotypes, list)
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
        the four frequencies summing to 1.

    Raises:
        ValueError: If rows are not diploid two-locus pairs or alleles are
            not 0/1.
    """
    for row in genotypes:
        if (
            len(row) != 2
            or not isinstance(row[0], (list, tuple))
            or not isinstance(row[1], (list, tuple))
        ):
            raise ValueError(
                "Unphased diploid rows must be [[a1, a2], [b1, b2]]: "
                "the unordered allele pair at each of the two loci"
            )
        if len(row[0]) != 2 or len(row[1]) != 2:
            raise ValueError(
                "Unphased diploid rows must hold exactly two alleles per locus"
            )
        for allele in (row[0][0], row[0][1], row[1][0], row[1][1]):
            if allele not in (0, 1):
                raise ValueError("Unphased diploid mode requires 0/1 allele coding")

    n = len(genotypes)
    locus1 = [row[0] for row in genotypes]
    locus2 = [row[1] for row in genotypes]
    p1 = sum(sum(pair) for pair in locus1) / (2 * n)
    p2 = sum(sum(pair) for pair in locus2) / (2 * n)

    # Deterministic product-of-marginals initialization.
    freqs = {
        (a, b): (p1 if a else 1 - p1) * (p2 if b else 1 - p2)
        for a in (0, 1)
        for b in (0, 1)
    }

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


def ld_decay_r2(
    distances: List[float] | float,
    r_squared_values: List[float] | None = None,
    max_distance: Optional[float] = None,
    *,
    recombination_rate: float | None = None,
    generations: int | None = None,
) -> Dict[str, float] | float:
    """Analyze LD decay with distance.

    Args:
        distances: Physical distances between SNP pairs
        r_squared_values: Corresponding r² values
        max_distance: Maximum distance to consider

    Returns:
        Dictionary with LD decay statistics
    """
    if isinstance(distances, (int, float)):
        if recombination_rate is None or generations is None:
            raise ValueError(
                "recombination_rate and generations are required for scalar LD decay"
            )
        return float(distances) * ((1.0 - recombination_rate) ** (2 * generations))

    if r_squared_values is None or len(distances) != len(r_squared_values):
        raise ValueError("Distances and r² values must have same length")

    # Filter by max distance if specified
    if max_distance is not None:
        filtered = [
            (d, r) for d, r in zip(distances, r_squared_values) if d <= max_distance
        ]
        distances = [d for d, _ in filtered]
        r_squared_values = [r for _, r in filtered]

    if not distances:
        return {"decay_rate": 0.0, "half_decay_distance": float("inf")}

    # Fit exponential decay: r² = r²₀ * exp(-d/d₀)
    # Use simple binning approach
    bins = np.logspace(0, np.log10(max(distances)), 20)
    bin_means = []

    for i in range(len(bins) - 1):
        mask = (np.array(distances) >= bins[i]) & (np.array(distances) < bins[i + 1])
        if np.any(mask):
            mean_r2 = np.mean(np.array(r_squared_values)[mask])
            bin_means.append((bins[i], mean_r2))

    if len(bin_means) < 3:
        return {"decay_rate": 0.0, "half_decay_distance": float("inf")}

    # Estimate half-decay distance (where r² drops to 0.5 of initial)
    initial_r2 = bin_means[0][1]
    half_decay_value = initial_r2 * 0.5

    half_decay_distance = float("inf")
    for dist, r2 in bin_means:
        if r2 <= half_decay_value:
            half_decay_distance = dist
            break

    # Estimate decay rate (rough approximation)
    if len(bin_means) >= 2:
        decay_rate = (bin_means[0][1] - bin_means[-1][1]) / (
            bin_means[-1][0] - bin_means[0][0]
        )
    else:
        decay_rate = 0.0

    return {
        "decay_rate": decay_rate,
        "half_decay_distance": half_decay_distance,
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
        genetic_distance: Genetic distance in Morgans

    Returns:
        Recombination fraction (0 <= c <= 0.5)
    """
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

    return 0.25 * math.log(
        (1 + 2 * recombination_fraction) / (1 - 2 * recombination_fraction)
    )


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
