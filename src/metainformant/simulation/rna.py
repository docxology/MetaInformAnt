from __future__ import annotations

import math
import random


def _sample_negative_binomial(rng: random.Random, mean: float, dispersion: float) -> int:
    """Sample from NB parameterized by mean and dispersion (size = 1/dispersion).

    Uses Poisson-Gamma mixture.
    """
    if mean <= 0:
        return 0
    if dispersion <= 0:
        # approx Poisson
        lam = mean
        # Knuth algorithm for Poisson
        L = math.exp(-lam)
        k = 0
        p = 1.0
        while p > L:
            k += 1
            p *= rng.random()
        return k - 1
    size = 1.0 / dispersion
    p = size / (size + mean)
    # gamma-poisson mixture
    gamma_shape = size
    gamma_scale = (1 - p) / p
    # sample gamma via Marsaglia and Tsang (k >= 1) with fallback
    if gamma_shape < 1:
        # boost shape
        u = rng.random()
        return int(_sample_negative_binomial(rng, mean * u ** (1.0 / gamma_shape), dispersion))
    d = gamma_shape - 1.0 / 3.0
    c = 1.0 / math.sqrt(9 * d)
    while True:
        x = rng.normal(0, 1)
        v = (1 + c * x) ** 3
        if v <= 0:
            continue
        u = rng.random()
        if u < 1 - 0.0331 * (x**4) or math.log(u) < 0.5 * x * x + d * (1 - v + math.log(v)):
            lam = d * v * gamma_scale
            # Poisson draw with rate lam
            L = math.exp(-lam)
            k = 0
            p2 = 1.0
            while p2 > L:
                k += 1
                p2 *= rng.random()
            return k - 1


def simulate_counts_negative_binomial(
    num_genes: int,
    num_samples: int,
    *,
    mean_expression: float = 100.0,
    dispersion: float = 0.1,
    rng: random.Random | None = None,
) -> list[list[int]]:
    """Simulate RNA-seq count matrix (genes x samples) from NB.

    Returns a list of rows, each a list of counts per sample.
    """
    r = rng or random
    matrix: list[list[int]] = []
    for _ in range(max(0, num_genes)):
        row = [_sample_negative_binomial(r, mean_expression, dispersion) for _ in range(max(0, num_samples))]
        matrix.append(row)
    return matrix
