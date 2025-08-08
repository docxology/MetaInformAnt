from __future__ import annotations

from typing import Iterable, Sequence, Tuple


def expectation(values: Iterable[float], weights: Iterable[float] | None = None) -> float:
    vs = list(values)
    if not vs:
        return 0.0
    if weights is None:
        return sum(vs) / len(vs)
    ws = list(weights)
    total_w = sum(ws)
    if total_w == 0:
        return 0.0
    return sum(v * w for v, w in zip(vs, ws)) / total_w


def covariance(x: Sequence[float], y: Sequence[float]) -> float:
    if not x or not y or len(x) != len(y):
        return 0.0
    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    return sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y)) / n


def price_equation(
    fitness: Sequence[float],
    trait_parent: Sequence[float],
    trait_offspring: Sequence[float] | None = None,
) -> Tuple[float, float, float]:
    """Compute a basic Price equation decomposition.

    Returns (covariance_term, transmission_term, total_change)

    - fitness: relative fitness (w)
    - trait_parent: parental trait values (z)
    - trait_offspring: offspring trait values (z') for transmission term; if None, term is 0.
    """
    if not fitness or len(fitness) != len(trait_parent):
        return 0.0, 0.0, 0.0
    w = list(fitness)
    z = list(trait_parent)
    cov_term = covariance(w, z)

    trans_term = 0.0
    if trait_offspring is not None and len(trait_offspring) == len(z):
        dz = [z_prime - zi for z_prime, zi in zip(trait_offspring, z)]
        trans_term = expectation(dz, weights=w)

    total = cov_term + trans_term
    return cov_term, trans_term, total


