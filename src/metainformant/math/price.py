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


def variance(values: Sequence[float]) -> float:
    """Population variance of a sequence (dividing by n).

    Returns 0.0 for empty input.
    """
    if not values:
        return 0.0
    n = len(values)
    mu = sum(values) / n
    return sum((v - mu) ** 2 for v in values) / n


def correlation(x: Sequence[float], y: Sequence[float]) -> float:
    """Pearson correlation coefficient between two sequences.

    Returns 0.0 if inputs are empty, different lengths, or either variance is 0.
    Uses population definitions (divide by n).
    """
    cov = covariance(x, y)
    var_x = variance(x)
    var_y = variance(y)
    if var_x <= 0.0 or var_y <= 0.0:
        return 0.0
    return cov / ((var_x**0.5) * (var_y**0.5))


def standard_deviation(values: Sequence[float]) -> float:
    """Population standard deviation (sqrt of population variance)."""
    return variance(values) ** 0.5


def weighted_variance(values: Sequence[float], weights: Sequence[float]) -> float:
    """Weighted population variance with weights normalized to sum to 1.

    If total weight is 0 or inputs are empty/mismatched, returns 0.0.
    """
    if not values or not weights or len(values) != len(weights):
        return 0.0
    w_sum = float(sum(weights))
    if w_sum == 0.0:
        return 0.0
    w_norm = [w / w_sum for w in weights]
    mu = sum(v * w for v, w in zip(values, w_norm))
    return sum(w * (v - mu) ** 2 for v, w in zip(values, w_norm))


def weighted_covariance(x: Sequence[float], y: Sequence[float], weights: Sequence[float]) -> float:
    """Weighted population covariance with weights normalized to sum to 1.

    If total weight is 0 or inputs are empty/mismatched, returns 0.0.
    """
    if not x or not y or not weights or len(x) != len(y) or len(x) != len(weights):
        return 0.0
    w_sum = float(sum(weights))
    if w_sum == 0.0:
        return 0.0
    w_norm = [w / w_sum for w in weights]
    mean_x = sum(xi * wi for xi, wi in zip(x, w_norm))
    mean_y = sum(yi * wi for yi, wi in zip(y, w_norm))
    return sum(wi * (xi - mean_x) * (yi - mean_y) for xi, yi, wi in zip(x, y, w_norm))


def weighted_correlation(x: Sequence[float], y: Sequence[float], weights: Sequence[float]) -> float:
    """Weighted Pearson correlation based on weighted population moments.

    Returns 0.0 if either weighted variance is 0 or inputs invalid.
    """
    cov_w = weighted_covariance(x, y, weights)
    var_x = weighted_variance(x, weights)
    var_y = weighted_variance(y, weights)
    if var_x <= 0.0 or var_y <= 0.0:
        return 0.0
    return cov_w / ((var_x**0.5) * (var_y**0.5))


def relative_fitness(fitness: Sequence[float]) -> list[float]:
    """Return fitness normalized by its mean. If mean is 0, returns zeros.

    The mean of the returned values is 1.0 when mean(fitness) > 0.
    """
    if not fitness:
        return []
    mean_w = sum(fitness) / len(fitness)
    if mean_w == 0.0:
        return [0.0 for _ in fitness]
    return [wi / mean_w for wi in fitness]


def selection_differential(
    fitness: Sequence[float],
    trait: Sequence[float],
    normalize_by_mean_fitness: bool = True,
) -> float:
    """Selection differential S = Cov(w_rel, z).

    If ``normalize_by_mean_fitness`` is True (default), relative fitness is
    defined as w_rel = w / mean(w). If False, uses raw w (appropriate when w
    are already relative and mean(w) ≈ 1).
    """
    if not fitness or len(fitness) != len(trait):
        return 0.0
    w = list(fitness)
    z = list(trait)
    if normalize_by_mean_fitness:
        mean_w = sum(w) / len(w)
        if mean_w == 0.0:
            return 0.0
        w_rel = [wi / mean_w for wi in w]
    else:
        w_rel = w
    return covariance(w_rel, z)


def selection_gradient(
    fitness: Sequence[float],
    trait: Sequence[float],
    normalize_by_mean_fitness: bool = True,
) -> float:
    """Lande–Arnold directional selection gradient β = Cov(w_rel, z) / Var(z).

    Uses population moments (divide by n). Returns 0.0 if variance is 0.
    """
    var_z = variance(trait)
    if var_z == 0.0:
        return 0.0
    S = selection_differential(fitness, trait, normalize_by_mean_fitness)
    return S / var_z


def selection_intensity(
    fitness: Sequence[float],
    trait: Sequence[float],
    normalize_by_mean_fitness: bool = True,
) -> float:
    """Standardized selection differential i = S / sd(z).

    Returns 0.0 if sd(z) == 0.
    """
    sd_z = standard_deviation(trait)
    if sd_z == 0.0:
        return 0.0
    S = selection_differential(fitness, trait, normalize_by_mean_fitness)
    return S / sd_z


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

    w_raw = list(fitness)
    z = list(trait_parent)

    # Use standard Price form with relative fitness normalized by mean(w)
    mean_w = sum(w_raw) / len(w_raw)
    if mean_w == 0.0:
        return 0.0, 0.0, 0.0
    w_rel = [wi / mean_w for wi in w_raw]

    cov_term = covariance(w_rel, z)

    trans_term = 0.0
    if trait_offspring is not None and len(trait_offspring) == len(z):
        dz = [z_prime - zi for z_prime, zi in zip(trait_offspring, z)]
        trans_term = expectation(dz, weights=w_rel)

    total = cov_term + trans_term
    return cov_term, trans_term, total


def delta_mean_trait(
    fitness: Sequence[float],
    trait_parent: Sequence[float],
    trait_offspring: Sequence[float] | None = None,
) -> float:
    """Change in mean trait Δz̄ given fitness and trait values.

    Equivalent to the total term from ``price_equation``.
    If ``trait_offspring`` is None, returns Cov(w_rel, z).
    """
    _, _, total = price_equation(fitness, trait_parent, trait_offspring)
    return total
