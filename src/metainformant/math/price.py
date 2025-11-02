from __future__ import annotations

from typing import Iterable, Sequence, Tuple


def expectation(values: Iterable[float], weights: Iterable[float] | None = None) -> float:
    """Calculate weighted or unweighted expectation (mean) of values.
    
    Args:
        values: Iterable of numeric values
        weights: Optional weights for each value. If None, computes simple mean.
        
    Returns:
        Weighted mean if weights provided, else arithmetic mean. Returns 0.0 for empty input.
        
    Examples:
        >>> expectation([1.0, 2.0, 3.0])
        2.0
        >>> expectation([1.0, 2.0, 3.0], weights=[1.0, 2.0, 1.0])
        2.0
    """
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
    """Calculate population covariance between two sequences.
    
    Uses population formula: Cov(X,Y) = E[(X - E[X])(Y - E[Y])] = E[XY] - E[X]E[Y]
    
    Args:
        x: First sequence of values
        y: Second sequence of values (must match length of x)
        
    Returns:
        Population covariance. Returns 0.0 if inputs are empty, mismatched lengths, or invalid.
        
    Examples:
        >>> covariance([1, 2, 3], [2, 4, 6])
        1.333...
        >>> covariance([1, 1, 1], [2, 3, 4])
        0.0
    """
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
    """Normalize fitness values by their mean to obtain relative fitness.
    
    Relative fitness is used in evolutionary analysis where the mean fitness
    defines the baseline. After normalization, mean relative fitness is 1.0.
    
    Args:
        fitness: Sequence of absolute fitness values
        
    Returns:
        List of relative fitness values (fitness / mean(fitness)). 
        Returns zeros if mean fitness is 0 or input is empty.
        The mean of returned values is 1.0 when mean(fitness) > 0.
        
    Examples:
        >>> relative_fitness([1.0, 2.0, 3.0])
        [0.5, 1.0, 1.5]
        >>> sum(relative_fitness([1.0, 2.0, 3.0])) / 3
        1.0
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
    """Calculate selection differential S = Cov(w_rel, z).
    
    The selection differential measures the change in mean trait due to selection
    and is the covariance between relative fitness and trait values.
    
    Args:
        fitness: Sequence of fitness values
        trait: Sequence of trait values (must match length of fitness)
        normalize_by_mean_fitness: If True (default), normalize fitness by mean.
            If False, use raw fitness values (assumes they are already relative).
            
    Returns:
        Selection differential. Returns 0.0 if inputs are empty or mismatched.
        
    Examples:
        >>> selection_differential([1.0, 1.2, 0.9], [0.2, 0.4, 0.1])
        0.0333...
        
    References:
        Lande, R., & Arnold, S. J. (1983). The measurement of selection on
        correlated characters. Evolution, 37(6), 1210-1226.
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
    """Calculate Lande–Arnold directional selection gradient β = Cov(w_rel, z) / Var(z).
    
    The selection gradient measures the strength of directional selection per unit
    of phenotypic variance. Standardized by trait variance to remove scale dependence.
    
    Args:
        fitness: Sequence of fitness values
        trait: Sequence of trait values (must match length of fitness)
        normalize_by_mean_fitness: If True (default), normalize fitness by mean.
            
    Returns:
        Selection gradient. Returns 0.0 if variance is 0 or inputs are invalid.
        
    Examples:
        >>> selection_gradient([1.0, 1.2, 0.9], [0.2, 0.4, 0.1])
        0.25...
        
    References:
        Lande, R., & Arnold, S. J. (1983). The measurement of selection on
        correlated characters. Evolution, 37(6), 1210-1226.
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
    """Compute Price equation decomposition of evolutionary change.
    
    The Price equation decomposes change in mean trait into selection (covariance)
    and transmission components: Δz̄ = Cov(w, z) + E(wΔz)
    
    Args:
        fitness: Sequence of fitness values (w)
        trait_parent: Sequence of parental trait values (z)
        trait_offspring: Optional sequence of offspring trait values (z').
            If None, transmission term is set to 0.
            
    Returns:
        Tuple of (covariance_term, transmission_term, total_change):
        - covariance_term: Selection component Cov(w_rel, z)
        - transmission_term: Transmission bias E(w_rel * (z' - z))
        - total_change: Total change in mean trait
        
    Examples:
        >>> fitness = [1.0, 1.2, 0.9]
        >>> parent = [0.2, 0.4, 0.1]
        >>> offspring = [0.25, 0.35, 0.15]
        >>> cov, trans, total = price_equation(fitness, parent, offspring)
        >>> abs(total - (cov + trans)) < 1e-10
        True
        
    References:
        Price, G. R. (1970). Selection and covariance. Nature, 227(5257), 520-521.
        Frank, S. A. (2012). Natural selection. IV. The Price equation.
        Journal of Evolutionary Biology, 25(5), 1002-1019.
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
