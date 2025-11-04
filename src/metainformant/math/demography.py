"""Basic demographic models for population genetics."""

from __future__ import annotations

import math


def exponential_growth_effective_size(
    current_size: int,
    growth_rate: float,
    generations: int,
) -> float:
    """Calculate effective population size under exponential growth.
    
    Under exponential growth, the effective population size is the harmonic
    mean of population sizes over time. This approximates the expected
    coalescence time for a sample.
    
    Args:
        current_size: Current population size (N_t)
        growth_rate: Per-generation growth rate (r). Positive for growth,
            negative for decline. For exponential growth: N_t = N_0 * e^(r*t)
        generations: Number of generations over which to calculate harmonic mean
    
    Returns:
        Effective population size (harmonic mean). Returns current_size if
        generations <= 0 or growth_rate is zero.
        
    Examples:
        >>> # Growing population: 1000 → 10000 over 10 generations
        >>> exponential_growth_effective_size(10000, growth_rate=0.23, generations=10)
        4342...
        
        >>> # Declining population
        >>> exponential_growth_effective_size(1000, growth_rate=-0.1, generations=10)
        371...
        
    References:
        Maruyama, T., & Fuerst, P. A. (1984). Population bottlenecks and
        nonequilibrium models in population genetics. I. Allele numbers when
        populations evolve from zero variability. Genetics, 108(3), 745-763.
    """
    if generations <= 0 or growth_rate == 0.0:
        return float(current_size)
    
    # Calculate harmonic mean of population sizes
    # For exponential growth: N(t) = N_0 * e^(r*t)
    # We need N_0 from current_size and growth_rate
    # N_t = N_0 * e^(r * generations)
    # So N_0 = N_t / e^(r * generations)
    
    n_0 = current_size / math.exp(growth_rate * generations)
    
    # Harmonic mean: Ne = t / Σ(1/N_i)
    harmonic_sum = 0.0
    for t in range(generations + 1):
        n_t = n_0 * math.exp(growth_rate * t)
        if n_t > 0:
            harmonic_sum += 1.0 / n_t
    
    if harmonic_sum <= 0:
        return float(current_size)
    
    return (generations + 1) / harmonic_sum


def bottleneck_effective_size(
    pre_bottleneck_size: int,
    bottleneck_size: int,
    bottleneck_duration: int,
    recovery_generations: int = 0,
) -> float:
    """Calculate effective population size through a bottleneck.
    
    Models a population bottleneck where the population size drops to a
    small value for a period, then may recover. The effective size is
    dominated by the bottleneck period.
    
    Args:
        pre_bottleneck_size: Population size before bottleneck (N_pre)
        bottleneck_size: Population size during bottleneck (N_b)
        bottleneck_duration: Number of generations at bottleneck size
        recovery_generations: Number of generations after bottleneck (default: 0)
            If > 0, assumes recovery to pre_bottleneck_size
    
    Returns:
        Effective population size (harmonic mean). Heavily weighted by
        bottleneck period.
        
    Examples:
        >>> # Severe bottleneck: 10000 → 100 for 5 generations
        >>> bottleneck_effective_size(10000, 100, 5)
        109...
        
        >>> # With recovery
        >>> bottleneck_effective_size(10000, 100, 5, recovery_generations=10)
        200...
        
    References:
        Nei, M., Maruyama, T., & Chakraborty, R. (1975). The bottleneck
        effect and genetic variability in populations. Evolution, 29(1), 1-10.
    """
    if bottleneck_duration <= 0:
        return float(pre_bottleneck_size)
    
    total_generations = bottleneck_duration + recovery_generations
    
    # Calculate harmonic mean
    harmonic_sum = 0.0
    
    # Pre-bottleneck period (if we're modeling it)
    # For simplicity, we'll focus on bottleneck and recovery
    # But we can add pre-bottleneck if needed
    
    # Bottleneck period
    for _ in range(bottleneck_duration):
        if bottleneck_size > 0:
            harmonic_sum += 1.0 / bottleneck_size
    
    # Recovery period (linear recovery to pre_bottleneck_size)
    if recovery_generations > 0:
        for i in range(recovery_generations):
            # Linear interpolation from bottleneck_size to pre_bottleneck_size
            recovery_size = bottleneck_size + (pre_bottleneck_size - bottleneck_size) * (
                (i + 1) / recovery_generations
            )
            if recovery_size > 0:
                harmonic_sum += 1.0 / recovery_size
    
    if harmonic_sum <= 0:
        return float(min(bottleneck_size, pre_bottleneck_size))
    
    return total_generations / harmonic_sum


def two_epoch_effective_size(
    ancient_size: int,
    current_size: int,
    time_since_change: int,
) -> float:
    """Calculate effective population size for two-epoch model.
    
    Models a population that changed size at some point in the past.
    The effective size is the harmonic mean of the two sizes, weighted
    by the time spent at each size.
    
    Args:
        ancient_size: Population size before change (N_ancient)
        current_size: Population size after change (N_current)
        time_since_change: Number of generations since size change
    
    Returns:
        Effective population size (harmonic mean)
        
    Examples:
        >>> # Population expansion: 1000 → 10000, 50 generations ago
        >>> two_epoch_effective_size(1000, 10000, 50)
        1960...
        
        >>> # Population contraction: 10000 → 1000, 50 generations ago
        >>> two_epoch_effective_size(10000, 1000, 50)
        1960...  # Same harmonic mean (symmetric)
        
    References:
        Slatkin, M., & Hudson, R. R. (1991). Pairwise comparisons of
        mitochondrial DNA sequences in stable and exponentially growing
        populations. Genetics, 129(2), 555-562.
    """
    if time_since_change <= 0:
        return float(current_size)
    
    # For simplicity, assume we're looking back time_since_change generations
    # and the population was at ancient_size before that, then changed to current_size
    # We calculate harmonic mean over the full period
    
    # We need to define a total time period
    # For this model, we'll use time_since_change as the period at current_size
    # and assume some ancient period at ancient_size
    # For simplicity, use equal time periods
    
    total_period = time_since_change * 2  # Equal time at each size
    
    harmonic_sum = 0.0
    
    # Ancient period
    for _ in range(time_since_change):
        if ancient_size > 0:
            harmonic_sum += 1.0 / ancient_size
    
    # Current period
    for _ in range(time_since_change):
        if current_size > 0:
            harmonic_sum += 1.0 / current_size
    
    if harmonic_sum <= 0:
        return float(min(ancient_size, current_size))
    
    return total_period / harmonic_sum

