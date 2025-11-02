from __future__ import annotations

from typing import Sequence, Tuple


def kin_selection_response(relatedness: float, benefit: float, cost: float) -> float:
    """Calculate Hamilton's rule for kin selection.
    
    Hamilton's rule states that a trait will be favored by natural selection when
    the genetic relatedness (r) times the benefit to the recipient (b) exceeds
    the cost to the actor (c): r * b > c
    
    Args:
        relatedness: Genetic relatedness between actor and recipient (r)
        benefit: Benefit to recipient (b)
        cost: Cost to actor (c)
        
    Returns:
        Hamilton's rule value: r * b - c. Positive values indicate selection
        favors the trait, negative values indicate selection against.
        
    Examples:
        >>> kin_selection_response(relatedness=0.5, benefit=0.4, cost=0.1)
        0.1
        >>> kin_selection_response(relatedness=0.25, benefit=0.3, cost=0.2)
        -0.125
        
    References:
        Hamilton, W. D. (1964). The genetical evolution of social behaviour.
        Journal of Theoretical Biology, 7(1), 1-52.
    """
    return relatedness * benefit - cost


def multilevel_selection_decomposition(
    group_means: Sequence[float],
    individual_deviations: Sequence[float],
    selection_strength_group: float,
    selection_strength_individual: float,
) -> Tuple[float, float, float]:
    """Decompose multilevel selection into between-group and within-group components.
    
    Separates selection acting at different levels: selection between groups
    (group-level fitness differences) and selection within groups (individual-level
    fitness differences). Useful for understanding the balance of group vs.
    individual selection.
    
    Args:
        group_means: Sequence of mean trait values for each group
        individual_deviations: Sequence of individual deviations from their
            group means (individual_value - group_mean for each individual)
        selection_strength_group: Selection coefficient for group-level selection (s_g)
        selection_strength_individual: Selection coefficient for individual-level
            selection (s_i)
            
    Returns:
        Tuple of (between_group_selection, within_group_selection, total_selection):
        - between_group: Selection due to group-level variance = s_g × Var(group_means)
        - within_group: Selection due to individual variance = s_i × Var(deviations)
        - total: Sum of both components
        
    Examples:
        >>> group_means = [0.5, 0.7, 0.3]
        >>> deviations = [0.1, -0.1, 0.05, -0.05, 0.2, -0.2]
        >>> between, within, total = multilevel_selection_decomposition(
        ...     group_means, deviations, selection_strength_group=0.5, selection_strength_individual=0.3
        ... )
        >>> total == between + within
        True
        
    References:
        Price, G. R. (1972). Extension of covariance selection mathematics.
        Annals of Human Genetics, 35(4), 485-490.
    """
    if not group_means or not individual_deviations:
        return 0.0, 0.0, 0.0

    def variance(xs: Sequence[float]) -> float:
        if not xs:
            return 0.0
        n = len(xs)
        mu = sum(xs) / n
        return sum((x - mu) ** 2 for x in xs) / n

    between = selection_strength_group * variance(group_means)
    within = selection_strength_individual * variance(individual_deviations)
    return between, within, between + within
