from __future__ import annotations

from typing import Sequence, Tuple


def kin_selection_response(relatedness: float, benefit: float, cost: float) -> float:
    """Hamilton's rule response r*b - c.

    Positive values imply selection favors the trait.
    """
    return relatedness * benefit - cost


def multilevel_selection_decomposition(
    group_means: Sequence[float],
    individual_deviations: Sequence[float],
    selection_strength_group: float,
    selection_strength_individual: float,
) -> Tuple[float, float, float]:
    """Very small multilevel selection partition.

    Returns (between_group, within_group, total), where
    - between_group = s_g * Var(group_means)
    - within_group = s_i * Var(individual_deviations)
    - total = sum
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


