from __future__ import annotations

from typing import Sequence, Tuple


def narrow_sense_heritability(additive_variance: float, phenotypic_variance: float) -> float:
    """h^2 = VA / VP, clamped to [0, 1] when inputs valid."""
    VA = max(0.0, additive_variance)
    VP = max(0.0, phenotypic_variance)
    if VP <= 0:
        return 0.0
    h2 = VA / VP
    if h2 < 0:
        return 0.0
    if h2 > 1:
        return 1.0
    return h2


def breeders_equation_response(selection_differential: float, heritability: float) -> float:
    """R = h^2 * S."""
    h2 = max(0.0, min(1.0, heritability))
    return h2 * selection_differential


def lande_equation_response(gradient: Sequence[float], G_matrix: Sequence[Sequence[float]]) -> Tuple[float, ...]:
    """Multivariate response Δz = G β.

    G is a symmetric positive semidefinite matrix; minimal checks are performed.
    """
    beta = list(gradient)
    G = [list(row) for row in G_matrix]
    if not G or not beta or any(len(row) != len(beta) for row in G):
        return tuple()
    response = [sum(G[i][j] * beta[j] for j in range(len(beta))) for i in range(len(beta))]
    return tuple(response)


def realized_heritability(response: float, selection_differential: float) -> float:
    """h^2_realized = R / S (guarding S=0)."""
    if abs(selection_differential) < 1e-12:
        return 0.0
    h2 = response / selection_differential
    return max(0.0, min(1.0, h2))
