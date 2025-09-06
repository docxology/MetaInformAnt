from __future__ import annotations

from typing import Iterable, List


def replicator_step(frequencies: Iterable[float], payoff_matrix: Iterable[Iterable[float]]) -> list[float]:
    """Discrete replicator update x_i' = x_i * ( (A x)_i ) / (x^T A x).

    Frequencies are clamped to be non-negative and renormalized to sum to 1.
    """
    x = [max(0.0, float(v)) for v in frequencies]
    total = sum(x)
    if total <= 0:
        return [0.0 for _ in x]
    x = [v / total for v in x]
    A = [list(map(float, row)) for row in payoff_matrix]
    if not A or any(len(row) != len(x) for row in A):
        return x
    Ax = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(x))]
    xTAx = sum(x[i] * Ax[i] for i in range(len(x)))
    if xTAx <= 0:
        return x
    x_next = [x[i] * (Ax[i] / xTAx) for i in range(len(x))]
    # renormalize to guard against numeric drift
    s = sum(max(0.0, v) for v in x_next)
    return [max(0.0, v) / s if s > 0 else 0.0 for v in x_next]


def replicator_derivative(frequencies: Iterable[float], payoff_matrix: Iterable[Iterable[float]]) -> list[float]:
    """Continuous-time replicator equation derivative: dx_i/dt = x_i[(A x)_i - x^T A x]."""
    x = [max(0.0, float(v)) for v in frequencies]
    total = sum(x)
    x = [v / total for v in x] if total > 0 else [0.0 for _ in x]
    A = [list(map(float, row)) for row in payoff_matrix]
    if not A or any(len(row) != len(x) for row in A):
        return [0.0 for _ in x]
    Ax = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(x))]
    xTAx = sum(x[i] * Ax[i] for i in range(len(x)))
    return [x[i] * (Ax[i] - xTAx) for i in range(len(x))]
