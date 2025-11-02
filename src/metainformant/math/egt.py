from __future__ import annotations

from typing import Iterable, List


def replicator_step(frequencies: Iterable[float], payoff_matrix: Iterable[Iterable[float]]) -> list[float]:
    """Perform one discrete-time step of the replicator equation.
    
    The replicator equation models frequency-dependent selection in evolutionary
    game theory. Strategies with above-average payoffs increase in frequency.
    
    Args:
        frequencies: Current strategy frequencies (must sum to 1)
        payoff_matrix: Payoff matrix A where A[i][j] is payoff to strategy i
            when playing against strategy j
            
    Returns:
        Updated frequencies after one generation, renormalized to sum to 1.
        Formula: x_i' = x_i × (Ax)_i / (x^T A x)
        
    Examples:
        >>> freqs = [0.5, 0.3, 0.2]
        >>> payoff = [[1.0, 0.5, 0.3], [0.5, 1.0, 0.7], [0.3, 0.7, 1.0]]
        >>> new_freqs = replicator_step(freqs, payoff)
        >>> abs(sum(new_freqs) - 1.0) < 1e-10
        True
        
    References:
        Hofbauer, J., & Sigmund, K. (1998). Evolutionary games and population
        dynamics. Cambridge University Press.
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
    """Calculate continuous-time replicator equation derivative.
    
    Returns the rate of change for each strategy frequency under continuous-time
    replicator dynamics. Used for continuous-time simulations and stability analysis.
    
    Args:
        frequencies: Current strategy frequencies
        payoff_matrix: Payoff matrix A where A[i][j] is payoff to strategy i
            when playing against strategy j
            
    Returns:
        List of derivatives dx_i/dt for each strategy. Formula:
        dx_i/dt = x_i × [(Ax)_i - x^T A x]
        
    Examples:
        >>> freqs = [0.5, 0.3, 0.2]
        >>> payoff = [[1.0, 0.5], [0.5, 1.0]]
        >>> derivatives = replicator_derivative(freqs[:2], payoff)
        >>> len(derivatives)
        2
        
    References:
        Hofbauer, J., & Sigmund, K. (1998). Evolutionary games and population
        dynamics. Cambridge University Press.
    """
    x = [max(0.0, float(v)) for v in frequencies]
    total = sum(x)
    x = [v / total for v in x] if total > 0 else [0.0 for _ in x]
    A = [list(map(float, row)) for row in payoff_matrix]
    if not A or any(len(row) != len(x) for row in A):
        return [0.0 for _ in x]
    Ax = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(x))]
    xTAx = sum(x[i] * Ax[i] for i in range(len(x)))
    return [x[i] * (Ax[i] - xTAx) for i in range(len(x))]
