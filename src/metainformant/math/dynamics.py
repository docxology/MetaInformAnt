from __future__ import annotations

import math
from typing import Callable, List, Tuple


def logistic_map(r: float, x0: float, steps: int) -> list[float]:
    """Iterate the logistic map: x_{t+1} = r * x_t * (1 - x_t).
    
    The logistic map is a classic example of discrete-time population dynamics
    that can exhibit chaotic behavior for certain parameter values. Used to model
    population growth with density dependence.
    
    Args:
        r: Growth rate parameter (typically in [0, 4])
        x0: Initial population value (clamped to [0, 1])
        steps: Number of iterations to perform
        
    Returns:
        List of population values over time, starting with x0.
        Values are clamped to [0, 1] at each step.
        
    Examples:
        >>> logistic_map(r=2.5, x0=0.5, steps=5)
        [0.5, 0.625, 0.5859375, 0.606536865234375, 0.5966247408685684, ...]
        
    References:
        May, R. M. (1976). Simple mathematical models with very complicated dynamics.
        Nature, 261(5560), 459-467.
    """
    x = max(0.0, min(1.0, x0))
    result = [x]
    for _ in range(max(0, steps)):
        x = r * x * (1.0 - x)
        x = max(0.0, min(1.0, x))
        result.append(x)
    return result


def lotka_volterra_step(
    prey: float,
    predator: float,
    alpha: float,
    beta: float,
    delta: float,
    gamma: float,
    dt: float = 0.01,
) -> tuple[float, float]:
    """Perform single Euler integration step for Lotka–Volterra predator-prey model.
    
    The Lotka–Volterra equations model population dynamics of interacting species:
    - Prey growth: dx/dt = alpha*x - beta*x*y
    - Predator growth: dy/dt = delta*x*y - gamma*y
    
    Args:
        prey: Current prey population size
        predator: Current predator population size
        alpha: Prey intrinsic growth rate
        beta: Predation rate (predator efficiency at consuming prey)
        delta: Predator growth efficiency (conversion of prey to predators)
        gamma: Predator mortality rate
        dt: Time step size for Euler integration (default 0.01)
        
    Returns:
        Tuple of (next_prey, next_predator) after one time step.
        Populations are clamped to non-negative values.
        
    Examples:
        >>> prey, pred = lotka_volterra_step(
        ...     prey=100.0, predator=10.0,
        ...     alpha=1.0, beta=0.1, delta=0.075, gamma=1.5, dt=0.01
        ... )
        >>> prey > 0 and pred > 0
        True
        
    References:
        Lotka, A. J. (1925). Elements of physical biology. Williams & Wilkins.
        Volterra, V. (1926). Variazioni e fluttuazioni del numero d'individui
        in specie animali conviventi. Mem. Acad. Lincei Roma, 2, 31-113.
    """
    x = max(0.0, prey)
    y = max(0.0, predator)
    x_next = x + dt * (alpha * x - beta * x * y)
    y_next = y + dt * (delta * x * y - gamma * y)
    return max(0.0, x_next), max(0.0, y_next)
