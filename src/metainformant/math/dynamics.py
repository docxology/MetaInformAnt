from __future__ import annotations

import math
from typing import Callable, List, Tuple


def logistic_map(r: float, x0: float, steps: int) -> list[float]:
    """Iterate the logistic map x_{t+1} = r x_t (1 - x_t)."""
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
    """Single Euler step for Lotka–Volterra predator-prey dynamics."""
    x = max(0.0, prey)
    y = max(0.0, predator)
    x_next = x + dt * (alpha * x - beta * x * y)
    y_next = y + dt * (delta * x * y - gamma * y)
    return max(0.0, x_next), max(0.0, y_next)



