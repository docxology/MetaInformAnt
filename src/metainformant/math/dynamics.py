"""Population dynamics functions.

This module provides mathematical functions for population dynamics and chaos theory.
"""

from __future__ import annotations

from typing import List

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def logistic_map(r: float, x0: float, n_iterations: int) -> List[float]:
    """Generate logistic map sequence.

    Args:
        r: Growth parameter (0 < r <= 4)
        x0: Initial population (0 < x0 < 1)
        n_iterations: Number of iterations to generate

    Returns:
        List of population values over time
    """
    if not (0 < r <= 4):
        raise ValueError("Growth parameter r must be in (0, 4]")
    if not (0 < x0 < 1):
        raise ValueError("Initial population x0 must be in (0, 1)")

    sequence = [x0]
    x = x0

    for _ in range(n_iterations):
        x = r * x * (1 - x)
        sequence.append(x)

    return sequence


def lotka_volterra_step(prey: float, predator: float, alpha: float = 1.0,
                       beta: float = 0.1, gamma: float = 1.5, delta: float = 0.075,
                       dt: float = 0.01) -> tuple[float, float]:
    """Single step of Lotka-Volterra predator-prey model.

    Args:
        prey: Current prey population
        predator: Current predator population
        alpha: Prey reproduction rate
        beta: Predation rate
        gamma: Predator death rate
        delta: Predator reproduction efficiency
        dt: Time step size

    Returns:
        Tuple of (new_prey, new_predator) after one time step
    """
    if prey < 0 or predator < 0:
        raise ValueError("Population sizes cannot be negative")

    # Lotka-Volterra equations:
    # d(prey)/dt = alpha * prey - beta * prey * predator
    # d(predator)/dt = -gamma * predator + delta * prey * predator

    d_prey = alpha * prey - beta * prey * predator
    d_predator = -gamma * predator + delta * prey * predator

    new_prey = prey + d_prey * dt
    new_predator = predator + d_predator * dt

    # Ensure populations don't go negative
    new_prey = max(0, new_prey)
    new_predator = max(0, new_predator)

    return new_prey, new_predator
