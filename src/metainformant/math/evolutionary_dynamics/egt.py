"""Evolutionary game theory functions.

This module provides functions for evolutionary game theory simulations,
including replicator dynamics and fitness calculations.
"""

from __future__ import annotations

from typing import List

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def replicator_derivative(fitnesses: List[float], frequencies: List[float]) -> List[float]:
    """Calculate the derivative of replicator dynamics.

    dx_i/dt = x_i * (f_i - mean fitness), with mean fitness the frequency-weighted mean of the
    ``fitnesses``. The derivative is defined whatever the mean fitness is, including zero (or negative)
    fitnesses.

    Args:
        fitnesses: Fitness values for each strategy
        frequencies: Current frequencies of each strategy

    Returns:
        List of derivatives for each strategy frequency

    Examples:
        >>> fitnesses = [1.1, 0.9, 1.0]
        >>> frequencies = [0.4, 0.3, 0.3]
        >>> derivs = replicator_derivative(fitnesses, frequencies)
    """
    if len(fitnesses) != len(frequencies):
        raise ValueError("Fitnesses and frequencies must have same length")

    if not frequencies:
        return []

    # Mean fitness
    mean_fitness = sum(f * freq for f, freq in zip(fitnesses, frequencies))

    # Replicator equation: dx_i/dt = x_i * (f_i - mean_f)
    derivatives = []
    for fitness, freq in zip(fitnesses, frequencies):
        derivative = freq * (fitness - mean_fitness)
        derivatives.append(derivative)

    return derivatives


def replicator_step(fitnesses: List[float], frequencies: List[float], time_step: float = 0.01) -> List[float]:
    """Perform one step of replicator dynamics.

    Args:
        fitnesses: Fitness values for each strategy
        frequencies: Current frequencies of each strategy
        time_step: Time step for Euler integration

    Returns:
        Updated frequencies after one time step: Euler step, clipped to [0, 1], renormalised to sum to 1

    Raises:
        ValueError: If a frequency is negative, the frequencies sum to zero, or the step drives every
            frequency to zero (a time step too large for frequencies that do not sum to 1)

    Examples:
        >>> fitnesses = [1.1, 0.9, 1.0]
        >>> frequencies = [0.4, 0.3, 0.3]
        >>> new_freqs = replicator_step(fitnesses, frequencies, 0.01)
    """
    if any(freq < 0 for freq in frequencies):
        raise ValueError("Frequencies must be non-negative")
    if frequencies and sum(frequencies) <= 0:
        raise ValueError("Frequencies must sum to a positive value")

    derivatives = replicator_derivative(fitnesses, frequencies)

    # Euler integration
    new_frequencies = []
    for freq, deriv in zip(frequencies, derivatives):
        new_freq = freq + time_step * deriv
        # Ensure frequencies stay in [0, 1]
        new_freq = max(0.0, min(1.0, new_freq))
        new_frequencies.append(new_freq)

    # Renormalize to ensure they sum to 1
    total = sum(new_frequencies)
    if new_frequencies and total <= 0:
        raise ValueError("The step drove every frequency to zero: reduce the time step or normalise the frequencies")
    return [f / total for f in new_frequencies]
