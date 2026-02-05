"""Quantitative genetics functions.

This module provides functions for quantitative genetics analysis,
including heritability estimation and selection response.
"""

from __future__ import annotations

import statistics
from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def narrow_sense_heritability(additive_genetic_variance: float, phenotypic_variance: float) -> float:
    """Calculate narrow-sense heritability.

    Args:
        additive_genetic_variance: Additive genetic variance (VA)
        phenotypic_variance: Total phenotypic variance (VP)

    Returns:
        Narrow-sense heritability (h² = VA/VP)
    """
    if phenotypic_variance <= 0:
        raise ValueError("Phenotypic variance must be positive")
    if additive_genetic_variance < 0:
        raise ValueError("Additive genetic variance cannot be negative")

    return additive_genetic_variance / phenotypic_variance


def realized_heritability(selection_response: float, selection_differential: float) -> float:
    """Calculate realized heritability from selection experiment.

    Args:
        selection_response: Response to selection (R)
        selection_differential: Selection differential (S)

    Returns:
        Realized heritability (h²)
    """
    if selection_differential == 0:
        return 0.0

    return selection_response / selection_differential


def lande_equation_response(selection_gradient: List[float], genetic_variance_covariance: List[List[float]]) -> float:
    """Calculate evolutionary response using Lande's equation.

    Lande's equation predicts the evolutionary response to selection as:
    R = G * β, where G is the genetic variance-covariance matrix and β is the selection gradient.

    Args:
        selection_gradient: Vector of selection gradients
        genetic_variance_covariance: Genetic variance-covariance matrix (G matrix)

    Returns:
        Predicted evolutionary response

    Examples:
        >>> gradient = [0.1, -0.05]
        >>> G = [[1.0, 0.3], [0.3, 0.8]]
        >>> response = lande_equation_response(gradient, G)
    """
    import numpy as np

    if not selection_gradient or not genetic_variance_covariance:
        raise ValueError("Selection gradient and G matrix cannot be empty")

    n_traits = len(selection_gradient)
    if len(genetic_variance_covariance) != n_traits:
        raise ValueError("G matrix dimensions must match selection gradient length")

    for row in genetic_variance_covariance:
        if len(row) != n_traits:
            raise ValueError("G matrix must be square")

    # Convert to numpy arrays
    beta = np.array(selection_gradient)
    G = np.array(genetic_variance_covariance)

    # Calculate response: R = G * β
    response_vector = np.dot(G, beta)

    # Return the magnitude of the response vector
    return float(np.linalg.norm(response_vector))
