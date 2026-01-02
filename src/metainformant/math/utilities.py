"""Mathematical utility functions.

This module provides general mathematical utility functions for statistical analysis.
"""

from __future__ import annotations

from typing import List

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def correlation(x: List[float], y: List[float]) -> float:
    """Calculate Pearson correlation coefficient.

    Args:
        x: First variable values
        y: Second variable values

    Returns:
        Pearson correlation coefficient
    """
    if len(x) != len(y):
        raise ValueError("Input lists must have equal length")

    if len(x) < 2:
        return 0.0

    # Convert to numpy arrays
    x_arr = np.array(x)
    y_arr = np.array(y)

    # Calculate correlation
    return np.corrcoef(x_arr, y_arr)[0, 1]


def correlation_coefficient(x: List[float], y: List[float]) -> float:
    """Calculate Pearson correlation coefficient (alias for correlation).

    Args:
        x: First variable values
        y: Second variable values

    Returns:
        Pearson correlation coefficient
    """
    return correlation(x, y)


def linear_regression(x: List[float], y: List[float]) -> tuple[float, float, float]:
    """Perform linear regression and return slope, intercept, and r-squared.

    Args:
        x: Independent variable values
        y: Dependent variable values

    Returns:
        Tuple of (slope, intercept, r_squared)
    """
    if len(x) != len(y):
        raise ValueError("Input lists must have equal length")

    if len(x) < 2:
        raise ValueError("Need at least 2 data points for regression")

    # Convert to numpy arrays
    x_arr = np.array(x)
    y_arr = np.array(y)

    # Perform linear regression
    slope, intercept = np.polyfit(x_arr, y_arr, 1)

    # Calculate R-squared
    y_pred = slope * x_arr + intercept
    ss_res = np.sum((y_arr - y_pred) ** 2)
    ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)

    if ss_tot == 0:
        r_squared = 1.0 if ss_res == 0 else 0.0
    else:
        r_squared = 1 - (ss_res / ss_tot)

    return slope, intercept, r_squared


def r_squared(x: List[float], y: List[float]) -> float:
    """Calculate R-squared for linear regression of x vs y.

    Args:
        x: Independent variable values
        y: Dependent variable values

    Returns:
        R-squared value (0 to 1)
    """
    _, _, r_squared = linear_regression(x, y)
    return r_squared
