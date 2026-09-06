"""Shared numeric helpers for ML modules that treat NumPy as optional.

Several ML subpackages (automl, interpretability) accept either NumPy arrays
or plain nested lists and must keep working when NumPy is not installed.
These helpers centralize that conversion logic.
"""

from __future__ import annotations

from typing import Any

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:  # pragma: no cover - exercised only without numpy
    HAS_NUMPY = False
    np = None


def to_2d_list(X: Any) -> list[list[float]]:
    """Convert an input matrix to a list of lists of floats.

    Args:
        X: Input matrix (NumPy 2D array or sequence of rows).

    Returns:
        Matrix as list of lists of floats.
    """
    if HAS_NUMPY and isinstance(X, np.ndarray):
        return [[float(X[i, j]) for j in range(X.shape[1])] for i in range(X.shape[0])]
    return [[float(v) for v in row] for row in X]


def to_1d_list(y: Any) -> list[float]:
    """Convert an input vector to a list of floats.

    Args:
        y: Input vector (NumPy array or sequence).

    Returns:
        Vector as list of floats.
    """
    if HAS_NUMPY and isinstance(y, np.ndarray):
        return [float(v) for v in y.ravel()]
    return [float(v) for v in y]


def get_shape(X: Any) -> tuple[int, int]:
    """Get the shape of a 2D matrix.

    Args:
        X: Input matrix.

    Returns:
        Tuple of (n_rows, n_cols).
    """
    if HAS_NUMPY and isinstance(X, np.ndarray):
        return int(X.shape[0]), int(X.shape[1])
    n_rows = len(X)
    n_cols = len(X[0]) if n_rows > 0 else 0
    return n_rows, n_cols
