"""Shared numeric helpers for ML modules that treat NumPy as optional.

Several ML subpackages (automl, interpretability) accept either NumPy arrays
or plain nested lists and must keep working when NumPy is not installed.
These helpers centralize that conversion logic.
"""

from __future__ import annotations

import numbers
import random
from typing import Any, cast

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:  # pragma: no cover - exercised only without numpy
    HAS_NUMPY = False
    np = cast("Any", None)


def resolve_python_rng(random_state: Any = None) -> random.Random:
    """Return the ``random.Random`` generator that ``random_state`` stands for.

    The pure-Python ML routines draw all their randomness from the generator this returns, so they never seed (or
    consume) the global ``random`` state. An int gives a fresh, reproducible generator, ``None`` a generator seeded from
    system entropy, and a ``random.Random`` (or a NumPy ``Generator``, from which a seed is drawn) is used as given.

    Args:
        random_state: ``None``, an integer seed, a ``random.Random`` or a ``numpy.random.Generator``.

    Returns:
        A ``random.Random`` instance.
    """
    if isinstance(random_state, random.Random):
        return random_state
    if HAS_NUMPY and isinstance(random_state, np.random.Generator):
        return random.Random(int(random_state.integers(0, 2**63 - 1)))
    if isinstance(random_state, numbers.Integral):
        return random.Random(int(random_state))
    return random.Random(random_state)


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
