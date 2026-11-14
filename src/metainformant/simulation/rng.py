"""Shared RNG coercion for simulation models.

Every simulator in ``metainformant.simulation.models`` accepts an optional
``random.Random`` and needs a deterministic numpy seed derived from it. This
module holds that boilerplate so the seeding contract has one implementation:
``coerce_rng(None)`` returns a fresh unseeded ``random.Random``, and numpy is
seeded with ``np.random.seed(rng.randint(0, 2**32))`` (2**32 mapped to 0).
"""

from __future__ import annotations

import random

import numpy as np


def coerce_rng(rng: random.Random | None) -> random.Random:
    """Return *rng*, or a fresh unseeded ``random.Random`` if None."""
    if rng is None:
        return random.Random()
    return rng


def _numpy_seed(value: int) -> int:
    """Map a draw from ``randint(0, 2**32)`` onto numpy's seed range [0, 2**32).

    ``randint`` includes its upper bound, which ``np.random.seed`` rejects. Every
    other value passes through unchanged.
    """
    return value % 2**32


def seed_numpy_from_rng(rng: random.Random) -> None:
    """Seed the legacy global numpy RNG deterministically from *rng*.

    Draws ``rng.randint(0, 2**32)``, so a given rng state always gives the same
    simulated stream. The one draw numpy cannot accept, ``2**32``, maps to 0.
    """
    np.random.seed(_numpy_seed(rng.randint(0, 2**32)))


def coerce_and_seed(rng: random.Random | None) -> random.Random:
    """Coerce *rng* (None -> fresh Random) and seed numpy from it.

    Returns the coerced rng: ``if rng is None: rng = random.Random()`` followed by
    ``np.random.seed(rng.randint(0, 2**32))``.
    """
    rng = coerce_rng(rng)
    seed_numpy_from_rng(rng)
    return rng
