"""Shared RNG coercion for simulation models.

Every simulator in ``metainformant.simulation.models`` accepts an optional
``random.Random`` and needs a deterministic numpy seed derived from it. This
module centralizes that boilerplate so the seeding contract has one
implementation (duplication across 18 call sites was the drift risk).

Behavior is API-preserving: ``coerce_rng(None)`` returns a fresh unseeded
``random.Random`` exactly as the previous inline ``if rng is None`` blocks did,
and the numpy seed derivation is byte-identical to the previous
``np.random.seed(rng.randint(0, 2**32))`` idiom.
"""

from __future__ import annotations

import random

import numpy as np


def coerce_rng(rng: random.Random | None) -> random.Random:
    """Return *rng*, or a fresh unseeded ``random.Random`` if None."""
    if rng is None:
        return random.Random()
    return rng


def seed_numpy_from_rng(rng: random.Random) -> None:
    """Seed the legacy global numpy RNG deterministically from *rng*.

    Uses ``rng.randint(0, 2**32)`` — identical to the previous inline idiom at
    every call site — so simulated streams are unchanged for a given rng state.
    """
    np.random.seed(rng.randint(0, 2**32))


def coerce_and_seed(rng: random.Random | None) -> random.Random:
    """Coerce *rng* (None -> fresh Random) and seed numpy from it.

    Returns the coerced rng. Equivalent to the previous two-line inline idiom
    ``if rng is None: rng = random.Random()`` followed by
    ``np.random.seed(rng.randint(0, 2**32))``.
    """
    rng = coerce_rng(rng)
    seed_numpy_from_rng(rng)
    return rng
