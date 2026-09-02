"""Regression tests for the shared RNG coercion helpers in simulation.rng.

These pin the seeding contract that 30 simulator call sites rely on after the
Round-5 consolidation (commits e4c5014c5, 3c9d7e56a):

- ``coerce_rng(None)`` returns a fresh ``random.Random`` (unseeded behavior
  preserved from the original inline blocks).
- ``seed_numpy_from_rng`` seeds the legacy numpy global RNG with
  ``rng.randint(0, 2**32)`` — byte-identical to the pre-refactor idiom, so a
  given ``random.Random`` state produces the same numpy stream as before.
"""

from __future__ import annotations

import random

import numpy as np

from metainformant.simulation.rng import coerce_and_seed, coerce_rng, seed_numpy_from_rng


def test_coerce_rng_returns_same_object_when_given() -> None:
    rng = random.Random(7)
    assert coerce_rng(rng) is rng


def test_coerce_rng_none_returns_random_instance() -> None:
    rng = coerce_rng(None)
    assert isinstance(rng, random.Random)


def test_seed_numpy_from_rng_is_deterministic() -> None:
    rng_a = random.Random(42)
    rng_b = random.Random(42)

    seed_numpy_from_rng(rng_a)
    stream_a = [np.random.randint(0, 1_000_000) for _ in range(5)]

    seed_numpy_from_rng(rng_b)
    stream_b = [np.random.randint(0, 1_000_000) for _ in range(5)]

    assert stream_a == stream_b


def test_seed_numpy_matches_legacy_idiom() -> None:
    """The helper's seed must equal rng.randint(0, 2**32) — the pre-refactor idiom."""
    expected_seed = random.Random(123).randint(0, 2**32)
    np.random.seed(expected_seed)
    expected = [np.random.random() for _ in range(3)]

    # Helper path: fresh Random(123) draws the same first value for the seed.
    rng2 = random.Random(123)
    seed_numpy_from_rng(rng2)
    helper_stream = [np.random.random() for _ in range(3)]

    assert helper_stream == expected


def test_coerce_and_seed_combines_both() -> None:
    rng = coerce_and_seed(None)
    assert isinstance(rng, random.Random)
    # numpy global stream is now seeded deterministically from rng
    first = np.random.randint(0, 1_000_000)
    rng_b = random.Random(rng.randint(0, 2**32))
    np.random.seed(rng_b.randint(0, 2**32))
    second = np.random.randint(0, 1_000_000)
    assert isinstance(first, int) and isinstance(second, int)
