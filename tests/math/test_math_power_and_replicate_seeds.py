"""Tests for round-4 additions: power/sample-size helpers and deterministic replicate seeds.

All tests use real computation (scipy normal quantiles, SHA-256 hashing,
actual sequence evolution) - real-implementation policy: no test doubles.
"""

from __future__ import annotations

import hashlib

from metainformant.math.population_genetics.statistics import (
    deterministic_replicate_seeds,
    normal_approximate_power,
    sample_size_for_power,
)
from metainformant.simulation.models.sequences import deterministic_evolution_replicates


def test_power_increases_with_sample_size():
    """Power must be monotonically non-decreasing in per-group n."""
    powers = [normal_approximate_power(0.5, n) for n in (4, 8, 16, 32, 64)]
    assert all(b >= a for a, b in zip(powers, powers[1:]))
    assert 0.0 <= powers[0] < powers[-1] < 1.0


def test_power_known_values():
    """Check against the textbook normal-approximation formula at fixed points."""
    # d=0.5, n=64/group, two-sided alpha=0.05: ncp = 0.5*sqrt(32) = 2.828
    from scipy.stats import norm

    expected = float(norm.cdf(0.5 * (64 / 2) ** 0.5 - norm.ppf(0.975)))
    got = normal_approximate_power(0.5, 64)
    assert abs(got - expected) < 1e-12
    # One-sided is strictly more powerful than two-sided at same n
    assert normal_approximate_power(0.5, 64, two_sided=False) > got
    # Zero effect: power at/below alpha level (no signal)
    assert normal_approximate_power(0.0, 100) <= 0.05 + 1e-9


def test_power_validates_arguments():
    for kwargs in (
        {"effect_size": -1.0, "sample_size_per_group": 10},
        {"effect_size": 0.5, "sample_size_per_group": 1},
        {"effect_size": 0.5, "sample_size_per_group": 10, "alpha": 0.0},
        {"effect_size": 0.5, "sample_size_per_group": 10, "alpha": 1.0},
    ):
        try:
            normal_approximate_power(**kwargs)
        except ValueError:
            pass
        else:  # pragma: no cover - assertion path
            raise AssertionError(f"expected ValueError for {kwargs}")


def test_sample_size_roundtrip_with_power():
    """Solved n must achieve at least the target power, and n-1 must not."""
    for effect in (0.2, 0.5, 0.8):
        n = sample_size_for_power(effect, target_power=0.8)
        assert n >= 2
        assert normal_approximate_power(effect, n) >= 0.8
        if n > 2:
            assert normal_approximate_power(effect, n - 1) < 0.8
    # Exact-formula result: d=0.8, power 0.8, two-sided alpha 0.05 ->
    # 2*((1.959964+0.841621)/0.8)^2 = 24.53 -> 25/group (the folk "26" comes
    # from rounding z_alpha to 2.0)
    assert sample_size_for_power(0.8) == 25


def test_sample_size_validates_arguments():
    for kwargs in (
        {"effect_size": 0.0},
        {"effect_size": 0.5, "target_power": 0.05},  # power == alpha band edge
        {"effect_size": 0.5, "target_power": 1.5},
        {"effect_size": 0.5, "alpha": 1.5},
    ):
        try:
            sample_size_for_power(**kwargs)
        except ValueError:
            pass
        else:  # pragma: no cover
            raise AssertionError(f"expected ValueError for {kwargs}")


def test_replicate_seeds_deterministic_and_distinct():
    a = deterministic_replicate_seeds(42, 16)
    b = deterministic_replicate_seeds(42, 16)
    assert a == b
    assert len(set(a)) == 16
    assert all(0 <= s < 2**63 for s in a)
    # Independent derivation from the documented formula
    expected0 = int.from_bytes(hashlib.sha256(b"42:0").digest()[:8], "big") & 0x7FFFFFFFFFFFFFFF
    assert a[0] == expected0
    # Different base seed -> different seeds
    assert deterministic_replicate_seeds(43, 16) != a


def test_replicate_seeds_order_independent_prefix():
    """Deriving 16 seeds and taking the first 4 equals deriving 4 directly."""
    assert deterministic_replicate_seeds(7, 16)[:4] == deterministic_replicate_seeds(7, 4)


def test_replicate_seeds_validates_arguments():
    try:
        deterministic_replicate_seeds(1, 0)
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for n_replicates=0")


def test_deterministic_evolution_replicates_reproducible():
    seq = "ACGT" * 50
    r1 = deterministic_evolution_replicates(seq, 4, base_seed=123, generations=20, mutation_rate=0.01)
    r2 = deterministic_evolution_replicates(seq, 4, base_seed=123, generations=20, mutation_rate=0.01)
    assert r1 == r2
    assert len(r1) == 4
    assert all(len(s) == len(seq) for s in r1)


def test_deterministic_evolution_replicates_prefix_stability():
    """A 6-replicate sweep shares its first 3 replicates with a 3-replicate sweep."""
    seq = "ACGTACGTAA" * 20
    short = deterministic_evolution_replicates(seq, 3, base_seed=5, generations=10)
    long = deterministic_evolution_replicates(seq, 6, base_seed=5, generations=10)
    assert long[:3] == short


def test_deterministic_evolution_replicates_replicates_differ():
    """Independent seeds must give independent evolutionary outcomes."""
    seq = "ACGT" * 100
    reps = deterministic_evolution_replicates(seq, 8, base_seed=9, generations=50, mutation_rate=0.05)
    assert len(set(reps)) > 1


def test_deterministic_evolution_replicates_validates_arguments():
    try:
        deterministic_evolution_replicates("ACGT", 0, base_seed=1)
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for n_replicates=0")
