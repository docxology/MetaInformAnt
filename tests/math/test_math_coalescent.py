from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.math.population_genetics.coalescent import (
    coalescent_time_to_mrca,
    expected_coalescent_waiting_times,
    expected_pairwise_diversity,
    expected_pairwise_diversity_from_theta,
    expected_segregating_sites,
    expected_sfs_counts,
    expected_time_to_mrca,
    expected_total_branch_length,
    simulate_coalescent,
    tajima_constants,
    tajimas_D,
    watterson_theta,
)


def test_basic_expectations_monotonic():
    assert expected_time_to_mrca(2, 1000) > 0
    assert expected_time_to_mrca(10, 1000) > expected_time_to_mrca(2, 1000)
    assert expected_total_branch_length(10, 1000) > expected_total_branch_length(2, 1000)


def test_pairwise_diversity():
    assert expected_pairwise_diversity(1000, 1e-8) == 4 * 1000 * 1e-8
    assert expected_pairwise_diversity_from_theta(0.01) == 0.01


def test_watterson_and_expected_S():
    n = 10
    theta = 0.01
    a1 = tajima_constants(n)["a1"]
    # per-site expectations
    assert math.isclose(expected_segregating_sites(n, theta), a1 * theta, rel_tol=1e-9)
    # with sequence length
    L = 1_000
    assert math.isclose(expected_segregating_sites(n, theta, sequence_length=L), a1 * theta * L, rel_tol=1e-9)
    # invert via Watterson
    S = int(round(a1 * theta * L))
    est = watterson_theta(S, n, sequence_length=L)
    assert abs(est - theta) < 1e-2  # coarse due to rounding S


def test_sfs_counts_shape_and_values():
    n = 6
    theta = 0.02
    xs = expected_sfs_counts(n, theta)
    assert len(xs) == n - 1
    # decreasing ~ theta/i
    for i in range(1, n - 1):
        assert xs[i - 1] > xs[i]


def test_waiting_times_and_tajimas_D_stability():
    w = expected_coalescent_waiting_times(5, 1000.0)
    assert len(w) == 4
    # Tajima constants defined and positive where expected
    const = tajima_constants(10)
    assert const["a1"] > 0 and const["a2"] > 0
    # Tajima's D handles zero-variance gracefully
    assert tajimas_D(0, 0.0, 10) == 0.0


def test_expected_time_to_mrca_simple():
    # n=2: sum_{k=2}^2 1/(k(k-1)) = 1/2 => 4N * 1/2 = 2N
    assert abs(expected_time_to_mrca(2, 1000) - 2000.0) < 1e-9


def test_expected_total_branch_length_harmonic():
    # n=4: H_{3} = 1 + 1/2 + 1/3 = 1.833333...
    val = expected_total_branch_length(4, 100)
    assert abs(val - (4 * 100 * (1.0 + 0.5 + 1.0 / 3.0))) < 1e-9


def test_time_to_mrca_validation():
    with pytest.raises(ValueError, match="at least 2"):
        expected_time_to_mrca(1, 1000)
    # Non-positive effective size degrades gracefully to zero
    assert expected_time_to_mrca(5, 0.0) == 0.0


def test_watterson_theta_alias_conventions():
    assert watterson_theta(10, 5) == pytest.approx(watterson_theta(n_sites=5, n_segregating_sites=10))
    assert watterson_theta(10, 5) == pytest.approx(watterson_theta(num_segregating_sites=10, sample_size=5))


def test_watterson_theta_validation():
    with pytest.raises(ValueError, match="> 1"):
        watterson_theta(10, 1)
    with pytest.raises(ValueError, match="negative"):
        watterson_theta(-1, 5)
    with pytest.raises(ValueError, match="must be provided"):
        watterson_theta(None, 5)


def test_waiting_times_validation():
    with pytest.raises(ValueError, match="at least 2"):
        expected_coalescent_waiting_times(1, 1000.0)
    with pytest.raises(ValueError, match="positive"):
        expected_coalescent_waiting_times(5, 0.0)


def test_coalescent_time_to_mrca_alias():
    assert coalescent_time_to_mrca(6, 500.0) == pytest.approx(expected_time_to_mrca(6, 500.0))


def test_simulate_coalescent_basic_shape():
    np.random.seed(1234)
    res = simulate_coalescent(4, effective_size=1000.0, mutation_rate=1e-6)
    assert res["n_samples"] == 4
    assert len(res["coalescence_times"]) == 3
    assert res["coalescence_times"] == sorted(res["coalescence_times"])
    assert res["total_time"] > 0
    assert res["total_branch_length"] >= res["total_time"]  # at least 2 lineages per interval
    assert res["n_segregating_sites"] >= 0
    assert res["theta_watterson"] >= 0.0


def test_simulate_coalescent_two_samples():
    # Regression: n_samples=2 previously crashed via watterson_theta(1000, 1)
    np.random.seed(7)
    res = simulate_coalescent(2, effective_size=1000.0, mutation_rate=0.1)
    assert len(res["coalescence_times"]) == 1
    # Single interval for k=2: E[T] = 2Ne
    assert res["total_time"] > 0


def test_simulate_coalescent_validation():
    with pytest.raises(ValueError, match="at least 2"):
        simulate_coalescent(1)
    with pytest.raises(ValueError, match="positive"):
        simulate_coalescent(4, effective_size=0.0)
    with pytest.raises(ValueError, match="negative"):
        simulate_coalescent(4, mutation_rate=-1e-9)
