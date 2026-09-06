from __future__ import annotations

import pytest

from metainformant.math.evolutionary_dynamics.core import (
    logistic_map,
    lotka_volterra_step,
    replicator_derivative,
    replicator_step,
)
from metainformant.math.evolutionary_dynamics.egt import (
    replicator_derivative as egt_replicator_derivative,
    replicator_step as egt_replicator_step,
)


def test_logistic_map_basic():
    seq = logistic_map(2.5, 0.2, 5)
    assert len(seq) == 6
    assert all(0.0 <= x <= 1.0 for x in seq)


def test_lotka_volterra_step_nonnegative():
    nx, ny = lotka_volterra_step(prey=10.0, predator=5.0, alpha=1.0, beta=0.1, delta=0.1, gamma=1.5, dt=0.05)
    assert nx >= 0.0 and ny >= 0.0


class TestLogisticMapErrors:
    def test_missing_iterations_raises(self):
        with pytest.raises(ValueError, match="n_iterations or steps"):
            logistic_map(2.5, 0.2)

    def test_r_out_of_range_raises(self):
        with pytest.raises(ValueError, match="r must be"):
            logistic_map(4.5, 0.2, 5)

    def test_x0_out_of_range_raises(self):
        with pytest.raises(ValueError, match="x0 must be"):
            logistic_map(2.5, 1.5, 5)


class TestLotkaVolterraErrors:
    def test_negative_populations_raise(self):
        with pytest.raises(ValueError, match="negative"):
            lotka_volterra_step(prey=-1.0, predator=5.0)


class TestCoreReplicator:
    def test_derivative_zero_at_uniform_with_equal_fitness(self):
        # Equal fitnesses: no strategy has an advantage, all derivatives are 0
        derivs = replicator_derivative([0.25, 0.25, 0.25, 0.25], [[1.0, 1.0, 1.0, 1.0]] * 4)
        assert all(d == pytest.approx(0.0) for d in derivs)

    def test_step_preserves_unit_sum(self):
        A = [[1.0, 0.0], [0.0, 1.5]]
        new_freqs = replicator_step([0.5, 0.5], A, dt=0.05)
        assert sum(new_freqs) == pytest.approx(1.0)
        assert all(f >= 0.0 for f in new_freqs)


class TestEgtReplicator:
    def test_derivative_direction(self):
        # Strategy 2 has higher fitness, so its frequency should increase
        derivs = egt_replicator_derivative([1.0, 2.0], [0.5, 0.5])
        assert derivs[0] < 0.0
        assert derivs[1] > 0.0

    def test_derivative_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="same length"):
            egt_replicator_derivative([1.0, 2.0], [0.5])

    def test_derivative_empty(self):
        assert egt_replicator_derivative([], []) == []

    def test_zero_mean_fitness_returns_zeros(self):
        assert egt_replicator_derivative([1.0, -1.0], [0.5, 0.5]) == [0.0, 0.0]

    def test_step_renormalizes(self):
        new_freqs = egt_replicator_step([1.0, 2.0], [0.5, 0.5], time_step=0.1)
        assert sum(new_freqs) == pytest.approx(1.0)
        assert all(0.0 <= f <= 1.0 for f in new_freqs)

    def test_step_all_zero_frequencies_resets_uniform(self):
        # Zero frequencies leave nothing to renormalize; the step resets to uniform
        new_freqs = egt_replicator_step([1.0, 2.0], [0.0, 0.0], time_step=0.1)
        assert new_freqs == [0.5, 0.5]
