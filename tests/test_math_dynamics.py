from __future__ import annotations

from metainformant.math import logistic_map, lotka_volterra_step


def test_logistic_map_basic():
    seq = logistic_map(2.5, 0.2, 5)
    assert len(seq) == 6
    assert all(0.0 <= x <= 1.0 for x in seq)


def test_lotka_volterra_step_nonnegative():
    nx, ny = lotka_volterra_step(prey=10.0, predator=5.0, alpha=1.0, beta=0.1, delta=0.1, gamma=1.5, dt=0.05)
    assert nx >= 0.0 and ny >= 0.0
