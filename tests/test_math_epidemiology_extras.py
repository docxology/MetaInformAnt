from __future__ import annotations

from metainformant.math import herd_immunity_threshold, seir_step


def test_seir_step_and_herd_immunity_threshold():
    S, E, I, R = 0.99, 0.0, 0.01, 0.0
    beta, sigma, gamma = 0.5, 0.2, 0.25
    dt = 0.1
    Sn, En, In, Rn = seir_step(S, E, I, R, beta, sigma, gamma, dt)
    for v in (Sn, En, In, Rn):
        assert v >= 0.0

    R0 = 3.0
    hit = herd_immunity_threshold(R0)
    assert abs(hit - (1.0 - 1.0 / R0)) < 1e-12
