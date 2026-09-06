from __future__ import annotations

import pytest

from metainformant.math.epidemiology.models import herd_immunity_threshold, seir_step, sir_step


def test_seir_step_and_herd_immunity_threshold():
    susceptible, exposed, infected, recovered = 0.99, 0.0, 0.01, 0.0
    beta, sigma, gamma = 0.5, 0.2, 0.25
    dt = 0.1
    Sn, En, In, Rn = seir_step(susceptible, exposed, infected, recovered, beta, sigma, gamma, dt)
    for v in (Sn, En, In, Rn):
        assert v >= 0.0

    R0 = 3.0
    hit = herd_immunity_threshold(R0)
    assert abs(hit - (1.0 - 1.0 / R0)) < 1e-12


def test_sir_step_missing_rates_raise():
    with pytest.raises(ValueError, match="required"):
        sir_step(0.99, 0.01, 0.0)


def test_sir_step_beta_gamma_keywords_accepted():
    Sn, In, Rn = sir_step(0.99, 0.01, 0.0, beta=0.5, gamma=0.25, dt=0.1)
    assert Sn >= 0.0 and In >= 0.0 and Rn >= 0.0
    # Population is conserved by the step (up to clipping)
    assert abs((Sn + In + Rn) - 1.0) < 1e-9
