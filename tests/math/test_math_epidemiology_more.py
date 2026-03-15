from __future__ import annotations

from metainformant.math.epidemiology.models import effective_reproduction_number, sis_step


def test_sis_step_and_effective_R():
    S, I = 0.9, 0.1
    beta, gamma = 0.6, 0.2
    dt = 0.05
    Sn, In = sis_step(S, I, beta, gamma, dt)
    assert Sn >= 0.0 and In >= 0.0
    R0 = beta / gamma
    Re = effective_reproduction_number(R0, Sn)
    assert abs(Re - R0 * Sn) < 1e-12
