"""Tests for epidemiology models: SIR, SEIR, SIS, R0, herd immunity.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import pytest

from metainformant.math.epidemiology.models import (
    basic_reproduction_number,
    effective_reproduction_number,
    herd_immunity_threshold,
    seir_step,
    sir_step,
    sis_step,
)


# ---------------------------------------------------------------------------
# basic_reproduction_number
# ---------------------------------------------------------------------------


class TestBasicReproductionNumber:
    def test_r0_calculation(self):
        r0 = basic_reproduction_number(transmission_rate=0.3, recovery_rate=0.1)
        assert r0 == pytest.approx(3.0)

    def test_r0_less_than_one(self):
        r0 = basic_reproduction_number(transmission_rate=0.05, recovery_rate=0.1)
        assert r0 < 1.0

    def test_contact_rate_ignored(self):
        r0 = basic_reproduction_number(transmission_rate=0.5, recovery_rate=0.25, contact_rate=10.0)
        assert r0 == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# effective_reproduction_number
# ---------------------------------------------------------------------------


class TestEffectiveReproductionNumber:
    def test_full_susceptibility(self):
        re = effective_reproduction_number(R0=3.0, susceptible_fraction=1.0)
        assert re == pytest.approx(3.0)

    def test_partial_susceptibility(self):
        re = effective_reproduction_number(R0=3.0, susceptible_fraction=0.5)
        assert re == pytest.approx(1.5)

    def test_zero_susceptibility(self):
        re = effective_reproduction_number(R0=3.0, susceptible_fraction=0.0)
        assert re == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# herd_immunity_threshold
# ---------------------------------------------------------------------------


class TestHerdImmunityThreshold:
    def test_measles_like(self):
        # R0 ~ 15 => threshold ~ 0.933
        threshold = herd_immunity_threshold(R0=15.0)
        assert threshold == pytest.approx(1.0 - 1.0 / 15.0, abs=1e-6)

    def test_r0_equals_2(self):
        threshold = herd_immunity_threshold(R0=2.0)
        assert threshold == pytest.approx(0.5)

    def test_low_r0(self):
        threshold = herd_immunity_threshold(R0=1.5)
        assert 0.0 < threshold < 1.0


# ---------------------------------------------------------------------------
# sir_step
# ---------------------------------------------------------------------------


class TestSIRStep:
    def test_conservation_of_population(self):
        S, I, R = sir_step(990.0, 10.0, 0.0, 0.3, 0.1, dt=1.0)
        assert S + I + R == pytest.approx(1000.0, abs=0.01)

    def test_infection_decreases_susceptible(self):
        S, I, R = sir_step(990.0, 10.0, 0.0, 0.5, 0.1, dt=1.0)
        assert S < 990.0

    def test_recovery_increases_recovered(self):
        S, I, R = sir_step(990.0, 10.0, 0.0, 0.3, 0.1, dt=1.0)
        assert R > 0.0

    def test_no_infection_when_no_infected(self):
        S, I, R = sir_step(1000.0, 0.0, 0.0, 0.3, 0.1, dt=1.0)
        assert S == pytest.approx(1000.0)
        assert I == pytest.approx(0.0)

    def test_small_timestep(self):
        S, I, R = sir_step(990.0, 10.0, 0.0, 0.3, 0.1, dt=0.01)
        # Should change less than dt=1.0
        assert abs(S - 990.0) < 1.0


# ---------------------------------------------------------------------------
# seir_step
# ---------------------------------------------------------------------------


class TestSEIRStep:
    def test_conservation_of_population(self):
        S, E, I, R = seir_step(980.0, 10.0, 10.0, 0.0, 0.3, 0.2, 0.1, dt=1.0)
        assert S + E + I + R == pytest.approx(1000.0, abs=0.01)

    def test_exposed_compartment(self):
        S, E, I, R = seir_step(990.0, 0.0, 10.0, 0.0, 0.5, 0.2, 0.1, dt=1.0)
        # Some should move to exposed
        assert E > 0.0 or S < 990.0

    def test_no_infection_no_change(self):
        S, E, I, R = seir_step(1000.0, 0.0, 0.0, 0.0, 0.3, 0.2, 0.1, dt=1.0)
        assert S == pytest.approx(1000.0)


# ---------------------------------------------------------------------------
# sis_step
# ---------------------------------------------------------------------------


class TestSISStep:
    def test_conservation_of_population(self):
        S, I = sis_step(990.0, 10.0, 0.3, 0.1, dt=1.0)
        assert S + I == pytest.approx(1000.0, abs=0.01)

    def test_recovery_returns_to_susceptible(self):
        S, I = sis_step(500.0, 500.0, 0.0, 0.5, dt=1.0)
        # With no transmission, infected recover to susceptible
        assert S > 500.0
        assert I < 500.0

    def test_non_negative_values(self):
        S, I = sis_step(10.0, 990.0, 0.9, 0.01, dt=1.0)
        assert S >= 0
        assert I >= 0
