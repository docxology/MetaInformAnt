"""Epidemiology and disease modeling functions.

This module provides mathematical functions for epidemiological modeling and analysis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def basic_reproduction_number(transmission_rate: float, recovery_rate: float, contact_rate: float = 0.0) -> float:
    """Calculate basic reproduction number (R₀) for SIR model.

    Args:
        transmission_rate: Rate of transmission (β)
        recovery_rate: Rate of recovery (γ)
        contact_rate: Contact rate (c) - not used in simple SIR

    Returns:
        Basic reproduction number (R₀ = β/γ)
    """
    return transmission_rate / recovery_rate


def effective_reproduction_number(R0: float, susceptible_fraction: float) -> float:
    """Calculate effective reproduction number (Rₑ).

    Args:
        R0: Basic reproduction number
        susceptible_fraction: Fraction of population susceptible

    Returns:
        Effective reproduction number (Rₑ = R₀ * S/N)
    """
    return R0 * susceptible_fraction


def herd_immunity_threshold(R0: float) -> float:
    """Calculate herd immunity threshold.

    Args:
        R0: Basic reproduction number

    Returns:
        Herd immunity threshold (1 - 1/R₀)
    """
    return 1 - 1 / R0


def sir_step(
    susceptible: float,
    infected: float,
    recovered: float,
    transmission_rate: float,
    recovery_rate: float,
    dt: float = 1.0,
) -> tuple[float, float, float]:
    """Single time step of SIR model.

    Args:
        susceptible: Number of susceptible individuals
        infected: Number of infected individuals
        recovered: Number of recovered individuals
        transmission_rate: Transmission rate (β)
        recovery_rate: Recovery rate (γ)
        dt: Time step size

    Returns:
        Tuple of (new_susceptible, new_infected, new_recovered)
    """
    total_population = susceptible + infected + recovered

    # SIR equations:
    # dS/dt = -β * S * I / N
    # dI/dt = β * S * I / N - γ * I
    # dR/dt = γ * I

    infection_rate = transmission_rate * susceptible * infected / total_population

    dS = -infection_rate
    dI = infection_rate - recovery_rate * infected
    dR = recovery_rate * infected

    new_susceptible = susceptible + dS * dt
    new_infected = infected + dI * dt
    new_recovered = recovered + dR * dt

    # Ensure non-negative values
    new_susceptible = max(0, new_susceptible)
    new_infected = max(0, new_infected)
    new_recovered = max(0, new_recovered)

    return new_susceptible, new_infected, new_recovered


def seir_step(
    susceptible: float,
    exposed: float,
    infected: float,
    recovered: float,
    transmission_rate: float,
    incubation_rate: float,
    recovery_rate: float,
    dt: float = 1.0,
) -> tuple[float, float, float, float]:
    """Single time step of SEIR model.

    Args:
        susceptible: Number of susceptible individuals
        exposed: Number of exposed individuals
        infected: Number of infected individuals
        recovered: Number of recovered individuals
        transmission_rate: Transmission rate (β)
        incubation_rate: Incubation rate (σ)
        recovery_rate: Recovery rate (γ)
        dt: Time step size

    Returns:
        Tuple of (new_susceptible, new_exposed, new_infected, new_recovered)
    """
    total_population = susceptible + exposed + infected + recovered

    # SEIR equations:
    # dS/dt = -β * S * I / N
    # dE/dt = β * S * I / N - σ * E
    # dI/dt = σ * E - γ * I
    # dR/dt = γ * I

    infection_rate = transmission_rate * susceptible * infected / total_population

    dS = -infection_rate
    dE = infection_rate - incubation_rate * exposed
    dI = incubation_rate * exposed - recovery_rate * infected
    dR = recovery_rate * infected

    new_susceptible = susceptible + dS * dt
    new_exposed = exposed + dE * dt
    new_infected = infected + dI * dt
    new_recovered = recovered + dR * dt

    # Ensure non-negative values
    new_susceptible = max(0, new_susceptible)
    new_exposed = max(0, new_exposed)
    new_infected = max(0, new_infected)
    new_recovered = max(0, new_recovered)

    return new_susceptible, new_exposed, new_infected, new_recovered


def sis_step(
    susceptible: float, infected: float, transmission_rate: float, recovery_rate: float, dt: float = 1.0
) -> tuple[float, float]:
    """Single time step of SIS model.

    Args:
        susceptible: Number of susceptible individuals
        infected: Number of infected individuals
        transmission_rate: Transmission rate (β)
        recovery_rate: Recovery rate (γ)
        dt: Time step size

    Returns:
        Tuple of (new_susceptible, new_infected)
    """
    total_population = susceptible + infected

    # SIS equations:
    # dS/dt = -β * S * I / N + γ * I
    # dI/dt = β * S * I / N - γ * I

    infection_rate = transmission_rate * susceptible * infected / total_population

    dS = -infection_rate + recovery_rate * infected
    dI = infection_rate - recovery_rate * infected

    new_susceptible = susceptible + dS * dt
    new_infected = infected + dI * dt

    # Ensure non-negative values
    new_susceptible = max(0, new_susceptible)
    new_infected = max(0, new_infected)

    return new_susceptible, new_infected
