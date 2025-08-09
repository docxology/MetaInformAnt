from __future__ import annotations


def sir_step(S: float, I: float, R: float, beta: float, gamma: float, dt: float = 0.01) -> tuple[float, float, float]:
    """Single Euler step of the SIR model.

    dS/dt = -beta S I, dI/dt = beta S I - gamma I, dR/dt = gamma I.
    Clamps negative values to 0.
    """
    S = max(0.0, S)
    I = max(0.0, I)
    R = max(0.0, R)
    dS = -beta * S * I
    dI = beta * S * I - gamma * I
    dR = gamma * I
    Sn = max(0.0, S + dt * dS)
    In = max(0.0, I + dt * dI)
    Rn = max(0.0, R + dt * dR)
    return Sn, In, Rn


def basic_reproduction_number(beta: float, gamma: float) -> float:
    """Basic reproduction number R0 for SIR: R0 = beta / gamma, guarded."""
    if gamma <= 0:
        return 0.0
    return max(0.0, beta / gamma)



def seir_step(
    S: float,
    E: float,
    I: float,
    R: float,
    beta: float,
    sigma: float,
    gamma: float,
    dt: float = 0.01,
) -> tuple[float, float, float, float]:
    """Single Euler step of the SEIR model.

    dS/dt = -beta S I, dE/dt = beta S I - sigma E, dI/dt = sigma E - gamma I, dR/dt = gamma I
    """
    S = max(0.0, S)
    E = max(0.0, E)
    I = max(0.0, I)
    R = max(0.0, R)
    dS = -beta * S * I
    dE = beta * S * I - sigma * E
    dI = sigma * E - gamma * I
    dR = gamma * I
    Sn = max(0.0, S + dt * dS)
    En = max(0.0, E + dt * dE)
    In = max(0.0, I + dt * dI)
    Rn = max(0.0, R + dt * dR)
    return Sn, En, In, Rn


def herd_immunity_threshold(R0: float) -> float:
    """Herd immunity threshold ≈ 1 - 1/R0 (non-negative and ≤1)."""
    if R0 <= 0:
        return 0.0
    val = 1.0 - 1.0 / R0
    return max(0.0, min(1.0, val))


def sis_step(S: float, I: float, beta: float, gamma: float, dt: float = 0.01) -> tuple[float, float]:
    """Single Euler step for SIS model (no removed class)."""
    S = max(0.0, S)
    I = max(0.0, I)
    dS = -beta * S * I + gamma * I
    dI = beta * S * I - gamma * I
    Sn = max(0.0, S + dt * dS)
    In = max(0.0, I + dt * dI)
    return Sn, In


def effective_reproduction_number(R0: float, susceptible_fraction: float) -> float:
    """Effective reproduction number R_e = R0 * S."""
    S = max(0.0, min(1.0, susceptible_fraction))
    return max(0.0, R0 * S)


