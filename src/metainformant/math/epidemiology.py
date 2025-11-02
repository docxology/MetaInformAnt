from __future__ import annotations


def sir_step(S: float, I: float, R: float, beta: float, gamma: float, dt: float = 0.01) -> tuple[float, float, float]:
    """Perform single Euler integration step of the SIR epidemiological model.
    
    The SIR model divides population into Susceptible (S), Infected (I), and
    Recovered (R) compartments. Used to model disease spread dynamics.
    
    Differential equations:
    - dS/dt = -beta * S * I
    - dI/dt = beta * S * I - gamma * I
    - dR/dt = gamma * I
    
    Args:
        S: Current susceptible population
        I: Current infected population
        R: Current recovered population
        beta: Transmission rate (contact rate * probability of transmission)
        gamma: Recovery rate (1 / average infectious period)
        dt: Time step size for Euler integration (default 0.01)
        
    Returns:
        Tuple of (next_S, next_I, next_R) after one time step.
        Values are clamped to non-negative.
        
    Examples:
        >>> S, I, R = sir_step(S=990.0, I=10.0, R=0.0, beta=0.3, gamma=0.1, dt=0.01)
        >>> S + I + R  # Total population should be preserved
        1000.0
        
    References:
        Kermack, W. O., & McKendrick, A. G. (1927). A contribution to the
        mathematical theory of epidemics. Proceedings of the royal society of
        London. Series A, 115(772), 700-721.
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
    """Calculate basic reproduction number R0 for SIR model.
    
    R0 represents the expected number of secondary infections from a single
    infected individual in a fully susceptible population. Critical threshold:
    R0 > 1 indicates epidemic potential, R0 < 1 indicates outbreak extinction.
    
    Args:
        beta: Transmission rate (contact rate * transmission probability)
        gamma: Recovery rate (1 / average infectious period)
        
    Returns:
        Basic reproduction number R0 = beta / gamma. Returns 0.0 if gamma <= 0.
        
    Examples:
        >>> basic_reproduction_number(beta=0.3, gamma=0.1)
        3.0
        >>> basic_reproduction_number(beta=0.1, gamma=0.2)
        0.5
        
    References:
        Diekmann, O., Heesterbeek, J. A. P., & Roberts, M. G. (2010). The construction
        of next-generation matrices for compartmental epidemic models. Journal of the
        Royal Society Interface, 7(47), 873-885.
    """
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
