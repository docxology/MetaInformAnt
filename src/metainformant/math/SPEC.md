# Math Module Technical Specification

## Overview

The `metainformant.math` module provides the mathematical backbone for theoretical biology models, including population genetics, epidemiology, and evolutionary game theory. It separates pure mathematical logic from data I/O and visualization.

## Architecture

The module is structured into domain-specific submodules:

- **`population_genetics`**: Implements theories of allele frequency change, drift, and coalescence.
    - *Key Algorithms*: Coalescent simulation (Kingman's), Tajima's D, Fst (Weir & Cockerham).
- **`epidemiology`**: Deterministic compartmental models.
    - *Models*: SIR, SEIR, SIS (ODEs).
- **`quantitative_genetics`**: Continuous trait evolution.
    - *Key Equations*: Price Equation, Breeder's Equation.
- **`evolutionary_dynamics`**: Time-series dynamics.
    - *Models*: Logistic Map (chaos), Lotka-Volterra (predator-prey), Replicator Dynamics.
- **`decision_theory`**: Behavioral modeling.
    - *Models*: Drift-Diffusion Model (DDM).
- **`core`**: Shared statistical and visualization utilities.

## Dependencies

- **Required**: `numpy` (>=1.20) for vectorization and array operations.
- **Optional**: `scipy` (>=1.7) for advanced optimization (DDM fitting) and statistical functions (digamma, special functions).
    - The module uses graceful degradation: if `scipy` is missing, slower or simplified implementations are used.

## Design Principles

1.  **Functional Purity**: Mathematical functions should generally be pure, taking parameters as input and returning results without side effects.
2.  **Vectorization**: Prefer `numpy` array operations over loops for performance.
3.  **Type Hints**: Use `Sequence[float] | np.ndarray` to accept flexible inputs.
4.  **Graceful Fallbacks**: Always provide a fallback for optional dependencies (e.g., `scipy`).
5.  **Triple Play Standards**: Every significant directory must have `README.md`, `AGENTS.md`, and `SPEC.md`.

## Numerical Stability

- **Log-space calculations**: Use log-probabilities for likelihoods to avoid underflow.
- **Bounds checking**: Validate inputs (e.g., probabilities in [0, 1], positive population sizes).

## Future Roadmap

- **Stochastic Differential Equations (SDEs)**: Add support for SDEs in epidemiology and dynamics.
- **Spatial Models**: Add spatial structure to population genetics (Stepping Stone model).
