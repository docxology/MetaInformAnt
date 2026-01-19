# SPEC: Math Scripts

Execution scripts for mathematical biology models, including population genetics theory and selection experiments.

## Core Models

### 1. Selection Dynamics
Scripts for simulating and fitting selection models.
- `fit_selection_parameters.py`: Uses optimization routines to estimate selection coefficients.

### 2. Decision Theory
Scripts for drift-diffusion models (DDM) and perception experiments.
- `run_ddm_analysis.py`: Fits DDM parameters to experimental data.

## Implementation Standards

- **Numerical Precision**: Use `numpy` and `scipy` for all intensive calculations.
- **Reporting**: Always output parameter estimates with associated confidence intervals or standard errors.
