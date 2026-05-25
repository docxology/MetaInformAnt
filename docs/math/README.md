# MATH

## Overview
Mathematical biology and theoretical modeling module for METAINFORMANT.

## Contents
- **[core/](index.md)**
- **[decision_theory/](ddm.md)**
- **[epidemiology/](epidemiology.md)**
- **[evolutionary_dynamics/](dynamics.md)**
- **perception/**
- **[population_genetics/](popgen.md)**
- **[quantitative_genetics/](popgen_stats.md)**
- **bayesian/** — Bayesian inference methods

## Structure

```mermaid
graph TD
    subgraph "Math Module"
        PG[population_genetics/] --> |coalescent.py| CO[Coalescent Theory]
        PG --> |demography.py| DM[Growth Models]
        PG --> |fst.py, ld.py| POP[Fst, LD, Selection]
        PG --> |effective_size.py| NE[Effective Population Size]

        ED[evolutionary_dynamics/] --> |core.py| LV[Lotka-Volterra, Logistic Map]
        ED --> |egt.py| EGT[Replicator Dynamics]

        QG[quantitative_genetics/] --> |core.py| H2[Heritability, Lande Equation]
        QG --> |price.py| PR[Price Equation]

        EP[epidemiology/] --> |models.py| SIR[SIR/SIS, R0, Herd Immunity]

        DT[decision_theory/] --> |ddm.py| DDM[Drift-Diffusion Models]

        BA[bayesian/] --> |inference.py| MCMC[MCMC, ABC, Bayes Factors]

        PE[perception/] --> PSY[Weber-Fechner, Signal Detection]

        CR[core/] --> |utilities.py| UT[Correlation, Regression]
    end
```

## Usage
Import module:
```python-snippet
from metainformant.math import ...
```
