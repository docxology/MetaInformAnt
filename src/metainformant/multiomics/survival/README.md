# Multi-Omic Survival Analysis

Cox proportional hazards regression, Kaplan-Meier estimation, log-rank tests, and regularized multi-omic survival models with risk stratification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `analysis` module |
| `analysis.py` | Cox PH, Kaplan-Meier, log-rank, Lasso-Cox, concordance index, risk stratification |

## Key Functions

| Function | Description |
|----------|-------------|
| `cox_regression()` | Cox proportional hazards regression via Newton-Raphson with Breslow ties |
| `kaplan_meier()` | Product-limit survival curve estimator with Greenwood variance |
| `log_rank_test()` | Mantel-Haenszel chi-squared test comparing survival curves |
| `multi_omic_survival_model()` | L1-penalized Lasso-Cox model for multi-omic feature selection |
| `risk_stratification()` | Stratify subjects into risk groups from Cox model predictions |
| `compute_concordance_index()` | Harrell's C-index for discrimination evaluation |

## Usage

```python
from metainformant.multiomics.survival import analysis

cox = analysis.cox_regression(time, event, covariates, covariate_names=names)
km = analysis.kaplan_meier(time, event)
lr = analysis.log_rank_test(time, event, groups)
model = analysis.multi_omic_survival_model(time, event, omic_features)
c_index = analysis.compute_concordance_index(predicted_risk, time, event)
```
