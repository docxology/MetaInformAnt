# Multi-Omic Survival Analysis

Cox proportional hazards regression, Kaplan-Meier survival estimation,
log-rank tests, L1-regularised multi-omic survival models, risk
stratification, and concordance index computation.

## Key Concepts

**Cox proportional hazards** models the hazard as
`h_i(t) = h_0(t) * exp(X_i @ beta)`. Coefficients are estimated by
maximising the partial likelihood via Newton-Raphson; tied events use the
Breslow approximation.

**Kaplan-Meier** is the product-limit estimator for the survival function.
Confidence intervals use Greenwood's variance formula:
`Var(S(t)) = S(t)^2 * sum(d_i / (n_i * (n_i - d_i)))`.

**Log-rank test** (Mantel-Haenszel) compares observed vs expected events
across groups at each distinct event time and produces a chi-squared
statistic.

**Lasso-Cox** uses coordinate descent on the L1-penalised partial
log-likelihood for automatic feature selection from high-dimensional
multi-omic feature matrices.

**Concordance index** (Harrell's C) measures discrimination: the probability
that a subject who experienced an event earlier has a higher predicted risk
score. C = 0.5 is random; C = 1.0 is perfect.

## Function Reference

### cox_regression

```python
def cox_regression(
    time: list[float],
    event: list[int],
    covariates: list[list[float]] | Any,
    covariate_names: list[str] | None = None,
) -> dict[str, Any]
```

Fits a Cox PH model. Returns `coefficients`, `hazard_ratios` (exp(beta)),
`se` (standard errors), `p_values` (Wald test), `concordance_index`,
`log_likelihood`, and `covariate_names`. Requires numpy.

### kaplan_meier

```python
def kaplan_meier(
    time: list[float],
    event: list[int],
    groups: list[int] | None = None,
) -> dict[str, Any]
```

If `groups` is None, returns a single curve: `times`, `survival_prob`,
`confidence_lower`, `confidence_upper`, `n_at_risk`, `median_survival`.
If `groups` is provided, returns `{"groups": {group_label: curve_dict}}`.

### log_rank_test

```python
def log_rank_test(
    time: list[float],
    event: list[int],
    groups: list[int],
) -> dict[str, Any]
```

Mantel-Haenszel test comparing survival between >= 2 groups. Returns `chi2`,
`p_value`, `df`, and `observed_expected_per_group`.

### multi_omic_survival_model

```python
def multi_omic_survival_model(
    omic_features: dict[str, Any],
    time: list[float],
    event: list[int],
    method: str = "lasso_cox",
) -> dict[str, Any]
```

Concatenates multi-omic features, standardises, and fits a Lasso-Cox model
via coordinate descent. Lambda is set heuristically at `lambda_max / 10`.
Returns `selected_features` (non-zero coefficients), `coefficients`,
`c_index`, and `risk_scores`. Requires numpy.

### risk_stratification

```python
def risk_stratification(
    risk_scores: list[float],
    time: list[float],
    event: list[int],
    n_groups: int = 2,
) -> dict[str, Any]
```

Splits patients into risk groups by quantiles of predicted risk scores. Fits
Kaplan-Meier curves per group and computes a log-rank test. Returns
`group_labels`, `survival_curves`, `log_rank_p`, and `hazard_ratio`
(highest vs lowest group, from observed/expected ratios).

### compute_concordance_index

```python
def compute_concordance_index(
    risk_scores: list[float],
    time: list[float],
    event: list[int],
) -> dict[str, Any]
```

Harrell's C-index with Noether standard error estimate. Returns `c_index`,
`se`, `n_concordant`, `n_discordant`, `n_tied`.

## Usage Example

```python
import numpy as np
from metainformant.multiomics.survival.analysis import (
    cox_regression,
    kaplan_meier,
    log_rank_test,
    multi_omic_survival_model,
    risk_stratification,
)

# Cox regression
time = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
event = [1, 0, 1, 1, 0, 1]
covariates = [[0.5, 1.2], [1.0, 0.8], [0.3, 1.5], [0.7, 0.9], [1.2, 0.3], [0.9, 1.1]]
result = cox_regression(time, event, covariates, covariate_names=["age", "stage"])
print(f"C-index: {result['concordance_index']:.4f}")

# Kaplan-Meier with groups
km = kaplan_meier(time, event, groups=[0, 0, 0, 1, 1, 1])
for g, curve in km["groups"].items():
    print(f"Group {g}: median survival = {curve['median_survival']}")

# Multi-omic survival model
omic_data = {"rna": np.random.randn(50, 20), "protein": np.random.randn(50, 15)}
times = np.random.exponential(10, 50).tolist()
events = np.random.binomial(1, 0.6, 50).tolist()
model = multi_omic_survival_model(omic_data, times, events)
print(f"Selected features: {len(model['selected_features'])}/{35}")

# Risk stratification
strat = risk_stratification(model["risk_scores"], times, events, n_groups=3)
print(f"Log-rank p: {strat['log_rank_p']:.4e}, HR: {strat['hazard_ratio']:.2f}")
```

## Related Modules

- `metainformant.multiomics.methods` -- factorization and clustering
- `metainformant.multiomics.pathways` -- multi-omic pathway analysis
- `metainformant.gwas.analysis` -- genome-wide association analysis
