# Bayesian Inference

Bayesian inference methods for biological and statistical modelling, including MCMC sampling, Approximate Bayesian Computation, conjugate analysis, and model comparison criteria.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `inference` module |
| `inference.py` | Metropolis-Hastings, ABC, Bayes factors, conjugate priors, DIC, WAIC |

## Key Functions

| Function | Description |
|----------|-------------|
| `metropolis_hastings()` | MCMC sampling via Metropolis-Hastings algorithm |
| `abc_rejection()` | Approximate Bayesian Computation with rejection sampling |
| `compute_bayes_factor()` | Compute Bayes factor for model comparison |
| `conjugate_beta_binomial()` | Conjugate Beta-Binomial posterior analysis |
| `conjugate_normal()` | Conjugate Normal-Normal posterior analysis |
| `compute_dic()` | Deviance Information Criterion for model selection |
| `compute_waic()` | Widely Applicable Information Criterion for model selection |

## Usage

```python
from metainformant.math.bayesian import inference

samples = inference.metropolis_hastings(log_posterior, initial, n_samples=10000)
abc_result = inference.abc_rejection(simulator, observed, n_samples=1000)
bf = inference.compute_bayes_factor(model_a_evidence, model_b_evidence)
posterior = inference.conjugate_beta_binomial(alpha=1, beta=1, successes=7, trials=10)
dic = inference.compute_dic(log_likelihood_samples, theta_samples)
```
