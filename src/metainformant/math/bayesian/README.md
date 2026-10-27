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
| `metropolis_hastings()` | MCMC sampling via Metropolis-Hastings algorithm; keyword-only `random_state` makes a chain reproducible |
| `abc_rejection()` | Approximate Bayesian Computation with rejection sampling; needs a `prior_sampler` and a `simulator` |
| `compute_bayes_factor()` | Bayes factor and Jeffreys-scale reading from two log marginal likelihoods |
| `conjugate_beta_binomial()` | Conjugate Beta-Binomial posterior analysis (`prior_alpha` and `prior_beta` set the prior; the 95% interval uses a normal approximation to the Beta posterior) |
| `conjugate_normal()` | Conjugate Normal-Normal posterior analysis |
| `compute_dic()` | Deviance Information Criterion for model selection (p_D is estimated as half the variance of the deviance) |
| `compute_waic()` | Widely Applicable Information Criterion for model selection |

## Usage

```python
import random

from metainformant.math.bayesian import inference


data = [0.4, 0.7, 0.5]


def log_likelihood(mu: float) -> float:
    return sum(-0.5 * (x - mu) ** 2 for x in data)  # unit-variance normal, up to a constant


# MCMC (flat prior): an int seed (or a random.Random) as random_state reproduces the chain exactly
chain = inference.metropolis_hastings(
    lambda params: log_likelihood(params[0]), [0.0], n_samples=10000, burn_in=1000, random_state=42
)
log_likelihoods = [log_likelihood(sample[0]) for sample in chain["samples"]]
dic = inference.compute_dic(log_likelihoods, chain["samples"])

# ABC rejection: the prior sampler and the simulator draw the random numbers
rng = random.Random(1)
abc_result = inference.abc_rejection(
    simulator=lambda theta: [theta[0] + rng.gauss(0.0, 0.1)],
    observed_summary=[0.5],
    prior_sampler=lambda: [rng.uniform(-2.0, 2.0)],
    n_simulations=1000,
    tolerance=0.05,
)

bf = inference.compute_bayes_factor(log_marginal_1=-10.0, log_marginal_2=-12.5)
posterior = inference.conjugate_beta_binomial(successes=7, trials=10, prior_alpha=1.0, prior_beta=1.0)
```

## Reproducibility

- `metropolis_hastings(..., random_state=...)` accepts an `int` seed or a `random.Random`. An `int` seeds a private generator and a `random.Random` is used as given (and advances); in both cases Python's global `random` state is neither read nor changed, so the same seed always gives the same chain. The default `random_state=None` draws from the global generator, which `random.seed(...)` controls.
- `abc_rejection` draws no random numbers itself and has no `random_state`: seed `prior_sampler` and `simulator` instead.
- `compute_bayes_factor`, `conjugate_beta_binomial`, `conjugate_normal`, `compute_dic` and `compute_waic` are deterministic.
