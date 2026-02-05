"""Bayesian inference subpackage for mathematical biology.

Provides MCMC sampling (Metropolis-Hastings), Approximate Bayesian
Computation, Bayes factor computation, conjugate analysis, and
information criteria (DIC, WAIC).
"""

from __future__ import annotations

from .inference import (
    abc_rejection,
    compute_bayes_factor,
    compute_dic,
    compute_waic,
    conjugate_beta_binomial,
    conjugate_normal,
    metropolis_hastings,
)

__all__ = [
    "metropolis_hastings",
    "abc_rejection",
    "compute_bayes_factor",
    "conjugate_beta_binomial",
    "conjugate_normal",
    "compute_dic",
    "compute_waic",
]
