"""Bayesian inference subpackage for mathematical biology.

Provides MCMC sampling (Metropolis-Hastings), Approximate Bayesian
Computation, Bayes factor computation, conjugate analysis, and
information criteria (DIC, WAIC).
"""

from __future__ import annotations

from .inference import (
    metropolis_hastings,
    abc_rejection,
    compute_bayes_factor,
    conjugate_beta_binomial,
    conjugate_normal,
    compute_dic,
    compute_waic,
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
