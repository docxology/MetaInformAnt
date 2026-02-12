"""Bayesian inference subpackage for mathematical biology.

Provides MCMC sampling (Metropolis-Hastings), Approximate Bayesian
Computation, Bayes factor computation, conjugate analysis, and
information criteria (DIC, WAIC)."""
from __future__ import annotations

from . import inference

__all__ = ['inference']
