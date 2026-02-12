"""Network information flow subpackage.

Provides transfer entropy, Granger causality, network entropy,
and methods for constructing directed and undirected information
flow networks from multivariate time series."""
from __future__ import annotations

from . import information_flow

__all__ = ['information_flow']
