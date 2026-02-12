"""Core information theory metrics.

Provides foundational information-theoretic measures: Shannon entropy,
continuous/differential entropy, and entropy estimation methods.
"""

from . import continuous, estimation, syntactic

__all__ = [
    "continuous",
    "estimation",
    "syntactic",
]
