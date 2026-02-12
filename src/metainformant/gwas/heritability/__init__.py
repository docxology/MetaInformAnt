"""GWAS heritability estimation submodule.

Provides LD Score Regression, partitioned heritability, genetic correlation,
Haseman-Elston regression, GREML, and liability-scale transformations for
estimating and interpreting SNP heritability."""
from __future__ import annotations

from . import estimation

__all__ = ['estimation']
