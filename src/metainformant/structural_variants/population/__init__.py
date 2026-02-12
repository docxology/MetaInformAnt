"""Population-scale structural variant analysis submodule.

Provides population genotyping, allele frequency computation, association
testing, population structure analysis via PCA, LD computation between
SVs and SNPs, and multi-sample callset merging."""
from __future__ import annotations

from . import sv_population

__all__ = ['sv_population']
