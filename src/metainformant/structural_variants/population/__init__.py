"""Population-scale structural variant analysis submodule.

Provides population genotyping, allele frequency computation, association
testing, population structure analysis via PCA, LD computation between
SVs and SNPs, and multi-sample callset merging.
"""

from __future__ import annotations

from . import sv_population

from .sv_population import (
    genotype_sv_population,
    merge_sv_callsets,
    sv_allele_frequency,
    sv_association_test,
    sv_ld_analysis,
    sv_population_structure,
)

__all__ = [
    "sv_population",
    "genotype_sv_population",
    "merge_sv_callsets",
    "sv_allele_frequency",
    "sv_association_test",
    "sv_ld_analysis",
    "sv_population_structure",
]
