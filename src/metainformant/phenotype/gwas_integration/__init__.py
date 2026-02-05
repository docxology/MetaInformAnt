"""GWAS integration subpackage for phenotype analysis.

Provides phenome-wide association studies (PheWAS), phenotype correlation,
genetic risk score computation, and heritability screening.
"""

from __future__ import annotations

from .phewas import (
    run_phewas,
    phenotype_correlation_matrix,
    genetic_risk_score,
    phenotype_heritability_screen,
    categorize_phenotypes,
)

__all__ = [
    "run_phewas",
    "phenotype_correlation_matrix",
    "genetic_risk_score",
    "phenotype_heritability_screen",
    "categorize_phenotypes",
]
