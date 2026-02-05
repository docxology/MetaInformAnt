"""Metabolizer status prediction subpackage for pharmacogenomics.

Provides CPIC-style metabolizer phenotype prediction from genotype data,
activity score computation, and dose adjustment recommendations.
"""

from __future__ import annotations

from .metabolizer_status import (
    predict_metabolizer_status,
    compute_activity_score,
    classify_metabolizer,
    dose_adjustment,
    default_allele_function_table,
)

__all__ = [
    "predict_metabolizer_status",
    "compute_activity_score",
    "classify_metabolizer",
    "dose_adjustment",
    "default_allele_function_table",
]
