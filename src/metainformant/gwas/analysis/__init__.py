"""GWAS analysis modules."""

from metainformant.gwas.analysis.association import (
    association_test_linear,
    association_test_logistic,
)
from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction
from metainformant.gwas.analysis.utils import compute_r_squared

__all__ = [
    "association_test_linear",
    "association_test_logistic",
    "bonferroni_correction",
    "compute_r_squared",
    "fdr_correction",
]
