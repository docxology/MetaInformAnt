"""Phenotype integration with other -omics modules."""

from .cross_omic import (
    multi_phenotype_integration,
    phenotype_environment_interaction,
    phenotype_genotype_association,
    trait_expression_correlation,
)

__all__ = [
    "phenotype_genotype_association",
    "trait_expression_correlation",
    "multi_phenotype_integration",
    "phenotype_environment_interaction",
]
