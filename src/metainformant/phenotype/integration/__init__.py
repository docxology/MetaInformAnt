"""Phenotype integration with other -omics modules."""

from .cross_omic import (
    phenotype_genotype_association,
    trait_expression_correlation,
    multi_phenotype_integration,
    phenotype_environment_interaction,
)

__all__ = [
    "phenotype_genotype_association",
    "trait_expression_correlation",
    "multi_phenotype_integration",
    "phenotype_environment_interaction",
]
