"""GWAS statistical fine-mapping submodule.

Provides credible set analysis, colocalization, SuSiE regression,
and conditional analysis for identifying causal variants from GWAS signals.
"""

from metainformant.gwas.finemapping.colocalization import (
    compute_clpp,
    eqtl_coloc,
    multi_trait_coloc,
    regional_coloc,
)
from metainformant.gwas.finemapping.credible_sets import (
    annotate_credible_set,
    colocalization,
    compute_bayes_factors,
    compute_credible_set,
    conditional_analysis,
    susie_regression,
)
from metainformant.gwas.finemapping.eqtl import (
    cis_eqtl_scan,
    conditional_eqtl,
    eqtl_effect_sizes,
    eqtl_summary_stats,
    trans_eqtl_scan,
)

__all__ = [
    # Credible set analysis
    "compute_credible_set",
    "susie_regression",
    "compute_bayes_factors",
    "colocalization",
    "conditional_analysis",
    "annotate_credible_set",
    # Multi-trait colocalization
    "multi_trait_coloc",
    "eqtl_coloc",
    "compute_clpp",
    "regional_coloc",
    # eQTL analysis
    "cis_eqtl_scan",
    "trans_eqtl_scan",
    "conditional_eqtl",
    "eqtl_effect_sizes",
    "eqtl_summary_stats",
]
