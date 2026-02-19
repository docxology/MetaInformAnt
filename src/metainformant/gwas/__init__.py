"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, data, finemapping, heritability, visualization, workflow

# Convenience exports for common GWAS functions
from .analysis.association import association_test_linear, association_test_logistic
from .analysis.correction import bonferroni_correction, fdr_correction
from .analysis.heritability import estimate_heritability
from .analysis.quality import parse_vcf_full
from .analysis.structure import compute_kinship_matrix, compute_pca
from .data.metadata import (
    get_geographic_coordinates,
    get_population_labels,
    load_sample_metadata,
    validate_metadata,
)
from .visualization.config import THEMES, PlotStyle, apply_style, get_style
from .visualization.genomic.genome import genome_wide_ld_heatmap
from .visualization.genomic.ld import compute_ld_decay, ld_decay_plot
from .visualization.interactive.composite import (
    gwas_summary_panel,
    population_structure_panel,
)
from .visualization.interactive.finemapping import (
    compute_credible_set,
    credible_set_plot,
)
from .visualization.interactive.interactive import interactive_manhattan
from .visualization.interactive.phenotype import (
    genotype_phenotype_boxplot,
    phenotype_correlation_matrix,
    phenotype_distribution,
)
from .visualization.population.geography import (
    population_count_map,
    sample_map,
)

__all__ = [
    "analysis",
    "data",
    "finemapping",
    "heritability",
    "visualization",
    "workflow",
    # Analysis
    "association_test_linear",
    "association_test_logistic",
    "bonferroni_correction",
    "fdr_correction",
    "estimate_heritability",
    "parse_vcf_full",
    "compute_kinship_matrix",
    "compute_pca",
    # Data
    "get_geographic_coordinates",
    "get_population_labels",
    "load_sample_metadata",
    "validate_metadata",
    # Visualization
    "THEMES",
    "PlotStyle",
    "apply_style",
    "get_style",
    "genome_wide_ld_heatmap",
    "compute_ld_decay",
    "ld_decay_plot",
    "gwas_summary_panel",
    "population_structure_panel",
    "compute_credible_set",
    "credible_set_plot",
    "interactive_manhattan",
    "genotype_phenotype_boxplot",
    "phenotype_correlation_matrix",
    "phenotype_distribution",
    "population_count_map",
    "sample_map",
]
