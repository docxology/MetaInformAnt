"""Metagenomics analysis module for METAINFORMANT.

This module provides comprehensive tools for metagenomic analysis, including
amplicon-based community profiling (16S/ITS), shotgun metagenomics (assembly,
binning, profiling), functional annotation (gene prediction, pathway reconstruction),
specialized visualization for microbial ecology data, community diversity metrics,
and comparative/differential abundance analysis.

Subpackages:
    amplicon: OTU clustering, ASV denoising, taxonomic classification
    shotgun: Metagenome assembly, binning, community profiling
    functional: Gene annotation, ORF prediction, pathway analysis
    visualization: Krona charts, stacked bars, rarefaction, ordination plots
    diversity: Alpha/beta diversity, rarefaction, PERMANOVA, ordination
    comparative: Differential abundance, indicator species, biomarker discovery
"""

from __future__ import annotations

# Type checking imports
from typing import TYPE_CHECKING

# Import new subpackages
# Import subpackages
from . import amplicon, comparative, diversity, functional, shotgun, visualization
from .amplicon.asv_denoising import denoise_sequences, estimate_error_rates, merge_paired_reads

# Re-export key classes and functions for convenience
from .amplicon.otu_clustering import calculate_identity, cluster_otus, filter_chimeras
from .amplicon.taxonomy import build_taxonomy_tree, calculate_confidence, classify_taxonomy
from .comparative.differential_abundance import (
    biomarker_discovery,
    clr_transform,
    differential_abundance,
    effect_size_analysis,
    indicator_species,
)

# Re-export new submodule functions for convenience
from .diversity.metrics import (
    alpha_diversity,
    beta_diversity,
    ordination,
    permanova,
    rarefaction_curve,
    rarefy,
)
from .functional.annotation import annotate_genes, classify_gene_families, predict_orfs
from .functional.pathways import calculate_pathway_completeness, compare_pathway_profiles, reconstruct_pathways
from .shotgun.assembly import assemble_contigs, calculate_assembly_stats, scaffold_contigs
from .shotgun.binning import assess_bin_quality, bin_contigs, calculate_tetranucleotide_freq, refine_bins
from .shotgun.profiling import build_kmer_index, calculate_relative_abundance, profile_community
from .visualization.plots import (
    plot_alpha_diversity,
    plot_heatmap,
    plot_krona_chart,
    plot_ordination,
    plot_rarefaction_curves,
    plot_stacked_bar,
)

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "amplicon",
    "shotgun",
    "functional",
    "visualization",
    "diversity",
    "comparative",
    # Amplicon
    "cluster_otus",
    "calculate_identity",
    "filter_chimeras",
    "denoise_sequences",
    "estimate_error_rates",
    "merge_paired_reads",
    "classify_taxonomy",
    "build_taxonomy_tree",
    "calculate_confidence",
    # Shotgun
    "assemble_contigs",
    "scaffold_contigs",
    "calculate_assembly_stats",
    "bin_contigs",
    "calculate_tetranucleotide_freq",
    "refine_bins",
    "assess_bin_quality",
    "profile_community",
    "calculate_relative_abundance",
    "build_kmer_index",
    # Functional
    "annotate_genes",
    "predict_orfs",
    "classify_gene_families",
    "reconstruct_pathways",
    "calculate_pathway_completeness",
    "compare_pathway_profiles",
    # Visualization
    "plot_krona_chart",
    "plot_stacked_bar",
    "plot_rarefaction_curves",
    "plot_ordination",
    "plot_alpha_diversity",
    "plot_heatmap",
    # Diversity
    "alpha_diversity",
    "beta_diversity",
    "rarefaction_curve",
    "rarefy",
    "permanova",
    "ordination",
    # Comparative
    "differential_abundance",
    "clr_transform",
    "indicator_species",
    "effect_size_analysis",
    "biomarker_discovery",
]
