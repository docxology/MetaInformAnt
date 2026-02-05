"""Metagenomics analysis module for METAINFORMANT.

This module provides comprehensive tools for metagenomic analysis, including
amplicon-based community profiling (16S/ITS), shotgun metagenomics (assembly,
binning, profiling), functional annotation (gene prediction, pathway reconstruction),
and specialized visualization for microbial ecology data.

Subpackages:
    amplicon: OTU clustering, ASV denoising, taxonomic classification
    shotgun: Metagenome assembly, binning, community profiling
    functional: Gene annotation, ORF prediction, pathway analysis
    visualization: Krona charts, stacked bars, rarefaction, ordination plots
"""

from __future__ import annotations

# Import subpackages
from . import amplicon
from . import shotgun
from . import functional
from . import visualization

# Re-export key classes and functions for convenience
from .amplicon.otu_clustering import cluster_otus, calculate_identity, filter_chimeras
from .amplicon.asv_denoising import denoise_sequences, estimate_error_rates, merge_paired_reads
from .amplicon.taxonomy import classify_taxonomy, build_taxonomy_tree, calculate_confidence

from .shotgun.assembly import assemble_contigs, scaffold_contigs, calculate_assembly_stats
from .shotgun.binning import bin_contigs, calculate_tetranucleotide_freq, refine_bins, assess_bin_quality
from .shotgun.profiling import profile_community, calculate_relative_abundance, build_kmer_index

from .functional.annotation import annotate_genes, predict_orfs, classify_gene_families
from .functional.pathways import reconstruct_pathways, calculate_pathway_completeness, compare_pathway_profiles

from .visualization.plots import (
    plot_krona_chart,
    plot_stacked_bar,
    plot_rarefaction_curves,
    plot_ordination,
    plot_alpha_diversity,
    plot_heatmap,
)

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "amplicon",
    "shotgun",
    "functional",
    "visualization",
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
]
