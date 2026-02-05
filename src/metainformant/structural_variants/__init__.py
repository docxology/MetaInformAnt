"""Structural Variant (SV) analysis module for METAINFORMANT.

This module provides comprehensive tools for detecting, annotating, filtering,
and visualizing structural variants including copy number variations (CNVs),
deletions, duplications, inversions, translocations, and insertions.

Subpackages:
    detection: SV calling, CNV detection, breakpoint refinement
    annotation: Gene overlap, regulatory overlap, functional impact
    filtering: Quality filtering, multi-caller merging, consensus
    visualization: Circos plots, coverage tracks, SV summaries

Config prefix: SV_
"""

from __future__ import annotations

# Type checking imports
from typing import TYPE_CHECKING

# Import subpackages
from . import annotation, detection, filtering, population, visualization

# Annotation - Functional impact
from .annotation.functional_impact import (
    assess_dosage_sensitivity,
    predict_functional_impact,
    predict_tad_disruption,
    score_pathogenicity,
)

# Annotation - Overlap
from .annotation.overlap import (
    annotate_gene_overlap,
    annotate_regulatory_overlap,
    calculate_overlap_fraction,
    find_nearest_gene,
)

# Detection - Breakpoints
from .detection.breakpoints import (
    calculate_breakpoint_confidence,
    cluster_breakpoints,
    detect_microhomology,
    refine_breakpoints,
)

# Detection - CNV
from .detection.cnv import (
    calculate_log2_ratio,
    call_cnv_states,
    detect_cnv_from_depth,
    merge_adjacent_segments,
    segment_coverage,
)

# Detection - SV calling
from .detection.sv_calling import (
    call_structural_variants,
    classify_sv_type,
    detect_discordant_pairs,
    detect_split_reads,
    genotype_sv,
)

# Filtering - Merge
from .filtering.merge import (
    calculate_reciprocal_overlap,
    deduplicate_variants,
    merge_callsets,
    survivor_merge,
)

# Filtering - Quality
from .filtering.quality_filter import (
    apply_blacklist,
    filter_by_frequency,
    filter_by_quality,
    filter_by_size,
)

# Population submodule
from .population import (
    genotype_sv_population,
    merge_sv_callsets,
    sv_allele_frequency,
    sv_association_test,
    sv_ld_analysis,
    sv_population_structure,
)

# Visualization
from .visualization.plots import (
    plot_breakpoint_detail,
    plot_circos,
    plot_cnv_profile,
    plot_coverage_track,
    plot_sv_size_distribution,
    plot_sv_type_summary,
)

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "detection",
    "annotation",
    "filtering",
    "visualization",
    # Detection - CNV
    "detect_cnv_from_depth",
    "segment_coverage",
    "call_cnv_states",
    "merge_adjacent_segments",
    "calculate_log2_ratio",
    # Detection - SV calling
    "call_structural_variants",
    "detect_split_reads",
    "detect_discordant_pairs",
    "classify_sv_type",
    "genotype_sv",
    # Detection - Breakpoints
    "refine_breakpoints",
    "detect_microhomology",
    "cluster_breakpoints",
    "calculate_breakpoint_confidence",
    # Annotation - Overlap
    "annotate_gene_overlap",
    "annotate_regulatory_overlap",
    "calculate_overlap_fraction",
    "find_nearest_gene",
    # Annotation - Functional impact
    "predict_functional_impact",
    "assess_dosage_sensitivity",
    "predict_tad_disruption",
    "score_pathogenicity",
    # Filtering - Quality
    "filter_by_quality",
    "filter_by_size",
    "filter_by_frequency",
    "apply_blacklist",
    # Filtering - Merge
    "merge_callsets",
    "calculate_reciprocal_overlap",
    "survivor_merge",
    "deduplicate_variants",
    # Visualization
    "plot_circos",
    "plot_coverage_track",
    "plot_sv_size_distribution",
    "plot_sv_type_summary",
    "plot_breakpoint_detail",
    "plot_cnv_profile",
    # Population submodule
    "population",
    "genotype_sv_population",
    "merge_sv_callsets",
    "sv_allele_frequency",
    "sv_association_test",
    "sv_ld_analysis",
    "sv_population_structure",
]
