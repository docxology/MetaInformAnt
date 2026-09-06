"""Quality control data visualization functions.

Aggregates the quality-control plotting surface across all biological data types
(sequencing reads, genomic variants, protein structures, single-cell data, and
multi-omics layers) under a single namespace.

Implementations live in the focused submodules and are re-exported here:

- ``quality_sequencing``: read-level QC (quality metrics, GC, lengths, adapters,
  duplication, overrepresented sequences, k-mers)
- ``quality_omics``: VCF, single-cell, protein structure, and multi-omics QC
- ``quality_assessment``: coverage uniformity, error profiles, batch effects,
  data integrity
"""

from __future__ import annotations

from .quality_assessment import (
    plot_batch_effects_qc,
    plot_coverage_uniformity,
    plot_data_integrity_metrics,
    plot_error_profiles,
)
from .quality_omics import (
    plot_multiomics_quality_overview,
    plot_protein_structure_quality,
    plot_singlecell_qc_metrics,
    plot_vcf_quality_metrics,
)
from .quality_sequencing import (
    plot_adapter_content,
    plot_gc_distribution,
    plot_kmer_profiles,
    plot_length_distribution,
    plot_overrepresented_sequences,
    plot_per_base_quality_boxplot,
    plot_quality_metrics,
    plot_sequence_duplication_levels,
)

__all__ = [
    # sequencing
    "plot_quality_metrics",
    "plot_adapter_content",
    "plot_gc_distribution",
    "plot_length_distribution",
    "plot_per_base_quality_boxplot",
    "plot_sequence_duplication_levels",
    "plot_overrepresented_sequences",
    "plot_kmer_profiles",
    # omics
    "plot_vcf_quality_metrics",
    "plot_singlecell_qc_metrics",
    "plot_protein_structure_quality",
    "plot_multiomics_quality_overview",
    # assessment
    "plot_coverage_uniformity",
    "plot_error_profiles",
    "plot_batch_effects_qc",
    "plot_data_integrity_metrics",
]
