"""Long-read sequencing analysis module for METAINFORMANT.

This module provides comprehensive tools for analyzing long-read sequencing data
from Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) platforms,
including signal processing, quality control, modified base detection, structural
variant analysis, haplotype phasing, and de novo assembly.

Submodules:
- io: FAST5/POD5/BAM reading and format conversion
- quality: Read quality metrics and filtering
- analysis: Modified base detection, SV calling, haplotype phasing
- assembly: Overlap computation, consensus generation, hybrid assembly
- visualization: Read length histograms, dotplots, alignment views
- utils: Batch processing, run summaries, report generation
- workflow: End-to-end pipeline orchestration, configs, and reporting
"""

from __future__ import annotations

# Type checking imports
from typing import TYPE_CHECKING

# Import subpackages
from . import analysis, assembly, io, methylation, phasing, quality, utils, visualization, workflow

# Analysis
from .analysis.modified_bases import (
    aggregate_methylation,
    call_5mc,
    call_6ma,
    detect_methylation,
    differential_methylation,
)
from .analysis.phasing import (
    build_haplotype_blocks,
    calculate_phase_block_stats,
    phase_reads,
    tag_reads_by_haplotype,
)
from .analysis.structural import (
    detect_insertions,
    detect_inversions,
    detect_sv_from_long_reads,
    phase_structural_variants,
)
from .assembly.consensus import (
    calculate_consensus_quality,
    generate_consensus,
    multiple_sequence_alignment,
    polish_consensus,
)
from .assembly.hybrid import (
    correct_with_short_reads,
    hybrid_assemble,
    scaffold_with_long_reads,
)

# Assembly
from .assembly.overlap import (
    compute_overlap_graph,
    filter_contained_reads,
    find_overlaps,
    minimizer_sketch,
)
from .io.bam import (
    calculate_alignment_stats,
    extract_methylation_tags,
    get_supplementary_alignments,
    read_long_read_bam,
)

# IO
from .io.fast5 import (
    extract_basecalls,
    extract_signal,
    get_read_metadata,
    read_fast5,
)
from .io.formats import (
    convert_pod5_to_fast5,
    fast5_to_fastq,
    write_paf,
)

# Methylation submodule
from .methylation import aggregate_methylation as meth_aggregate_methylation
from .methylation import (
    call_methylation_from_signal,
    compute_methylation_stats,
    detect_dmrs,
    methylation_pattern_analysis,
)

# Phasing submodule
from .phasing import (
    allele_specific_analysis,
)
from .phasing import build_phase_blocks as phasing_build_phase_blocks
from .phasing import (
    compute_switch_errors,
    haplotag_reads,
)
from .phasing import phase_reads as phasing_phase_reads
from .quality.filtering import (
    detect_adapters,
    filter_by_length,
    filter_by_quality,
    split_chimeric_reads,
    trim_adapters,
)

# Quality
from .quality.metrics import (
    calculate_n50,
    calculate_nx,
    calculate_throughput,
    estimate_accuracy,
    quality_score_distribution,
    read_length_stats,
)

# Utils
from .utils.batch import (
    BatchResult,
    batch_compute_metrics,
    batch_filter_reads,
    process_batch,
)
from .utils.summary import (
    RunSummary,
    build_run_summary,
    compare_run_summaries,
    export_run_summary,
    generate_assembly_summary,
    generate_methylation_summary,
    generate_qc_summary,
    generate_sv_summary,
)

# Visualization
from .visualization.plots import (
    plot_alignment_view,
    plot_dotplot,
    plot_methylation_track,
    plot_phasing_blocks,
    plot_quality_vs_length,
    plot_read_length_histogram,
)

# Workflow
from .workflow.orchestrator import (
    LongReadOrchestrator,
    PipelineResult,
    PipelineStep,
)
from .workflow.reporting import (
    QCReport,
)

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "io",
    "quality",
    "analysis",
    "assembly",
    "visualization",
    "utils",
    # IO functions
    "read_fast5",
    "extract_signal",
    "extract_basecalls",
    "get_read_metadata",
    "read_long_read_bam",
    "extract_methylation_tags",
    "get_supplementary_alignments",
    "calculate_alignment_stats",
    "fast5_to_fastq",
    "convert_pod5_to_fast5",
    "write_paf",
    # Quality functions
    "calculate_n50",
    "calculate_nx",
    "read_length_stats",
    "quality_score_distribution",
    "estimate_accuracy",
    "calculate_throughput",
    "filter_by_length",
    "filter_by_quality",
    "trim_adapters",
    "detect_adapters",
    "split_chimeric_reads",
    # Analysis functions
    "detect_methylation",
    "call_5mc",
    "call_6ma",
    "aggregate_methylation",
    "differential_methylation",
    "detect_sv_from_long_reads",
    "detect_insertions",
    "detect_inversions",
    "phase_structural_variants",
    "phase_reads",
    "build_haplotype_blocks",
    "tag_reads_by_haplotype",
    "calculate_phase_block_stats",
    # Assembly functions
    "find_overlaps",
    "minimizer_sketch",
    "compute_overlap_graph",
    "filter_contained_reads",
    "generate_consensus",
    "polish_consensus",
    "multiple_sequence_alignment",
    "calculate_consensus_quality",
    "hybrid_assemble",
    "correct_with_short_reads",
    "scaffold_with_long_reads",
    # Visualization functions
    "plot_read_length_histogram",
    "plot_quality_vs_length",
    "plot_dotplot",
    "plot_alignment_view",
    "plot_methylation_track",
    "plot_phasing_blocks",
    # Utils functions
    "BatchResult",
    "process_batch",
    "batch_filter_reads",
    "batch_compute_metrics",
    "RunSummary",
    "generate_qc_summary",
    "generate_assembly_summary",
    "generate_methylation_summary",
    "generate_sv_summary",
    "build_run_summary",
    "export_run_summary",
    "compare_run_summaries",
    # Workflow
    "workflow",
    "LongReadOrchestrator",
    "PipelineStep",
    "PipelineResult",
    "QCReport",
    # Methylation submodule
    "methylation",
    "meth_aggregate_methylation",
    "call_methylation_from_signal",
    "compute_methylation_stats",
    "detect_dmrs",
    "methylation_pattern_analysis",
    # Phasing submodule
    "phasing",
    "allele_specific_analysis",
    "phasing_build_phase_blocks",
    "compute_switch_errors",
    "haplotag_reads",
    "phasing_phase_reads",
]
