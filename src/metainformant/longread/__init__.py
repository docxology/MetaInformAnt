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

# Import subpackages
from . import io
from . import quality
from . import analysis
from . import assembly
from . import visualization
from . import utils
from . import workflow

# IO
from .io.fast5 import (
    read_fast5,
    extract_signal,
    extract_basecalls,
    get_read_metadata,
)
from .io.bam import (
    read_long_read_bam,
    extract_methylation_tags,
    get_supplementary_alignments,
    calculate_alignment_stats,
)
from .io.formats import (
    fast5_to_fastq,
    convert_pod5_to_fast5,
    write_paf,
)

# Quality
from .quality.metrics import (
    calculate_n50,
    calculate_nx,
    read_length_stats,
    quality_score_distribution,
    estimate_accuracy,
    calculate_throughput,
)
from .quality.filtering import (
    filter_by_length,
    filter_by_quality,
    trim_adapters,
    detect_adapters,
    split_chimeric_reads,
)

# Analysis
from .analysis.modified_bases import (
    detect_methylation,
    call_5mc,
    call_6ma,
    aggregate_methylation,
    differential_methylation,
)
from .analysis.structural import (
    detect_sv_from_long_reads,
    detect_insertions,
    detect_inversions,
    phase_structural_variants,
)
from .analysis.phasing import (
    phase_reads,
    build_haplotype_blocks,
    tag_reads_by_haplotype,
    calculate_phase_block_stats,
)

# Assembly
from .assembly.overlap import (
    find_overlaps,
    minimizer_sketch,
    compute_overlap_graph,
    filter_contained_reads,
)
from .assembly.consensus import (
    generate_consensus,
    polish_consensus,
    multiple_sequence_alignment,
    calculate_consensus_quality,
)
from .assembly.hybrid import (
    hybrid_assemble,
    correct_with_short_reads,
    scaffold_with_long_reads,
)

# Visualization
from .visualization.plots import (
    plot_read_length_histogram,
    plot_quality_vs_length,
    plot_dotplot,
    plot_alignment_view,
    plot_methylation_track,
    plot_phasing_blocks,
)

# Utils
from .utils.batch import (
    BatchResult,
    process_batch,
    batch_filter_reads,
    batch_compute_metrics,
)
from .utils.summary import (
    RunSummary,
    generate_qc_summary,
    generate_assembly_summary,
    generate_methylation_summary,
    generate_sv_summary,
    build_run_summary,
    export_run_summary,
    compare_run_summaries,
)

# Workflow
from .workflow.orchestrator import (
    LongReadOrchestrator,
    PipelineStep,
    PipelineResult,
)
from .workflow.reporting import (
    QCReport,
)

# Type checking imports
from typing import TYPE_CHECKING

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
]
