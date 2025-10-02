"""Quality control tools for bioinformatics data.

This module provides comprehensive quality control analysis for various
types of bioinformatics data including sequencing reads, expression data,
and genomic features. All implementations use real computational methods
without mocking.

Key Features:
- FASTQ quality assessment (FastQC-like functionality)
- Sequencing quality metrics and statistics
- Contamination detection and filtering
- Adapter detection and trimming recommendations
- Quality score distributions and per-base quality
- GC content analysis and bias detection
- Duplication analysis
- Length distribution analysis
- Overrepresented sequence detection
- Quality control reporting and visualization

Real Implementation Policy:
All functions perform actual quality control computations without mocking.
External tools like FastQC may be used when available, with graceful
degradation to internal implementations when not.
"""

from __future__ import annotations

# Core quality control imports (only import existing modules)
try:
    from .fastq import (
        adapter_content,
        analyze_fastq_quality,
        duplication_levels,
        gc_content_distribution,
        overrepresented_sequences,
        per_base_quality,
        per_sequence_quality,
        sequence_length_distribution,
    )
except ImportError:
    pass

# Import implemented modules
from .contamination import (
    detect_adapter_contamination,
    detect_cross_species_contamination,
    detect_mycoplasma_contamination,
    detect_rrna_contamination,
    detect_vector_contamination,
    generate_contamination_report,
)

from .metrics import (
    calculate_complexity_metrics,
    calculate_coverage_metrics,
    calculate_duplication_metrics,
    calculate_gc_metrics,
    calculate_length_metrics,
    calculate_quality_metrics,
    generate_quality_report,
)

__all__ = [
    # FASTQ analysis (implemented)
    "analyze_fastq_quality",
    "per_base_quality",
    "per_sequence_quality",
    "sequence_length_distribution",
    "gc_content_distribution",
    "adapter_content",
    "overrepresented_sequences",
    "duplication_levels",

    # Contamination detection (implemented)
    "detect_cross_species_contamination",
    "detect_rrna_contamination",
    "detect_mycoplasma_contamination",
    "detect_adapter_contamination",
    "detect_vector_contamination",
    "generate_contamination_report",

    # Quality metrics (implemented)
    "calculate_quality_metrics",
    "calculate_gc_metrics",
    "calculate_length_metrics",
    "calculate_duplication_metrics",
    "calculate_complexity_metrics",
    "calculate_coverage_metrics",
    "generate_quality_report",
]
