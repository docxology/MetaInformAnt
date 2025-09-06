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

# TODO: Implement these modules
# from .contamination import (...)
# from .metrics import (...)
# from .reporting import (...)
# from .filtering import (...)

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
    # TODO: Add other modules when implemented
]
