"""RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.

This module provides comprehensive tools for RNA-seq analysis, including
amalgkit workflow orchestration, expression quantification, quality control,
differential expression analysis, multi-species comparative studies,
alternative splicing detection and isoform quantification, and
cell type deconvolution from bulk RNA-seq data.
"""

from __future__ import annotations

# Import subpackages
from . import amalgkit  # Amalgkit integration
from . import analysis  # Analysis modules (expression, QC, validation)
from . import core  # Base configs and utils
from . import deconvolution  # Cell type deconvolution from bulk RNA-seq
from . import engine  # Orchestration engine (depends on others)
from . import retrieval  # Data retrieval (ENA downloader)
from . import splicing  # Alternative splicing detection and isoform analysis
from . import steps  # Step runners shim
from .amalgkit import amalgkit as amalgkit_module
from .amalgkit import (
    genome_prep,
    metadata_filter,
    tissue_normalizer,
)
from .analysis import (
    protein_integration,
    validation,
)

# Import modules from subpackages for backward compatibility
from .core import (
    cleanup,
    configs,
    deps,
    environment,
    validate_environment,
)

# Direct imports of commonly used classes and functions
from .core.configs import AmalgkitRunLayout, RNAPipelineConfig
from .engine import (
    discovery,
    monitoring,
    orchestration,
    pipeline,
    progress_tracker,
    workflow,
)
from .engine.pipeline import summarize_curate_tables
from .engine.progress_tracker import ProgressTracker
from .engine.workflow import AmalgkitWorkflowConfig

# Conditional imports for new analysis modules
try:
    from .analysis.expression import (
        differential_expression,
        filter_low_expression,
        normalize_counts,
        pca_analysis,
    )
except ImportError:
    pass

try:
    from .analysis.qc import (
        compute_sample_metrics,
        generate_qc_report,
    )
except ImportError:
    pass

# Conditional imports for splicing submodule
try:
    from .splicing import (
        build_isoform_graph,
        classify_splicing_events,
        compare_isoform_usage,
        compute_isoform_diversity,
        compute_psi,
        compute_splice_site_strength,
        detect_splice_junctions,
        differential_splicing,
        enumerate_isoforms,
        find_novel_junctions,
        quantify_isoforms,
    )
except ImportError:
    pass

# Conditional imports for deconvolution submodule
try:
    from .deconvolution import (
        batch_deconvolve,
        build_signature_matrix,
        deconvolve_nnls,
        deconvolve_svr,
        select_marker_genes,
        validate_deconvolution,
    )
except ImportError:
    pass

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "amalgkit",
    "analysis",
    "core",
    "engine",
    "retrieval",
    "splicing",
    "deconvolution",
    # Core workflow management
    "workflow",
    "orchestration",
    "pipeline",
    "monitoring",
    "progress_tracker",
    # Amalgkit integration
    "amalgkit_module",
    "tissue_normalizer",
    # Configuration and setup
    "configs",
    # Environment and discovery
    "discovery",
    "environment",
    "validate_environment",
    "genome_prep",
    "protein_integration",
    # Utilities
    "cleanup",
    "deps",
    "metadata_filter",
    "validation",
    # Direct exports
    "RNAPipelineConfig",
    "AmalgkitRunLayout",
    "ProgressTracker",
    "summarize_curate_tables",
    "AmalgkitWorkflowConfig",
    "steps",
    # Expression analysis
    "normalize_counts",
    "differential_expression",
    "pca_analysis",
    "filter_low_expression",
    # QC
    "compute_sample_metrics",
    "generate_qc_report",
    # Splicing analysis
    "detect_splice_junctions",
    "classify_splicing_events",
    "compute_psi",
    "differential_splicing",
    "find_novel_junctions",
    "compute_splice_site_strength",
    "quantify_isoforms",
    "build_isoform_graph",
    "enumerate_isoforms",
    "compute_isoform_diversity",
    "compare_isoform_usage",
    # Deconvolution
    "deconvolve_nnls",
    "deconvolve_svr",
    "build_signature_matrix",
    "select_marker_genes",
    "validate_deconvolution",
    "batch_deconvolve",
]
