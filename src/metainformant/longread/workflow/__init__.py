"""Long-read workflow orchestration module.

Provides end-to-end pipeline orchestration for long-read sequencing analysis
including QC, assembly, methylation, and structural variant calling pipelines.
Ties together the io, quality, analysis, assembly, and visualization subpackages
into cohesive, configurable workflows.

Submodules:
- orchestrator: Pipeline execution engine with dependency resolution
- pipelines: Pre-defined pipeline configurations
- reporting: QC and analysis report generation
"""

from __future__ import annotations

from . import orchestrator
from . import pipelines
from . import reporting

from .orchestrator import LongReadOrchestrator, PipelineStep, PipelineResult
from .pipelines import (
    get_qc_pipeline_config,
    get_assembly_pipeline_config,
    get_methylation_pipeline_config,
    get_sv_pipeline_config,
    load_pipeline_config,
    validate_pipeline_config,
)
from .reporting import (
    QCReport,
    generate_qc_report,
    generate_assembly_report,
    generate_methylation_report,
    export_report,
    generate_run_summary,
)

__all__ = [
    # Submodules
    "orchestrator",
    "pipelines",
    "reporting",
    # Orchestrator
    "LongReadOrchestrator",
    "PipelineStep",
    "PipelineResult",
    # Pipeline configs
    "get_qc_pipeline_config",
    "get_assembly_pipeline_config",
    "get_methylation_pipeline_config",
    "get_sv_pipeline_config",
    "load_pipeline_config",
    "validate_pipeline_config",
    # Reporting
    "QCReport",
    "generate_qc_report",
    "generate_assembly_report",
    "generate_methylation_report",
    "export_report",
    "generate_run_summary",
]
