"""Pipeline orchestration engine for long-read analysis workflows.

This is a re-export module that aggregates the public API from:
- orchestrator_core: Core orchestrator class, dataclasses, and helpers
- pipeline_stages: Pipeline step builders for QC, assembly, methylation, SV

All classes and functions are re-exported here for backward compatibility.
"""

from __future__ import annotations

# Re-export core orchestrator components
from metainformant.longread.workflow.orchestrator_core import (
    LongReadOrchestrator,
    PipelineResult,
    PipelineStep,
    _extract_methylation_features,
    _write_json_safe,
)

# Re-export pipeline stage builders
from metainformant.longread.workflow.pipeline_stages import (
    _build_assembly_steps,
    _build_methylation_steps,
    _build_qc_steps,
    _build_sv_steps,
    _run_full_pipeline,
)

__all__ = [
    # Core
    "LongReadOrchestrator",
    "PipelineStep",
    "PipelineResult",
    "_extract_methylation_features",
    "_write_json_safe",
    # Pipeline stages
    "_build_qc_steps",
    "_build_assembly_steps",
    "_build_methylation_steps",
    "_build_sv_steps",
    "_run_full_pipeline",
]
