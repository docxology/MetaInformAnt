"""Core engine module for workflow orchestration.

This module provides:
- Generic pipeline primitives (``BasePipelineManager``, ``PipelinePhase``,
  ``PipelineItem``, ``Stage``) for building domain-agnostic multi-phase
  workflows with TUI visualization.
- RNA-seq backward-compatible types (``WorkflowManager``, ``SampleStage``,
  ``SampleState``) that wrap the generic pipeline for the download -> getfastq
  -> quant workflow.
"""

from metainformant.core.engine.workflow_manager import (  # Generic pipeline primitives; Backward-compatible RNA-seq types
    BasePipelineManager,
    PipelineItem,
    PipelinePhase,
    SampleStage,
    SampleState,
    Stage,
    WorkflowManager,
)

__all__ = [
    # Generic pipeline primitives
    "BasePipelineManager",
    "PipelineItem",
    "PipelinePhase",
    "Stage",
    # Backward-compatible RNA-seq types
    "WorkflowManager",
    "SampleStage",
    "SampleState",
]
