"""Core engine module for workflow orchestration.

This module provides:
- Generic pipeline primitives (``BasePipelineManager``, ``PipelinePhase``,
  ``PipelineItem``, ``Stage``) for building domain-agnostic multi-phase
  workflows with TUI visualization.
- RNA-seq backward-compatible types (``WorkflowManager``, ``SampleStage``,
  ``SampleState``) that wrap the generic pipeline for the download -> getfastq
  -> quant workflow."""
from __future__ import annotations

from . import workflow_manager

__all__ = ['workflow_manager']
