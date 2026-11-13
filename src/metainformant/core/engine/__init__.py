"""Core engine module for workflow orchestration.

This module provides generic pipeline primitives (``BasePipelineManager``,
``PipelinePhase``, ``PipelineItem``, ``Stage``) for building domain-agnostic
multi-phase workflows with TUI visualization."""

from __future__ import annotations

from . import workflow_manager

__all__ = ["workflow_manager"]
