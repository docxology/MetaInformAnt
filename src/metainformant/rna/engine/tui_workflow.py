"""RNA-owned TUI workflow manager.

The generic pipeline manager lives in :mod:`metainformant.core.engine`, but the
actual amalgkit runner is domain-specific.  This wrapper injects the RNA runner
so core does not depend on :mod:`metainformant.rna`.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core.engine.workflow_manager import (
    BasePipelineManager,
    PipelineItem,
    PipelinePhase,
    SampleStage,
    SampleState,
    Stage,
)
from metainformant.core.engine.workflow_manager import WorkflowManager as _CoreWorkflowManager
from metainformant.rna.amalgkit.amalgkit import run_amalgkit


class WorkflowManager(_CoreWorkflowManager):
    """RNA-seq workflow manager with the real amalgkit runner injected."""

    def __init__(self, config_path: Path, max_threads: int = 5) -> None:
        super().__init__(config_path=config_path, max_threads=max_threads, run_amalgkit=run_amalgkit)


__all__ = [
    "BasePipelineManager",
    "PipelineItem",
    "PipelinePhase",
    "SampleStage",
    "SampleState",
    "Stage",
    "WorkflowManager",
]
