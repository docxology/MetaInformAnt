"""Information workflows and pipelines.

Higher-level orchestration for information-theoretic analysis workflows,
combining metrics, integration, and reporting.
"""

from __future__ import annotations

from .workflows import (
    batch_entropy_analysis,
    compare_datasets,
    information_report,
    information_workflow,
)

__all__ = [
    "batch_entropy_analysis",
    "compare_datasets",
    "information_report",
    "information_workflow",
]
