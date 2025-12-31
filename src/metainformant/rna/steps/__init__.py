"""RNA-seq workflow steps for METAINFORMANT.

This module contains individual step implementations for the RNA-seq analysis
workflow, providing modular components that can be composed into complete
analysis pipelines.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


@dataclass
class StepResult:
    """Result of executing a workflow step.

    Attributes:
        success: Whether the step completed successfully
        outputs: Dictionary of output files/data produced by the step
        metadata: Additional metadata about the step execution
        error_message: Error message if step failed
        execution_time: Time taken to execute the step (seconds)
    """
    success: bool
    outputs: Dict[str, Any]
    metadata: Dict[str, Any]
    error_message: Optional[str] = None
    execution_time: Optional[float] = None


# Import all step modules for easy access
try:
    from . import metadata
except ImportError:
    metadata = None

try:
    from . import integrate
except ImportError:
    integrate = None

try:
    from . import config
except ImportError:
    config = None

try:
    from . import select
except ImportError:
    select = None

try:
    from . import getfastq
except ImportError:
    getfastq = None

try:
    from . import quant
except ImportError:
    quant = None

try:
    from . import merge
except ImportError:
    merge = None

try:
    from . import cstmm
except ImportError:
    cstmm = None

try:
    from . import curate
except ImportError:
    curate = None

try:
    from . import csca
except ImportError:
    csca = None

try:
    from . import sanity
except ImportError:
    sanity = None

try:
    from . import process_samples
except ImportError:
    process_samples = None

try:
    from . import download_progress
except ImportError:
    download_progress = None


__all__ = [
    "StepResult",
    "metadata",
    "integrate",
    "config",
    "select",
    "getfastq",
    "quant",
    "merge",
    "cstmm",
    "curate",
    "csca",
    "sanity",
    "process_samples",
    "download_progress",
]
