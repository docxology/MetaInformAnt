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
    from .metadata import run_metadata
except ImportError:
    metadata = None
    run_metadata = None

try:
    from . import integrate
    from .integrate import run_integrate
except ImportError:
    integrate = None
    run_integrate = None

try:
    from . import config
    from .config import run_config
except ImportError:
    config = None
    run_config = None

try:
    from . import select
    from .select import run_select
except ImportError:
    select = None
    run_select = None

try:
    from . import getfastq
    from .getfastq import run_getfastq
except ImportError:
    getfastq = None
    run_getfastq = None

try:
    from . import quant
    from .quant import run_quant
except ImportError:
    quant = None
    run_quant = None

try:
    from . import merge
    from .merge import run_merge
except ImportError:
    merge = None
    run_merge = None

try:
    from . import cstmm
    from .cstmm import run_cstmm
except ImportError:
    cstmm = None
    run_cstmm = None

try:
    from . import curate
    from .curate import run_curate
except ImportError:
    curate = None
    run_curate = None

try:
    from . import csca
    from .csca import run_csca
except ImportError:
    csca = None
    run_csca = None

try:
    from . import sanity
    from .sanity import run_sanity
except ImportError:
    sanity = None
    run_sanity = None

try:
    from . import process_samples
except ImportError:
    process_samples = None

try:
    from . import download_progress
except ImportError:
    download_progress = None


# Step runners registry
STEP_RUNNERS = {
    'metadata': metadata,
    'integrate': integrate,
    'config': config,
    'select': select,
    'getfastq': getfastq,
    'quant': quant,
    'merge': merge,
    'cstmm': cstmm,
    'curate': curate,
    'csca': csca,
    'sanity': sanity,
}


def _sample_already_quantified(sample_id: str, quant_dir: str) -> bool:
    """Check if a sample has already been quantified.

    Args:
        sample_id: Sample identifier
        quant_dir: Quantification output directory

    Returns:
        True if sample quantification files exist
    """
    import os
    from pathlib import Path

    quant_path = Path(quant_dir)
    if not quant_path.exists():
        return False

    # Check for common quantification output files
    expected_files = [
        f"{sample_id}_quant.sf",  # Salmon
        f"{sample_id}_abundance.tsv",  # Kallisto
        f"{sample_id}.genes.results",  # RSEM genes
        f"{sample_id}.isoforms.results",  # RSEM isoforms
    ]

    for expected_file in expected_files:
        if (quant_path / expected_file).exists():
            return True

    return False


__all__ = [
    "StepResult",
    "STEP_RUNNERS",
    "_sample_already_quantified",
    "metadata",
    "run_metadata",
    "integrate",
    "run_integrate",
    "config",
    "run_config",
    "select",
    "run_select",
    "getfastq",
    "run_getfastq",
    "quant",
    "run_quant",
    "merge",
    "run_merge",
    "cstmm",
    "run_cstmm",
    "curate",
    "run_curate",
    "csca",
    "run_csca",
    "sanity",
    "run_sanity",
    "process_samples",
    "download_progress",
]



