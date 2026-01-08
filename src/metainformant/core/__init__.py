"""Core utilities for METAINFORMANT bioinformatics toolkit.

This module provides essential infrastructure components for the METAINFORMANT
bioinformatics analysis platform, including configuration management, I/O operations,
parallel processing, caching, and workflow orchestration.
"""

from __future__ import annotations

# Import subpackages (reordered for dependency resolution)
# 1. Utils (Foundation)
from . import utils
from .utils import (
    config,
    errors,
    hash,
    logging,
    optional_deps,
    progress,
    symbols,
    text,
)

# 2. Data (Structures)
from . import data
from .data import (
    validation,
)

# 3. IO (Input/Output)
# Depends on Utils for logging/config
from . import io
from .io import (
    cache,
    disk,
    download,
    io as io_module,
    paths,
)

# 4. Execution (Workflow)
# Depends on IO and Utils
from . import execution
from .execution import (
    discovery,
    parallel,
    workflow,
)

# Optional dependencies
try:
    from .data import db
except ImportError:
    db = None

# Core functions - always available
from .execution.workflow import create_sample_config, download_and_process_data, run_config_based_workflow, validate_config_file

# Type system
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

__all__ = [
    # Subpackages
    "data",
    "execution",
    "io",
    "utils",

    # Core modules (backward compatibility)
    "cache",
    "config",
    "discovery",
    "download",
    "errors",
    "hash",
    "io_module",
    "logging",
    "optional_deps",
    "parallel",
    "paths",
    "progress",
    "symbols",
    "text",
    "validation",
    "workflow",

    # Core functions
    "create_sample_config",
    "download_and_process_data",
    "run_config_based_workflow",
    "validate_config_file",

    # Optional modules
    "db",
    "disk",
]