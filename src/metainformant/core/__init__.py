"""Core utilities for METAINFORMANT bioinformatics toolkit.

This module provides essential infrastructure components for the METAINFORMANT
bioinformatics analysis platform, including configuration management, I/O operations,
parallel processing, caching, and workflow orchestration.
"""

from __future__ import annotations

# Core infrastructure - always available
from . import (
    cache,
    config,
    discovery,
    errors,
    hash,
    io,
    logging,
    download,
    optional_deps,
    parallel,
    paths,
    progress,
    symbols,
    text,
    validation,
    workflow,
)

# Core functions - always available
from .workflow import create_sample_config, download_and_process_data, run_config_based_workflow, validate_config_file

# Optional dependencies - gracefully handle missing components
try:
    from . import db
except ImportError:
    db = None

try:
    from . import disk
except ImportError:
    disk = None

try:
    from . import filesystem
except ImportError:
    filesystem = None

# Type system
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

__all__ = [
    # Core modules
    "cache",
    "config",
    "discovery",
    "download",
    "errors",
    "hash",
    "io",
    "logging",
    "optional_deps",
    "parallel",
    "paths",
    "progress",
    "symbols",
    "text",
    "validation",
    "workflow",

    # Core functions (from workflow)
    "create_sample_config",
    "download_and_process_data",
    "run_config_based_workflow",
    "validate_config_file",

    # Optional modules
    "db",
    "disk",
    "filesystem",
]