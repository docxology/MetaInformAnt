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
    parallel,
    paths,
    progress,
    symbols,
    text,
    validation,
    workflow,
)

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
    "errors",
    "hash",
    "io",
    "logging",
    "parallel",
    "paths",
    "progress",
    "symbols",
    "text",
    "validation",
    "workflow",

    # Optional modules
    "db",
    "disk",
    "filesystem",
]