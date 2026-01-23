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
    validate_type,
    validate_range,
    validate_path_exists,
    validate_not_none,
    validate_not_empty,
    validate_schema,
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
    load_json,
    dump_json,
    load_yaml,
    download_file,
)

# 4. Execution (Workflow)
# Depends on IO and Utils
from . import execution
from .execution import (
    discovery,
    parallel,
    workflow,
    discover_functions,
    discover_configs,
    thread_map,
    parallel_batch,
    cpu_count,
)

# 5. Engine (Workflow Manager)
from . import engine
from .engine import (
    WorkflowManager,
    SampleStage,
    SampleState,
)

# 6. UI (Terminal Interface)
from . import ui
from .ui import (
    ProgressState,
    TerminalInterface,
)

# Optional dependencies
try:
    from .data import db
except ImportError:
    db = None

# Core functions - always available
from .execution.workflow import (
    create_sample_config,
    download_and_process_data,
    run_config_based_workflow,
    validate_config_file,
)

# Type system
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

__all__ = [
    # Subpackages
    "data",
    "engine",
    "execution",
    "io",
    "ui",
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
    # Key utility functions
    "validate_type",
    "validate_range",
    "validate_path_exists",
    "validate_not_none",
    "validate_not_empty",
    "validate_schema",
    "load_json",
    "dump_json",
    "load_yaml",
    "download_file",
    "discover_functions",
    "discover_configs",
    "thread_map",
    "parallel_batch",
    "cpu_count",
    # Engine
    "WorkflowManager",
    "SampleStage",
    "SampleState",
    # UI
    "ProgressState",
    "TerminalInterface",
    # Optional modules
    "db",
    "disk",
]
