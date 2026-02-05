"""Core utilities for METAINFORMANT bioinformatics toolkit.

This module provides essential infrastructure components for the METAINFORMANT
bioinformatics analysis platform, including configuration management, I/O operations,
parallel processing, caching, and workflow orchestration.
"""

from __future__ import annotations

# 6. UI (Terminal Interface)
# 5. Engine (Workflow Manager)
# 4. Execution (Workflow)
# Depends on IO and Utils
# 3. IO (Input/Output)
# Depends on Utils for logging/config
# 2. Data (Structures)
# Import subpackages (reordered for dependency resolution)
# 1. Utils (Foundation)
from . import data, engine, execution, io, ui, utils
from .data import (
    validate_not_empty,
    validate_not_none,
    validate_path_exists,
    validate_range,
    validate_schema,
    validate_type,
    validation,
)
from .engine import (
    SampleStage,
    SampleState,
    WorkflowManager,
)
from .execution import (
    cpu_count,
    discover_configs,
    discover_functions,
    discovery,
    parallel,
    parallel_batch,
    thread_map,
    workflow,
)
from .io import (
    atomic_replace,
    atomic_write,
    cache,
    compute_checksums_batch,
    compute_md5,
    compute_sha256,
    disk,
    download,
    download_file,
    dump_json,
)
from .io import io as io_module
from .io import (
    load_json,
    load_yaml,
    paths,
    safe_write_bytes,
    safe_write_text,
    temp_directory,
    verify_checksum,
    verify_checksum_file,
    write_checksum_file,
)
from .ui import (
    ProgressState,
    TerminalInterface,
)
from .utils import (
    Timer,
    config,
    errors,
    hash,
    logging,
    optional_deps,
    progress,
    rate_limiter,
    symbols,
    text,
    timed,
    timeout_after,
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
    # Timing & Performance
    "Timer",
    "rate_limiter",
    "timed",
    "timeout_after",
    # Atomic IO
    "atomic_replace",
    "atomic_write",
    "safe_write_bytes",
    "safe_write_text",
    "temp_directory",
    # Checksums
    "compute_checksums_batch",
    "compute_md5",
    "compute_sha256",
    "verify_checksum",
    "verify_checksum_file",
    "write_checksum_file",
]
