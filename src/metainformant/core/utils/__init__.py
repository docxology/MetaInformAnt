"""Core utilities for METAINFORMANT bioinformatics toolkit.

This module provides foundational utility functions used across all domain modules.
Import utilities directly from this package for convenient access.

Example:
    from metainformant.core.utils import get_logger, load_mapping_from_file

    logger = get_logger(__name__)
    config = load_mapping_from_file("settings.yaml")
"""

from __future__ import annotations

# Configuration Management
from .config import (
    PostgresConfig,
    load_postgres_config_from_env,
    load_config_file,
    get_env_or_default,
    load_typed_env,
    load_mapping_from_file,
    merge_configs,
    coerce_config_types,
    apply_env_overrides,
    discover_config_files,
    get_config_schema,
    find_configs_for_module,
    list_config_templates,
)

# Error Handling
from .errors import (
    METAINFORMANTError,
    ConfigError,
    IOError,
    ValidationError,
    NetworkError,
    CacheError,
    DownloadError,
    PipelineError,
    DependencyError,
    ResourceError,
    retry_with_backoff,
    error_context,
    safe_execute,
    validate_not_none,
    validate_type,
)

# Hashing Utilities
from .hash import (
    sha256_bytes,
    sha256_file,
    sha256_string,
    deterministic_seed,
    file_hash_comparison,
    hash_directory,
    verify_file_integrity,
)

# Logging Framework
from .logging import (
    get_logger,
    setup_logger,
    get_logger_with_level,
    configure_logging_from_env,
    log_with_metadata,
)

# Optional Dependencies
from .optional_deps import (
    suppress_optional_warnings,
    enable_optional_warnings,
    warn_optional_dependency,
    reset_warning_state,
    get_warning_state,
)

# Progress Tracking
from .progress import (
    progress_bar,
    task_context,
    log_progress,
)

# Text Processing
from .text import (
    normalize_whitespace,
    slugify,
    safe_filename,
    clean_whitespace,
    remove_control_chars,
    standardize_gene_name,
    format_species_name,
    clean_sequence_id,
    extract_numbers,
    truncate_text,
)

# Timing & Performance
from .timing import (
    Timer,
    rate_limiter,
    timed,
    timeout_after,
)

# Symbol Indexing (lazy import to avoid circular deps)
# from .symbols import index_functions, find_symbol, etc.

__all__ = [
    # Configuration
    "PostgresConfig",
    "load_postgres_config_from_env",
    "load_config_file",
    "get_env_or_default",
    "load_typed_env",
    "load_mapping_from_file",
    "merge_configs",
    "coerce_config_types",
    "apply_env_overrides",
    "discover_config_files",
    "get_config_schema",
    "find_configs_for_module",
    "list_config_templates",
    # Error Handling
    "METAINFORMANTError",
    "ConfigError",
    "IOError",
    "ValidationError",
    "NetworkError",
    "CacheError",
    "DownloadError",
    "PipelineError",
    "DependencyError",
    "ResourceError",
    "retry_with_backoff",
    "error_context",
    "safe_execute",
    "validate_not_none",
    "validate_type",
    # Hashing
    "sha256_bytes",
    "sha256_file",
    "sha256_string",
    "deterministic_seed",
    "file_hash_comparison",
    "hash_directory",
    "verify_file_integrity",
    # Logging
    "get_logger",
    "setup_logger",
    "get_logger_with_level",
    "configure_logging_from_env",
    "log_with_metadata",
    # Optional Dependencies
    "suppress_optional_warnings",
    "enable_optional_warnings",
    "warn_optional_dependency",
    "reset_warning_state",
    "get_warning_state",
    # Progress
    "progress_bar",
    "task_context",
    "log_progress",
    # Text Processing
    "normalize_whitespace",
    "slugify",
    "safe_filename",
    "clean_whitespace",
    "remove_control_chars",
    "standardize_gene_name",
    "format_species_name",
    "clean_sequence_id",
    "extract_numbers",
    "truncate_text",
    # Timing & Performance
    "Timer",
    "rate_limiter",
    "timed",
    "timeout_after",
]
