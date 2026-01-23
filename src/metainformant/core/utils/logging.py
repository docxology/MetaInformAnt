"""Logging utilities for METAINFORMANT.

Provides consistent, structured logging across all modules with support for
console and file output, environment-based configuration, and metadata logging.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path


def get_logger(name: str) -> logging.Logger:
    """Get or create a logger with default console handler.

    Args:
        name: Logger name (typically __name__)

    Returns:
        Configured logger with console handler if none exists
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


def setup_logger(name: str, log_file: str | None = None, level: str = "INFO") -> logging.Logger:
    """Set up a logger with file and/or console output.

    Args:
        name: Logger name
        log_file: Optional log file path
        level: Logging level

    Returns:
        Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))

    # Clear existing handlers
    logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler if specified
    if log_file:
        from pathlib import Path

        Path(log_file).parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger_with_level(name: str, level: str | int | None = None) -> logging.Logger:
    """Get or create a logger with specified log level.

    Args:
        name: Logger name (typically __name__)
        level: Logging level as string ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
            or integer. If None, uses INFO or CORE_LOG_LEVEL environment variable.

    Returns:
        Configured logger with console handler if none exists

    Examples:
        >>> logger = get_logger_with_level("my.module", "DEBUG")
        >>> logger.debug("This will be shown")
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Set level
    if level is None:
        # Check environment variable
        env_level = os.environ.get("CORE_LOG_LEVEL", "INFO")
        level = env_level
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)
    logger.setLevel(level)
    return logger


def configure_logging_from_env(default_level: str = "INFO") -> None:
    """Configure root logger from environment variables.

    Reads CORE_LOG_LEVEL environment variable and configures the root logger
    accordingly. This affects all loggers that inherit from root.

    Args:
        default_level: Default log level if CORE_LOG_LEVEL is not set

    Examples:
        >>> configure_logging_from_env()
        >>> # Now all loggers will use CORE_LOG_LEVEL if set
    """
    level_str = os.environ.get("CORE_LOG_LEVEL", default_level)
    level = getattr(logging, level_str.upper(), logging.INFO)
    logging.root.setLevel(level)


def log_with_metadata(
    logger: logging.Logger, message: str, metadata: dict, *, level: str = "INFO", structured: bool = False
) -> None:
    """Log message with structured metadata.

    Args:
        logger: Logger instance
        message: Log message
        metadata: Dictionary of metadata to include
        level: Log level to use ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
        structured: If True, writes metadata as separate JSON log entry for
            structured logging systems. If False, appends JSON to message.

    Examples:
        >>> logger = get_logger(__name__)
        >>> log_with_metadata(logger, "Processing complete", {
        ...     "samples": 100, "time_seconds": 45.2
        ... })
        >>> # With structured logging
        >>> log_with_metadata(logger, "Batch processed", {
        ...     "batch_id": "batch_001", "records": 1000
        ... }, structured=True)
    """
    log_level = getattr(logging, level.upper(), logging.INFO)

    if structured:
        # For structured logging, log metadata as separate JSON entry
        metadata_str = json.dumps(metadata, default=str)
        logger.log(log_level, f"{message}")
        logger.log(log_level, f"METADATA: {metadata_str}")
    else:
        # Format metadata as JSON appended to message
        metadata_str = json.dumps(metadata, default=str)
        full_message = f"{message} | {metadata_str}"
        logger.log(log_level, full_message)
