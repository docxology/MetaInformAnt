"""Common utilities for simulation scripts.

This module provides shared functionality for all simulation scripts including
parameter validation, progress reporting, and error handling.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)


def validate_positive_int(value: int, name: str, min_val: int = 1) -> None:
    """Validate that an integer is positive.

    Args:
        value: Value to validate
        name: Parameter name for error messages
        min_val: Minimum allowed value

    Raises:
        ValidationError: If value is invalid
    """
    validation.validate_type(value, int, name)
    validation.validate_range(value, min_val=min_val, name=name)


def validate_float_range(value: float, name: str, min_val: float = 0.0, max_val: float = 1.0) -> None:
    """Validate that a float is in a specified range.

    Args:
        value: Value to validate
        name: Parameter name for error messages
        min_val: Minimum allowed value
        max_val: Maximum allowed value

    Raises:
        ValidationError: If value is invalid
    """
    validation.validate_type(value, float, name)
    validation.validate_range(value, min_val=min_val, max_val=max_val, name=name)


def validate_output_dir(output: Path) -> Path:
    """Validate and ensure output directory exists.

    Args:
        output: Output directory path

    Returns:
        Validated and created output directory

    Raises:
        ValidationError: If path is invalid
    """
    validation.validate_type(output, (Path, str), "output")
    output_path = Path(output) if isinstance(output, str) else output
    return paths.ensure_directory(output_path)


def log_simulation_start(simulation_type: str, **params: Any) -> None:
    """Log simulation start with parameters.

    Args:
        simulation_type: Type of simulation
        **params: Simulation parameters to log
    """
    logger.info(f"Starting {simulation_type} simulation")
    for key, value in params.items():
        logger.debug(f"  {key}: {value}")


def log_simulation_complete(output_file: Path | str) -> None:
    """Log simulation completion.

    Args:
        output_file: Path to output file
    """
    logger.info(f"Simulation complete. Output: {output_file}")


def handle_simulation_error(error: Exception, context: str = "") -> None:
    """Handle simulation errors with context.

    Args:
        error: Exception that occurred
        context: Additional context string
    """
    context_msg = f" ({context})" if context else ""
    logger.error(f"Simulation failed{context_msg}: {error}", exc_info=True)
    raise
