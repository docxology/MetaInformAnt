"""Centralized optional dependency warning management.

This module provides centralized management for optional dependency warnings,
preventing duplicate warnings and allowing suppression of warnings during
environment setup (e.g., before venv activation).
"""

from __future__ import annotations

import logging
from typing import Set

# Module-level state for warning management
_optional_warnings_issued: Set[str] = set()
_suppress_warnings: bool = False


def suppress_optional_warnings() -> None:
    """Suppress optional dependency warnings (e.g., before venv activation).

    This is useful during environment setup when warnings about missing
    optional dependencies are not helpful to the user.
    """
    global _suppress_warnings
    _suppress_warnings = True


def enable_optional_warnings() -> None:
    """Enable optional dependency warnings.

    This should be called after environment setup is complete to allow
    warnings about missing optional dependencies to be displayed.
    """
    global _suppress_warnings
    _suppress_warnings = False


def warn_optional_dependency(module_name: str, functionality: str, fallback: str = "functionality disabled") -> None:
    """Warn about missing optional dependency (only once, and only if enabled).

    Args:
        module_name: Name of the missing module (e.g., "seaborn", "anndata")
        functionality: Description of affected functionality (e.g., "enhanced plots")
        fallback: Description of fallback behavior (e.g., "basic plots used")

    This function ensures:
    - Warnings are only issued once per session per module
    - Warnings can be suppressed during environment setup
    - Consistent warning format across all modules
    """
    if _suppress_warnings:
        return

    warning_key = f"{module_name}:{functionality}"
    if warning_key in _optional_warnings_issued:
        return

    # Use core logging to avoid circular imports
    logger = logging.getLogger("metainformant.core.optional_deps")
    logger.warning(f"{module_name} not available, {functionality} {fallback}")
    _optional_warnings_issued.add(warning_key)


def reset_warning_state() -> None:
    """Reset warning state (useful for testing)."""
    global _optional_warnings_issued, _suppress_warnings
    _optional_warnings_issued.clear()
    _suppress_warnings = False


def get_warning_state() -> dict[str, bool | Set[str]]:
    """Get current warning state (useful for debugging)."""
    return {"suppress_warnings": _suppress_warnings, "warnings_issued": _optional_warnings_issued.copy()}
