"""Test orchestrator for METAINFORMANT.

Provides programmatic entry points to run the pytest suite from the package CLI.
"""

from __future__ import annotations

from .runner import run_all_tests

__all__ = ["run_all_tests"]
