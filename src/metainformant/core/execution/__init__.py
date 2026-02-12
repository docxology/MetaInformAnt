"""Core execution utilities for METAINFORMANT.

This module provides workflow orchestration, parallel processing, and code discovery utilities."""
from __future__ import annotations

from . import discovery, parallel, workflow

__all__ = ['discovery', 'parallel', 'workflow']
