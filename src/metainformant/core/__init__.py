"""Core utilities for METAINFORMANT bioinformatics toolkit."""
from __future__ import annotations

# isort: skip_file

from . import data, engine, execution, io, ui, utils

# Optional dependencies
try:
    from .data import db
except ImportError:
    db = None

__all__ = ["data", "db", "engine", "execution", "io", "ui", "utils"]
