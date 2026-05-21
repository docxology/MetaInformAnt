"""Standalone MCP-adjacent command modules."""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import TYPE_CHECKING

__all__ = ["amalgkit_monitor"]


def __getattr__(name: str) -> ModuleType:
    """Lazily import standalone tool modules."""
    if name in __all__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


if TYPE_CHECKING:
    from . import amalgkit_monitor
