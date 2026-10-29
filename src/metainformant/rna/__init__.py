"""RNA transcriptomic analysis and workflow orchestration.

Subpackages load on first access so acquisition workers do not initialize
analysis plotting or platform font discovery during startup.
"""
from __future__ import annotations

import importlib
from types import ModuleType

__all__ = ["amalgkit", "analysis", "core", "deconvolution", "engine", "retrieval", "splicing", "steps"]


def __getattr__(name: str) -> ModuleType:
    """Load the existing public submodule only when it is requested."""
    if name in __all__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Keep discoverability of public modules before their first import."""
    return sorted(set(globals()) | set(__all__))
