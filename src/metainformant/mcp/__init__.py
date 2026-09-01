"""MCP (Model Context Protocol) helpers for METAINFORMANT.

The current public surface is a standalone amalgkit monitor module under
``metainformant.mcp.tools``.  A full MCP server is not implemented yet.
"""

from __future__ import annotations

__all__ = ["tools"]


def __getattr__(name):
    """Lazily expose the ``tools`` subpackage to honor ``__all__``."""
    if name in __all__:
        import importlib

        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
