"""MCP (Model Context Protocol) package for METAINFORMANT.

Provides a stdio JSON-RPC 2.0 MCP server (:mod:`metainformant.mcp.server`) and
a declarative tool registry (:mod:`metainformant.mcp.registry`).  Tools are
bundled tool adapters under ``metainformant.mcp.tools`` plus this package's
``tool_adapters`` module; every tool is registered through
``metainformant.mcp.registry.ToolRegistry``.
"""

from __future__ import annotations

#: Declarative list of tool modules registered into the default MCP registry.
TOOLS_MODULES = ("metainformant.mcp.tool_adapters", "metainformant.mcp.tools.catalog")

__all__ = ["TOOLS_MODULES", "registry", "server", "tool_adapters", "tools"]


def __getattr__(name: str):
    """Lazily expose submodules to honor ``__all__`` without import side effects."""

    if name in __all__:
        import importlib

        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
