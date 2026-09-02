"""Catalog of module-spec MCP tools and their registry bridges.

Two layers live here:

1. ``iter_tool_specs()`` / ``MODULE_TOOLS_SPECS``: the raw per-module
   ``TOOL_SPEC`` dictionaries (name, description, JSON-schema input, plain
   ``**kwargs`` handler, writes class). Kept independent of the MCP registry
   contract so both lanes' tools can be enumerated uniformly.

2. ``build_registry_tools()``: adapts every spec into a
   :class:`metainformant.mcp.registry.Tool` (single-Mapping handler), so the
   server can register the whole catalog via one module:

   .. code-block:: python

       registry.register_module_tools("metainformant.mcp.tools.catalog")

Handlers are invoked through the adapters in ``*_tools`` modules unchanged;
nothing in the underlying analysis code is modified.
"""

from __future__ import annotations

import importlib
from collections.abc import Iterator, Mapping
from typing import Any

from metainformant.mcp.registry import Tool

#: Tool-spec modules shipped with the toolkit (amalgkit_monitor is wired by the
#: registry lane's tool_adapters.py and is NOT duplicated here).
TOOL_SPEC_MODULES: tuple[str, ...] = (
    "metainformant.mcp.tools.core_tools",
    "metainformant.mcp.tools.dna_tools",
    "metainformant.mcp.tools.gwas_tools",
    "metainformant.mcp.tools.math_tools",
    "metainformant.mcp.tools.protein_tools",
    "metainformant.mcp.tools.rna_tools",
    "metainformant.mcp.tools.visualization_tools",
)


def iter_tool_specs() -> Iterator[dict[str, Any]]:
    """Yield every registered TOOL_SPEC across the catalog modules."""
    for module_name in TOOL_SPEC_MODULES:
        module = importlib.import_module(module_name)
        specs = getattr(module, "ALL_SPECS", None)
        if specs is None:
            spec = getattr(module, "TOOL_SPEC", None)
            if spec is not None:
                specs = [spec]
        yield from specs or []


def _adapt_handler(handler: Any) -> Any:
    """Wrap a **kwargs handler into the registry's Mapping-argument form."""

    def _call(arguments: Mapping[str, Any]) -> dict[str, Any]:
        return handler(**arguments)

    return _call


def build_registry_tools() -> list[Tool]:
    """Adapt every catalog spec into a validated registry Tool record."""
    tools: list[Tool] = []
    for spec in iter_tool_specs():
        tools.append(
            Tool(
                name=spec["name"],
                description=spec["description"],
                input_schema=spec["input_schema"],
                handler=_adapt_handler(spec["handler"]),
                metadata={"writes": spec.get("writes", "read-only")},
            )
        )
    return tools


TOOLS = build_registry_tools()
