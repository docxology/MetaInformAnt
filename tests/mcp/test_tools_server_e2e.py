"""End-to-end stdio server test with the catalog + adapter tools registered."""

from __future__ import annotations

import io
import json

from metainformant.mcp.registry import ToolRegistry
from metainformant.mcp.server import MCPServer
from metainformant.mcp.tools.catalog import build_registry_tools


def _build_server() -> MCPServer:
    registry = ToolRegistry()
    for tool in build_registry_tools():
        registry.register(tool)
    registry.register_module_tools("metainformant.mcp.tool_adapters")
    return MCPServer(registry)


def _serve(messages: list[dict]) -> list[dict]:
    server = _build_server()
    out = io.StringIO()
    payload = "\n".join(json.dumps(m) for m in messages) + "\n"
    server.serve(stdin=io.StringIO(payload), stdout=out, stderr=io.StringIO())
    return [json.loads(line) for line in out.getvalue().strip().splitlines()]


def test_stdio_tools_lists_merged_catalog() -> None:
    responses = _serve([{"jsonrpc": "2.0", "id": 1, "method": "tools/list", "params": {}}])
    names = [t["name"] for t in responses[0]["result"]["tools"]]
    assert "amalgkit_monitor" in names  # registry lane's adapter
    assert "dna_translate" in names  # this lane's catalog
    assert len(names) >= 15


def test_stdio_tools_call_dna_translate() -> None:
    responses = _serve(
        [
            {
                "jsonrpc": "2.0",
                "id": 3,
                "method": "tools/call",
                "params": {
                    "name": "dna_translate",
                    "arguments": {"sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"},
                },
            }
        ]
    )
    result = responses[0]["result"]
    assert result["isError"] is False
    assert result["structuredContent"]["protein"] == "MAIVMGR*KGAR*"


def test_stdio_unknown_tool_returns_error_response() -> None:
    responses = _serve(
        [
            {
                "jsonrpc": "2.0",
                "id": 9,
                "method": "tools/call",
                "params": {"name": "no_such_tool", "arguments": {}},
            }
        ]
    )
    # Unknown tools surface as a JSON-RPC error object, not a result envelope.
    assert "error" in responses[0]
