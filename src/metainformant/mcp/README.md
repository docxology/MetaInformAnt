# MCP Package

Model Context Protocol server for METAINFORMANT: a stdio JSON-RPC 2.0 server
with a declarative, schema-validated tool registry. Standard library only.

## Public Surface

```bash
# Run the MCP server over stdio
uv run python -m metainformant.mcp.server
```

```python
from metainformant.mcp.registry import Tool, ToolRegistry, SchemaError
from metainformant.mcp.server import MCPServer, build_default_registry
from metainformant.mcp.tools import amalgkit_monitor
```

## Available Modules

| Module | Purpose | Status |
|--------|---------|--------|
| `server` | stdio JSON-RPC 2.0 MCP server (initialize, tools/list, tools/call, resources) | Implemented |
| `registry` | Declarative `Tool`/`ToolRegistry` with schema validation | Implemented |
| `tool_adapters` | Bundled tool adapters wired into the default registry | Implemented |
| `tools.amalgkit_monitor` | Inspect local Amalgkit/RNA workflow progress | Implemented |

## Quick Start (client)

Launch the server as a subprocess from the repository root and speak
newline-delimited JSON-RPC 2.0 on its stdin/stdout:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05"}}
{"jsonrpc":"2.0","id":2,"method":"tools/list"}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"amalgkit_monitor","arguments":{}}}
```

Hermes/Claude-style clients register the same stdio command:
`uv run python -m metainformant.mcp.server` (working directory: repository
root, so the default data root `output/amalgkit` resolves, or pass
`data_root` per call).

## Adding a Tool

1. Write an adapter module (e.g. `metainformant/mcp/tool_adapters.py` pattern)
   declaring a `Tool` (name, description, JSON-schema `input_schema`, handler)
   as `TOOL` or `TOOLS`.
2. Append the module import path to `TOOLS_MODULES` in
   `metainformant/mcp/__init__.py`.
3. Add tests; the registry validates schemas at registration and arguments at
   call time.

## Related

- [SPEC.md](SPEC.md) — full protocol and registry contract
- [docs/mcp/](../../../docs/mcp/)
- [RNA package](../rna/)
