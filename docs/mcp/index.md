# MCP

METAINFORMANT ships a Model Context Protocol server: a stdio JSON-RPC 2.0
server with a declarative, schema-validated tool registry (standard library
only).

## Run the server

```bash
uv run python -m metainformant.mcp.server
```

The server speaks newline-delimited JSON-RPC 2.0 on stdin/stdout and
implements `initialize`, `tools/list`, `tools/call`, `resources/list`, and
`resources/read`, plus graceful shutdown via EOF or the `exit` notification.

Client example (from the repository root, so the default data root
`output/amalgkit` resolves):

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05"}}
{"jsonrpc":"2.0","id":2,"method":"tools/list"}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"amalgkit_monitor","arguments":{}}}
```

Read-only resources: `metainformant://capabilities` (registered tools and
schemas) and `metainformant://methods` (connection documentation).

## Python API

```python
from metainformant.mcp.registry import Tool, ToolRegistry, SchemaError
from metainformant.mcp.server import MCPServer, build_default_registry
from metainformant.mcp.tools import amalgkit_monitor
```

The default registry bundles the `amalgkit_monitor` adapter plus the 20
`TOOL_SPEC` tools from `metainformant.mcp.tools` (core, dna, gwas, math,
protein, rna, and visualization surfaces).

## Amalgkit monitor

The monitor also runs standalone for terminal use:

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor \
  --data-root "$AMALGKIT_DATA_ROOT" --no-process-scan
```

From Python, parse a workflow log file directly (`parse_log_progress` expects a
log file path, not a directory):

```python
from pathlib import Path
from metainformant.mcp.tools import amalgkit_monitor

progress = amalgkit_monitor.parse_log_progress(
    Path("output/amalgkit/run_all_species_incremental.log")
)
# {"processed": int, "total": int, "percent": float, "last_line": str}
# A missing file returns zeros with "No log file found".
```

The monitor is not a completion oracle: it separates executable readiness,
cohort readiness, descriptive analysis, and biological inference, and always
withholds the last field.

## Adding a tool

1. Write an adapter module (see `metainformant/mcp/tool_adapters.py`) declaring
   a `Tool` (name, description, JSON-schema `input_schema`, handler) as `TOOL`
   or `TOOLS`.
2. Append the module import path to `TOOLS_MODULES` in
   `metainformant/mcp/__init__.py`.
3. Add tests; the registry validates schemas at registration and arguments at
   call time.

## Related

- [MCP Notes](README.md)
- [SPEC](SPEC.md)
- [Source README](../../src/metainformant/mcp/README.md)
- [RNA workflow docs](../rna/)
