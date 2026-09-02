# MCP Package Specification

## Module

`metainformant.mcp`

## Status

Implemented: stdio JSON-RPC 2.0 MCP server with a validated declarative tool
registry. Standard library only (no MCP SDK dependency).

## Implemented Interface

### Server (`metainformant.mcp.server`)

Run:

```bash
uv run python -m metainformant.mcp.server
```

- Transport: newline-delimited JSON-RPC 2.0 over stdin/stdout; diagnostics on
  stderr only.
- Methods: `initialize`, `tools/list`, `tools/call`, `resources/list`,
  `resources/read`; `exit` notification shuts down cleanly (exit 0).
- Protocol version: `2024-11-05` (echoed back to the client on initialize).
- Error semantics: parse failures `-32700`; non-object requests `-32600`;
  unknown methods `-32601`; bad params / schema violations / unknown tool or
  resource `-32602`; handler failures `-32603`.
- Notifications (no `id`) produce no response.
- Nothing that is not a JSON-RPC response is ever written to stdout.

### Registry (`metainformant.mcp.registry`)

- `Tool(name, description, input_schema, handler, metadata)` — frozen record.
- `ToolRegistry.register` validates name, description, schema shape, and
  handler callability at registration time; duplicates rejected.
- `ToolRegistry.register_module_tools(module)` imports a tool module and
  registers its declared `TOOL` / `TOOLS` constants; missing modules or empty
  declarations raise loudly.
- `ToolRegistry.call(name, arguments)` validates arguments against the tool
  schema before invoking the handler (`SchemaError` on violation).
- Supported schema subset: `type` (object/string/integer/number/boolean/array),
  `required`, `properties`, `enum`, `additionalProperties: false`. `bool` is
  not accepted where `integer`/`number` is declared.

### Bundled tool adapter (`metainformant.mcp.tool_adapters`)

- Registers `amalgkit_monitor` by importing (never modifying)
  `metainformant.mcp.tools.amalgkit_monitor.build_status`.
- Arguments: `data_root` (string, default `$AMALGKIT_DATA_ROOT` or
  `output/amalgkit`), `log_file` (optional string), `inspect_processes`
  (boolean, default false).
- Output is the monitor's descriptive operational snapshot; biological
  inference is always `withheld`.

## Package Exports

```python
from metainformant.mcp import TOOLS_MODULES
from metainformant.mcp.registry import Tool, ToolRegistry, SchemaError
from metainformant.mcp.server import MCPServer, build_default_registry
```

`TOOLS_MODULES` is the declarative tuple of tool-adapter modules wired into the
default registry. Adding a tool = write an adapter module declaring `TOOL`/
`TOOLS`, append its import path to `TOOLS_MODULES`.

## Resources

- `metainformant://capabilities` — JSON index of server info, protocol version,
  methods, and live tool descriptors.
- `metainformant://methods` — Markdown connection/usage documentation.

## Client Connection (stdio)

Any MCP stdio client launches:

```
uv run python -m metainformant.mcp.server
```

from the repository root (or with `src/` on `PYTHONPATH`), then speaks
newline-delimited JSON-RPC 2.0 on the child's stdin/stdout. Example handshake
line:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05"}}
```

## Tests

`tests/mcp/test_registry.py` (schema/registration/validation) and
`tests/mcp/test_server.py` (handshake, tools/list completeness, tools/call
round-trip, schema-violation rejection, resources, and a real subprocess
stdio round-trip). Zero mocks; run with:

```bash
env -u VIRTUAL_ENV uv run pytest -q tests/mcp/test_registry.py tests/mcp/test_server.py
```

## Non-Goals For This Pass

- no SSE/HTTP transport;
- no prompts/list or sampling surfaces;
- no `run_workflow` or `list_outputs` tools (see `tools/` subpackage).
