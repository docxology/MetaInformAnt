# MCP Package Specification

## Module

`metainformant.mcp`

## Status

Implemented. `metainformant.mcp.server` provides a stdio JSON-RPC 2.0 MCP
server with a schema-validated tool registry; tests live under `tests/mcp/`.

## Implemented Interface

### `metainformant.mcp.server`

- Entry point: `python -m metainformant.mcp.server` (also exposed as
  `serve_stdio(...)` and `main()`).
- JSON-RPC methods: `initialize`, `tools/list`, `tools/call`,
  `resources/list`, `resources/read`; shutdown via EOF or the `exit`
  notification.
- Read-only resources: `metainformant://capabilities` (registered tools and
  schemas) and `metainformant://methods` (connection documentation).

### `metainformant.mcp.registry`

- `Tool` (name, description, JSON-schema `input_schema`, handler),
  `ToolRegistry.register` / `register_module_tools` / `get` / `names`, and
  `SchemaError`.
- Schemas are validated at registration; arguments at call time.

### `metainformant.mcp.tool_adapters`

- Bundled tool adapters wired into the default registry
  (`build_default_registry()`), alongside `metainformant.mcp.tools.catalog`
  which adapts the 20 `TOOL_SPEC` tools from `metainformant.mcp.tools`.

### `metainformant.mcp.tools.amalgkit_monitor`

Command:

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor \
  --data-root /path/to/campaign --no-process-scan
```

Responsibilities:

- inspect local Amalgkit/RNA pipeline processes when explicitly requested;
- parse an explicitly supplied log file for operational diagnostics
  (`parse_log_progress(log_file)` returns processed/total/percent/last_line);
- read the campaign `pipeline_progress.db` and provenance receipts;
- return separate executable, cohort, descriptive-analysis, and biological-
  inference readiness fields suitable for terminal use and adapters.

The default data root is `AMALGKIT_DATA_ROOT` or `output/amalgkit`. Use
`--no-process-scan` for a receipt/lock/database-only snapshot.

The helper never reports biological inference as ready. A current descriptive
receipt is not a completed scientific result.

## Package Exports

```python
from metainformant.mcp import tools
from metainformant.mcp.registry import Tool, ToolRegistry, SchemaError
from metainformant.mcp.server import MCPServer, build_default_registry
from metainformant.mcp.tools import amalgkit_monitor
```

## Adding a Tool

1. Write an adapter module declaring a `Tool` as `TOOL` or `TOOLS`.
2. Append the module import path to `TOOLS_MODULES` in
   `metainformant/mcp/__init__.py`.
3. Add tests; the registry validates schemas at registration and arguments at
   call time.
