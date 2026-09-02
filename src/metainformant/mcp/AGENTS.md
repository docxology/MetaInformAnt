# Agent Directives: mcp

**Context**: MCP server package for METAINFORMANT — stdio JSON-RPC 2.0
transport, declarative tool registry, bundled tool adapters.

## Capabilities

- `server.py` — newline-delimited JSON-RPC 2.0 over stdin/stdout: `initialize`,
  `tools/list`, `tools/call`, `resources/list`, `resources/read`; `exit`
  notification for clean shutdown. Stdlib only; deterministic error codes.
- `registry.py` — `Tool` records + `ToolRegistry`: schema validation at
  registration (`SchemaError` on malformed schemas), argument validation at
  dispatch, loud failure on missing tool modules.
- `tool_adapters.py` — registers `amalgkit_monitor` by importing
  `metainformant.mcp.tools.amalgkit_monitor` (never modifying it).

## Subpackages

### tools/

| File | Purpose |
|------|---------|
| `amalgkit_monitor.py` | Amalgkit pipeline status monitoring (owned by the mcp-tools lane; import, do not edit here) |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Stdio server writes JSON-RPC responses ONLY to stdout; diagnostics go to stderr
- Tools must be stateless and idempotent; handlers return JSON-serializable values
- New tools go through an adapter module declaring `TOOL`/`TOOLS`, appended to
  `TOOLS_MODULES` in `__init__.py`; the registry rejects undeclared or invalid tools
- Descriptive-only boundary: tool outputs must not make inferential
  cross-species claims; the amalgkit_monitor adapter always reports
  `biological_inference: withheld`
- Follow REAL IMPLEMENTATION policy — all tests use real processes/data (the
  server round-trip test spawns the real `python -m metainformant.mcp.server`)
- Use `uv` for dependency management
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/mcp/test_registry.py tests/mcp/test_server.py`
