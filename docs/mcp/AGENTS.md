# Agent Directives: mcp (Model Context Protocol)

## Role

Expose METAINFORMANT toolset as Model Context Protocol tools so LLM-based development environments (Claude Desktop, Cursor, Windsurf) can invoke bioinformatics workflows directly.

Think of this as a **translation layer**: LLM-friendly JSON-RPC → internal Python function calls.

## Module Scope

| File/Dir | Purpose |
|-----------|---------|
| `server.py` | JSON-RPC 2.0 server implementation (stdio + SSE transports) |
| `tools/` | Individual tool implementations using `@register_tool` decorator |
| `__init__.py` | Tool registry, MCP tool spec dataclasses |
| `transport.py` | STDIO and SSE transport mechanisms |

## Key Source Files

- **Server**: `src/metainformant/mcp/server.py` — JSON-RPC request router, tool dispatcher, error handling
- **Decorator**: `src/metainformant/mcp/__init__.py` — `@register_tool` and `ToolSpec` dataclass
- **Tools**: `src/metainformant/mcp/tools/amalgkit_monitor.py` — First-endemic tool (pipeline status)
- **Transport**: `src/metainformant/mcp/transport.py` — STDIO pipe and HTTP+SSE implementations

## Related Documentation

- **Module guide**: [index.md](index.md) — LLM integration examples, conversation flows, setup instructions
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **API spec**: [SPEC.md](SPEC.md) — JSON-RPC methods, transport protocols, error model
- **MCP spec**: https://spec.modelcontextprotocol.io — Official protocol documentation
- **Claude Desktop**: https://docs.anthropic.com/claude/docs/claude-desktop-mcp — Installation guide
- **Cursor integration**: https://docs.cursor.com/context/model-context-protocol — Configuration
- **Amalgkit**: [../rna/amalgkit/](../rna/amalgkit/) — The pipeline that MCP tools monitor

## Rules & Constraints

1. Tools MUST be stateless — no modification of global server state
2. Tools MUST be idempotent — same args → same result
3. Tool responses MUST be JSON-serializable (dataclasses → `asdict()`)
4. Timeouts enforced per-tool (default 300s, configurable via `METAINFORMANT_MCP_TIMEOUT`)
5. Never leave stdout/stderr pipes open (breaks Claude Desktop)

## Cross-Module Dependencies

- **rna**: [`metainformant.rna.engine.monitoring`](../../src/metainformant/rna/engine/monitoring.py) — `get_pipeline_status()` hooked by `amalgkit_status` tool
- **core.io**: [`metainformant.core.io`](../../src/metainformant/core/io/) — File operations (list_outputs, read_results)
- **visualization**: [`metainformant.visualization`](../../src/metainformant/visualization/) — Future plot-generation tools

## Development Notes

Current implementation (v0.3.0): STDIO transport + 1 tool. Plan: SSE transport, 8+ tools across rna/gwas/multiomics modules.
